#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using Libdl
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TRACE_DIR = joinpath(PROJECT_ROOT, "build", "trace")

@inline function _resolve_sysimage_path()
    if !isempty(ARGS)
        return normpath(ARGS[1])
    end
    if haskey(ENV, "JRT_SYSIMAGE_PATH") && !isempty(ENV["JRT_SYSIMAGE_PATH"])
        return normpath(ENV["JRT_SYSIMAGE_PATH"])
    end
    return joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime.$(Libdl.dlext)")
end

@inline function _resolve_label()
    if length(ARGS) >= 2 && !isempty(ARGS[2])
        return String(ARGS[2])
    end
    return "default"
end

const WORKLOAD_EXPR = """
if !isdefined(Main, :Models)
    Base.include(Main, joinpath("__PROJECT_ROOT__", "src", "models", "Models.jl"))
end
using .Models
using ForwardDiff
using Main.Constants_PNJL: ħc_MeV_fm
model = Models.create_model(:PNJL)
T_fm = 150.0 / ħc_MeV_fm
mu_fm = 250.0 / ħc_MeV_fm
residual! = Models.Conditions.build_residual!(
    model,
    Models.FixedMu(),
    T_fm,
    mu_fm;
    xi=0.0,
    p_num=6,
    t_num=3,
)
x0 = Float64.(Models.HADRON_SEED_5)
vec_residual = x -> begin
    F = similar(x)
    residual!(F, x)
    return F
end
ForwardDiff.jacobian(vec_residual, x0)
"""

function _workload_expr()
    return replace(WORKLOAD_EXPR, "__PROJECT_ROOT__" => replace(PROJECT_ROOT, "\\" => "\\\\"))
end

function _run_timed(cmd::Cmd)
    t0 = time_ns()
    run(cmd)
    return (time_ns() - t0) / 1e6
end

@inline function _count_lines(path::String)
    isfile(path) || return 0
    return length(readlines(path))
end

function _count_focus_lines(path::String)
    isfile(path) || return 0
    re = r"ForwardDiff|Models\.Conditions|Models\.var\"#residual!#140|_gap_conditions_with_model|gap_core_residual!"
    count = 0
    for line in eachline(path)
        occursin(re, line) && (count += 1)
    end
    return count
end

mkpath(TRACE_DIR)
const SYSIMAGE_PATH = _resolve_sysimage_path()
const LABEL = _resolve_label()
const SUMMARY_PATH = joinpath(TRACE_DIR, "solver_residual_$(LABEL)_summary.json")
const NO_SYS_TRACE = joinpath(TRACE_DIR, "solver_residual_$(LABEL)_no_sys.jl")
const WITH_SYS_TRACE = joinpath(TRACE_DIR, "solver_residual_$(LABEL)_with_sys.jl")
isfile(SYSIMAGE_PATH) || error("sysimage not found at $(SYSIMAGE_PATH)")

expr = _workload_expr()
no_sys_cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) --trace-compile=$(NO_SYS_TRACE) -e $(expr)`
with_sys_cmd = `$(Base.julia_cmd()) --sysimage=$(SYSIMAGE_PATH) --project=$(PROJECT_ROOT) --trace-compile=$(WITH_SYS_TRACE) -e $(expr)`

println("[solver-residual-sysimage] running no-sys workload...")
no_sys_ms = _run_timed(no_sys_cmd)

println("[solver-residual-sysimage] running with-sys workload...")
with_sys_ms = _run_timed(with_sys_cmd)

no_sys_lines = _count_lines(NO_SYS_TRACE)
with_sys_lines = _count_lines(WITH_SYS_TRACE)
no_sys_focus = _count_focus_lines(NO_SYS_TRACE)
with_sys_focus = _count_focus_lines(WITH_SYS_TRACE)

summary = (
    workload="conditions_build_residual_forwarddiff_jacobian",
    label=LABEL,
    sysimage_path=SYSIMAGE_PATH,
    no_sys=(wall_ms=no_sys_ms, trace_lines=no_sys_lines, focus_trace_lines=no_sys_focus),
    with_sys=(wall_ms=with_sys_ms, trace_lines=with_sys_lines, focus_trace_lines=with_sys_focus),
    improvement=(wall_ms=no_sys_ms - with_sys_ms, trace_lines=no_sys_lines - with_sys_lines, focus_trace_lines=no_sys_focus - with_sys_focus),
)

open(SUMMARY_PATH, "w") do io
    JSON3.pretty(io, summary)
end

println("[solver-residual-sysimage] summary")
println("  no-sys  : $(round(no_sys_ms; digits=1)) ms, trace=$(no_sys_lines), focus=$(no_sys_focus)")
println("  with-sys: $(round(with_sys_ms; digits=1)) ms, trace=$(with_sys_lines), focus=$(with_sys_focus)")
println("  delta   : $(round(no_sys_ms - with_sys_ms; digits=1)) ms, trace=$(no_sys_lines - with_sys_lines), focus=$(no_sys_focus - with_sys_focus)")
println("  summary : $(SUMMARY_PATH)")
