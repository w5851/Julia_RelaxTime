#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Libdl

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TRACE_DIR = joinpath(PROJECT_ROOT, "build", "trace")
const SYSIMAGE_PATH = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime.$(Libdl.dlext)")

const _TRACE_EXPR = "ENV[\"UNIT_FILES\"]=\"pnjl/test_conserved_charge_susceptibilities.jl\"; include(\"tests/unit/runtests.jl\")"
const NO_SYS_CMD = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) --trace-compile=$(joinpath(TRACE_DIR, "unit_conserved_no_sys.jl")) -e $(_TRACE_EXPR)`
const WITH_SYS_CMD = `$(Base.julia_cmd()) --sysimage=$(SYSIMAGE_PATH) --project=$(PROJECT_ROOT) --trace-compile=$(joinpath(TRACE_DIR, "unit_conserved_with_sys.jl")) -e $(_TRACE_EXPR)`

# Recent Ubuntu CI runs (Julia 1.12.5) observed up to 1047 lines without
# the sysimage and 402 with it. Cap the work left in each mode; covering
# more statements with the sysimage can increase both diagnostic deltas.
const MAX_NO_SYS_LINES = 1100
const MAX_WITH_SYS_LINES = 450

function trace_budget_errors(no_sys_count::Int, with_sys_count::Int)
    errors = String[]
    no_sys_count <= MAX_NO_SYS_LINES || push!(errors, "no-sys lines exceed budget: $no_sys_count > $MAX_NO_SYS_LINES")
    with_sys_count <= MAX_WITH_SYS_LINES || push!(errors, "with-sys lines exceed budget: $with_sys_count > $MAX_WITH_SYS_LINES")
    return errors
end

mkpath(TRACE_DIR)
isfile(SYSIMAGE_PATH) || error("sysimage not found at $(SYSIMAGE_PATH); run scripts/dev/build_sysimage.jl first")

println("[trace-budget] collecting no-sys trace...")
run(NO_SYS_CMD)

println("[trace-budget] collecting with-sys trace...")
run(WITH_SYS_CMD)

no_lines = readlines(joinpath(TRACE_DIR, "unit_conserved_no_sys.jl"))
with_lines = readlines(joinpath(TRACE_DIR, "unit_conserved_with_sys.jl"))

with_set = Set(with_lines)
delta = [line for line in no_lines if !(line in with_set)]

focus_re = r"ForwardDiff|NLSolversBase|Conserved|Thermo|chi|cumulant|baryon|StaticArrays"
focus = [line for line in delta if occursin(focus_re, line)]

println("[trace-budget] no-sys lines   = $(length(no_lines))")
println("[trace-budget] with-sys lines = $(length(with_lines))")
println("[trace-budget] delta lines    = $(length(delta))")
println("[trace-budget] focus delta    = $(length(focus))")

errors = trace_budget_errors(length(no_lines), length(with_lines))

open(joinpath(TRACE_DIR, "unit_conserved_delta_focus.jl"), "w") do io
    for line in focus
        println(io, line)
    end
end

if isempty(errors)
    println("[trace-budget] OK")
    exit(0)
end

println("[trace-budget] FAILED")
for err in errors
    println("  - " * err)
end
exit(1)
