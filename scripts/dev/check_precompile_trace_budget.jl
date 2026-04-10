#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TRACE_DIR = joinpath(PROJECT_ROOT, "build", "trace")

const _TRACE_EXPR = "ENV[\"UNIT_FILES\"]=\"pnjl/test_conserved_charge_susceptibilities.jl\"; include(\"tests/unit/runtests.jl\")"
const NO_SYS_CMD = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) --trace-compile=$(joinpath(TRACE_DIR, "unit_conserved_no_sys.jl")) -e $(_TRACE_EXPR)`
const WITH_SYS_CMD = `$(Base.julia_cmd()) --sysimage=$(joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime.dll")) --project=$(PROJECT_ROOT) --trace-compile=$(joinpath(TRACE_DIR, "unit_conserved_with_sys.jl")) -e $(_TRACE_EXPR)`

const MAX_WITH_SYS_LINES = 1000
const MAX_DELTA_LINES = 450
const MAX_FOCUS_DELTA_LINES = 120

mkpath(TRACE_DIR)

println("[trace-budget] collecting no-sys trace...")
run(NO_SYS_CMD)

println("[trace-budget] collecting with-sys trace...")
run(WITH_SYS_CMD)

no_lines = readlines(joinpath(TRACE_DIR, "unit_conserved_no_sys.jl"))
with_lines = readlines(joinpath(TRACE_DIR, "unit_conserved_with_sys.jl"))

no_set = Set(no_lines)
with_set = Set(with_lines)
delta = [line for line in no_lines if !(line in with_set)]

focus_re = r"ForwardDiff|NLSolversBase|Conserved|Thermo|chi|cumulant|baryon|StaticArrays"
focus = [line for line in delta if occursin(focus_re, line)]

println("[trace-budget] no-sys lines   = $(length(no_lines))")
println("[trace-budget] with-sys lines = $(length(with_lines))")
println("[trace-budget] delta lines    = $(length(delta))")
println("[trace-budget] focus delta    = $(length(focus))")

errors = String[]
length(with_lines) <= MAX_WITH_SYS_LINES || push!(errors, "with-sys lines exceed budget: $(length(with_lines)) > $(MAX_WITH_SYS_LINES)")
length(delta) <= MAX_DELTA_LINES || push!(errors, "delta lines exceed budget: $(length(delta)) > $(MAX_DELTA_LINES)")
length(focus) <= MAX_FOCUS_DELTA_LINES || push!(errors, "focus delta exceeds budget: $(length(focus)) > $(MAX_FOCUS_DELTA_LINES)")

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
