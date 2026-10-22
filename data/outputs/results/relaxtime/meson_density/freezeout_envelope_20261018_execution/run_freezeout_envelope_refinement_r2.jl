#!/usr/bin/env julia
# Two independent jobs share the existing runner; no new numerical implementation.
include(joinpath(@__DIR__, "run_charged_gbu_freezeout_comparison.jl"))
using JSON3

function run_envelope_refinement(args)
    length(args) == 3 || throw(ArgumentError("usage: SELECTION_JSON NUMERICAL_PLAN_JSON PART(0|1)"))
    selection_path, plan_path, part = args
    part in ("0", "1") || throw(ArgumentError("part must be 0 or 1"))
    selection = JSON3.read(read(selection_path, String))
    plan = JSON3.read(read(plan_path, String))
    names = [String(c.name) for c in selection.candidates]
    length(names) >= 2 || throw(ArgumentError("refinement plan requires at least two selected curves"))
    profile_dir = joinpath(dirname(abspath(selection_path)), "profiles")
    common = (; profile_dir=profile_dir, parameter_policy="literature_envelope_exploration")
    runner = Main.ChargedGBUFreezeoutComparison
    if part == "0"
        runner.run_scan(; common..., output=String(plan.standard_output), profiles=names,
            energies=Float64.(plan.standard_energies_GeV), resolution="standard")
        runner.run_scan(; common..., output=String(plan.fine_output_0), profiles=names[1:1],
            energies=Float64.(plan.fine_energies_GeV), resolution="fine")
    else
        runner.run_scan(; common..., output=String(plan.fine_output_1), profiles=names[2:end],
            energies=Float64.(plan.fine_energies_GeV), resolution="fine")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_envelope_refinement(ARGS)
end
