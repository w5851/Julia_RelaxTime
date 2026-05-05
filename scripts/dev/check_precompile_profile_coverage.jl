#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
using .Models

const REQUIRED_PROFILES = Dict(
    :smoke => Set([:gap_solver_ad]),
    :test => Set([:gap_solver_ad, :thermo_derivatives_ad]),
    :scan => Set([:gap_solver_ad, :solver_residual_ad, :thermo_derivatives_ad, :transport_point_api, :scan_pipeline_cli]),
    :core => Set([:gap_solver_ad, :solver_residual_ad, :thermo_derivatives_ad, :conserved_charge_highorder, :ad_shape_stabilization, :transport_point_api, :scan_pipeline_cli]),
    :full => Set([:gap_solver_ad, :solver_residual_ad, :thermo_derivatives_ad, :conserved_charge_highorder, :ad_shape_stabilization, :transport_point_api, :scan_pipeline_cli]),
)

all_capabilities = Set(Models.list_precompile_capabilities())
errors = String[]

for (profile, required) in REQUIRED_PROFILES
    listed = Set(Models.list_precompile_profile(profile))

    missing = setdiff(required, listed)
    isempty(missing) || push!(errors, "profile=$(profile) missing capabilities: $(collect(missing))")

    unknown = setdiff(listed, all_capabilities)
    isempty(unknown) || push!(errors, "profile=$(profile) has unknown capabilities: $(collect(unknown))")
end

if isempty(errors)
    println("[precompile-profile-coverage] OK")
    println("  capabilities = $(collect(all_capabilities))")
    for profile in sort!(collect(keys(REQUIRED_PROFILES)))
        println("  profile $(profile) -> $(Models.list_precompile_profile(profile))")
    end
    exit(0)
end

println("[precompile-profile-coverage] FAILED")
for err in errors
    println("  - " * err)
end
exit(1)
