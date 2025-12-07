#!/usr/bin/env julia

try
    using Plots
catch
    import Pkg
    Pkg.add("Plots")
    using Plots
end

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

include(joinpath(PROJECT_ROOT, "tests", "analysis", "relaxtime", "integrand_vs_s_average_rate.jl"))

params = build_params()
rows = sample_omega_terms(:ssbar_to_uubar; params=params, p_nodes=6, angle_nodes=4, phi_nodes=8, n_sigma_points=16)

s = [r.s for r in rows]
w = [r.weighted for r in rows]

outdir = joinpath(PROJECT_ROOT, "output", "figures")
mkpath(outdir)

plt = scatter(s, w; markersize=3, markerstrokewidth=0, xlabel="s", ylabel="weighted", title="weighted vs s")
savefig(joinpath(outdir, "weighted_vs_s.png"))
println("Saved plot to ", joinpath(outdir, "weighted_vs_s.png"))
