#!/usr/bin/env julia

"""Correct saved PALC target roots at a second production quadrature order."""

using JSON
using LinearAlgebra: norm
using StaticArrays

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(ROOT, "src", "models", "Models.jl"))
end
const HBARC = 197.3269804

function main(args=ARGS)
    length(args) == 3 || throw(ArgumentError("input JSON, p_num and output JSON required"))
    input, p_arg, output = args
    p_num = parse(Int, p_arg)
    data = JSON.parsefile(input)
    cases = Any[]
    for c in data["cases"]
        eB = c["eB_GeV2"] * (1000 / HBARC)^2
        model = Models.PNJLMagneticModel(eB_fm2=eB, p_num=p_num, n_max=32, pz_max=40.0)
        T = c["T_MeV"] / HBARC
        targets = Any[]
        for target in c["targets"]
            roots = Any[]
            for (i, source) in enumerate(target["palc_roots"])
                mu = Models.normalize_mu_vec(target["muB_MeV"] / (3 * HBARC))
                result = Models.solve_magnetic_gap(
                    model, T, mu;
                    p_num=p_num, t_num=4, pz_max=40.0, n_max=32,
                    initial_guess=source["x_state"], include_default_seeds=false,
                    method=:trust_region, fallback_method=:newton, iterations=260,
                    residual_norm_max=2e-7, root_merge_tol=2e-5,
                )
                converged = result.converged && !isempty(result.candidates)
                if converged
                    cand = result.candidates[result.selected_index]
                    push!(roots, (; source_index=i, converged=true,
                        x_state=collect(cand.x_state), residual_norm=cand.residual_norm,
                        omega=cand.omega, physical=cand.physical,
                        iterations=cand.iterations))
                else
                    push!(roots, (; source_index=i, converged=false,
                        x_state=nothing, residual_norm=nothing, omega=nothing,
                        physical=false, iterations=nothing))
                end
                println("eB=", c["eB_GeV2"], " muB=", target["muB_MeV"],
                    " source=", i, " converged=", converged,
                    " candidates=", length(result.candidates))
                flush(stdout)
            end
            unique_roots = Any[]
            for root in roots
                root.converged || continue
                any(norm(root.x_state - other.x_state) < 2e-5 for other in unique_roots) ||
                    push!(unique_roots, root)
            end
            push!(targets, (; muB_MeV=target["muB_MeV"],
                source_root_count=length(target["palc_roots"]),
                corrected_root_count=length(unique_roots), roots))
        end
        push!(cases, (; eB_GeV2=c["eB_GeV2"], T_MeV=c["T_MeV"], p_num,
            targets))
    end
    mkpath(dirname(output))
    open(output, "w") do io
        JSON.print(io, (; julia_version=string(VERSION), source=input, cases), 2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
