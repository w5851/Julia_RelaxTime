#!/usr/bin/env julia

"""Audit PALC states with the complete production magnetic kernel.

The input is a saved single-anchor PALC result.  This script does not solve a
new branch.  At selected target states it evaluates the production vacuum and
Landau thermal derivatives at several longitudinal quadrature orders, so the
condensate self-consistency gap is kept separate from the stationarity
residual used by the current solver.
"""

using JSON
using ForwardDiff
using LinearAlgebra: norm, svdvals
using StaticArrays

include(joinpath(@__DIR__, "MagneticReference.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(ROOT, "src", "models", "Models.jl"))
end
const MT = Models.MagneticThermodynamics
const MI = MT.MagneticIntegrals
const HBARC = 197.3269804
const CHARGES = (2 / 3, 1 / 3, 1 / 3)

function production_condensate_gap(x, T, mu, eB, p_num)
    conf = MT.MagneticConfig(eB_fm2=eB, n_max=32, p_num=p_num, pz_max=40.0)
    params = conf.params
    comp = MT.calculate_magnetic_omega_components(
        x, mu, T, conf; p_num=p_num, n_max=32, pz_max=40.0,
    )
    loops = zeros(3)
    for f in 1:3
        q = CHARGES[f]
        mass0 = comp.masses[f]
        loops[f] = ForwardDiff.derivative(mass0) do mass
            -2 * params.N_color * Models.PNJLCore.vacuum_integral_with_cutoff(
                mass, convert(typeof(mass), params.Λ_inv_fm),
            ) + MT.omega_magnetic_mfir(
                mass, q, eB; zeta_num=conf.zeta_num, N_c=params.N_color,
            ) + MI.omegat_flavor_landau(
                mass, mu[f], T, x[4], x[5], q, eB;
                n_max=32, p_num=p_num, pz_max=40.0,
            )
        end
    end
    return loops - collect(x[1:3])
end

function audit_state(c, target, root, p_num)
    eB = c["eB_GeV2"] * (1000 / HBARC)^2
    T = c["T_MeV"] / HBARC
    mu = SVector{3, Float64}(fill(target / (3 * HBARC), 3))
    x = SVector{5, Float64}(root["x_state"])
    conf = MT.MagneticConfig(eB_fm2=eB, n_max=32, p_num=p_num, pz_max=40.0)
    model = Models.PNJLMagneticModel(Models.PNJLModel(), conf)
    gap = production_condensate_gap(x, T, mu, eB, p_num)
    stationarity = Models.magnetic_gap_residual(
        model, x, T, mu; p_num=p_num, n_max=32, pz_max=40.0, t_num=4,
    )
    rho = MT.calculate_magnetic_rho(x, mu, T, conf)
    mass_map = phi -> Models.calculate_mass_vec(model, phi)
    Jm = ForwardDiff.jacobian(mass_map, collect(x[1:3]))
    reference = MagneticReference.evaluate(x, T, mu, eB,
        MagneticReference.Parameters(), MagneticReference.Controls(order=48))
    reference_high = MagneticReference.evaluate(x, T, mu, eB,
        MagneticReference.Parameters(), MagneticReference.Controls(order=96))
    return (; p_num, condensate_gap=collect(gap), condensate_gap_norm=norm(gap),
        stationarity_norm=norm(stationarity), rho=collect(rho),
        rhoB_over_rho0=sum(rho) / (3 * 0.16),
        mass_jacobian_sigma_min=minimum(svdvals(Jm)),
        chain_rule_error=norm(stationarity[1:3] - Jm' * gap),
        reference_condensate_gap_norm=norm(reference.condensate_gap),
        reference_stationarity_norm=norm(reference.stationarity),
        reference_rho=reference.rho,
        rho_reference_delta=norm(rho-reference.rho),
        reference_order_rho_delta=norm(reference.rho-reference_high.rho),
        reference_order_gap_delta=norm(reference.condensate_gap-reference_high.condensate_gap))
end

function main(args=ARGS)
    length(args) == 2 || throw(ArgumentError("input PALC JSON and output JSON required"))
    input, output = args
    data = JSON.parsefile(input)
    c = only(data["cases"])
    targets = Any[]
    for target in c["targets"]
        target["muB_MeV"] in (900.0, 1000.0, 1050.0, 1100.0, 1150.0, 1200.0) || continue
        roots = haskey(target, "palc_roots") ? target["palc_roots"] : target["roots"]
        for (i, root) in enumerate(roots)
            checks = Any[]
            for p_num in (128, 512, 1024)
                println("audit muB=", target["muB_MeV"], " root=", i, " p_num=", p_num)
                flush(stdout)
                push!(checks, audit_state(c, target["muB_MeV"], root, p_num))
            end
            push!(targets, (; muB_MeV=target["muB_MeV"], root_index=i,
                x_state=root["x_state"], checks))
        end
    end
    mkpath(dirname(output))
    open(output, "w") do io
        JSON.print(io, (; julia_version=string(VERSION), source=input,
            eB_GeV2=c["eB_GeV2"], T_MeV=c["T_MeV"], targets), 2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
