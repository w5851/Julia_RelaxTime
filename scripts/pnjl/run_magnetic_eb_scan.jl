#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

const PNJL = Models.pnjl_module()

function main(; T_MeV::Float64=150.0, mu_MeV::Float64=60.0, eB_min_MeV2::Float64=5.0e3, eB_max_MeV2::Float64=4.0e4, points::Int=8, output::String=joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl_magnetic", "scan", "pnjl_magnetic_eb_scan.csv"))
    points >= 2 || error("points must be >= 2")
    mkpath(dirname(output))

    T_fm = T_MeV / Constants_PNJL.ħc_MeV_fm
    mu_fm = mu_MeV / Constants_PNJL.ħc_MeV_fm
    mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)

    eB_values = range(eB_min_MeV2, eB_max_MeV2; length=points)

    open(output, "w") do io
        println(io, "T_MeV,mu_MeV,eB_MeV2,omega_fm4,pressure_fm4,rho_u,rho_d,rho_s,baryon_rho0,n_max,G_B,nmax_rel_diff,nmax_converged")
        for eB_MeV2 in eB_values
            eB_fm2 = eB_MeV2 / (Constants_PNJL.ħc_MeV_fm^2)
            conf = default_magnetic_config(eB_fm2=eB_fm2)

            comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
            rho = calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
            nd = calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
            conv = magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=5e-3)

            @printf(io, "%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%s\n",
                T_MeV, mu_MeV, eB_MeV2,
                comp.omega, -comp.omega,
                rho[1], rho[2], rho[3], nd.baryon,
                comp.n_max, comp.G_B,
                conv.rel_diff,
                string(conv.converged),
            )
        end
    end

    println("magnetic eB scan saved to: $output")
end

main()
