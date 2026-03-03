#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

const PNJL = Models.pnjl_module()

const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl_magnetic", "stability", "pnjl_magnetic_stability_scan_latest.csv")
const DEFAULT_FAILURES = joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl_magnetic", "stability", "pnjl_magnetic_stability_failures_latest.csv")

function main(; output::String=DEFAULT_OUTPUT, failures_output::String=DEFAULT_FAILURES)
    mkpath(dirname(output))
    mkpath(dirname(failures_output))

    T_values = [130.0, 145.0, 160.0, 175.0]             # MeV
    mu_values = [0.0, 60.0, 120.0, 180.0, 240.0]        # MeV
    eB_values = [8.0e3, 1.5e4, 2.5e4, 3.5e4]            # MeV^2

    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)

    rows = NamedTuple[]
    failures = NamedTuple[]

    for T_MeV in T_values, mu_MeV in mu_values, eB_MeV2 in eB_values
        T_fm = T_MeV / Constants_PNJL.ħc_MeV_fm
        mu_fm = mu_MeV / Constants_PNJL.ħc_MeV_fm
        eB_fm2 = eB_MeV2 / (Constants_PNJL.ħc_MeV_fm^2)
        mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)

        status = "pass"
        reason = ""
        omega = NaN
        pressure = NaN
        baryon = NaN
        nmax = -1
        rel_diff = NaN

        try
            conf = default_magnetic_config(eB_fm2=eB_fm2)
            comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
            nd = calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
            conv = magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=3e-2)

            omega = comp.omega
            pressure = -comp.omega
            baryon = nd.baryon
            nmax = comp.n_max
            rel_diff = conv.rel_diff

            finite_ok = isfinite(omega) && isfinite(pressure) && isfinite(baryon)
            bounded_ok = abs(omega) < 1e6 && abs(pressure) < 1e6 && abs(baryon) < 1e6
            conv_ok = conv.converged

            if !(finite_ok && bounded_ok && conv_ok)
                status = "fail"
                reasons = String[]
                finite_ok || push!(reasons, "nonfinite")
                bounded_ok || push!(reasons, "out_of_bound")
                conv_ok || push!(reasons, "nmax_not_converged")
                reason = join(reasons, ";")
            end
        catch err
            status = "fail"
            reason = "exception: " * sprint(showerror, err)
        end

        row = (
            T_MeV=T_MeV,
            mu_MeV=mu_MeV,
            eB_MeV2=eB_MeV2,
            status=status,
            reason=reason,
            omega_fm4=omega,
            pressure_fm4=pressure,
            baryon_rho0=baryon,
            n_max=nmax,
            nmax_rel_diff=rel_diff,
        )
        push!(rows, row)
        status == "fail" && push!(failures, row)
    end

    open(output, "w") do io
        println(io, "T_MeV,mu_MeV,eB_MeV2,status,reason,omega_fm4,pressure_fm4,baryon_rho0,n_max,nmax_rel_diff")
        for r in rows
            @printf(io, "%.6f,%.6f,%.6f,%s,%s,%.16e,%.16e,%.16e,%d,%.16e\n",
                r.T_MeV, r.mu_MeV, r.eB_MeV2,
                r.status, replace(r.reason, "," => "|"),
                r.omega_fm4, r.pressure_fm4, r.baryon_rho0, r.n_max, r.nmax_rel_diff)
        end
    end

    open(failures_output, "w") do io
        println(io, "T_MeV,mu_MeV,eB_MeV2,status,reason,omega_fm4,pressure_fm4,baryon_rho0,n_max,nmax_rel_diff")
        for r in failures
            @printf(io, "%.6f,%.6f,%.6f,%s,%s,%.16e,%.16e,%.16e,%d,%.16e\n",
                r.T_MeV, r.mu_MeV, r.eB_MeV2,
                r.status, replace(r.reason, "," => "|"),
                r.omega_fm4, r.pressure_fm4, r.baryon_rho0, r.n_max, r.nmax_rel_diff)
        end
    end

    println("stability scan done: total=$(length(rows)), failures=$(length(failures))")
    println("output: " * output)
    println("failures: " * failures_output)
end

main()
