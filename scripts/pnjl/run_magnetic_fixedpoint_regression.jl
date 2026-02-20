#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_magnetic_fixedpoints_nightly_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL

function load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV empty: $path")
    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) == 15 || error("invalid baseline row: $line")
        push!(rows, (
            T=parse(Float64, strip(cols[2])),
            mu=parse(Float64, strip(cols[3])),
            eB=parse(Float64, strip(cols[4])),
            omega=parse(Float64, strip(cols[6])),
            pressure=parse(Float64, strip(cols[7])),
            rho_u=parse(Float64, strip(cols[8])),
            rho_d=parse(Float64, strip(cols[9])),
            rho_s=parse(Float64, strip(cols[10])),
            baryon=parse(Float64, strip(cols[11])),
            n_max=parse(Int, strip(cols[12])),
            G_B=parse(Float64, strip(cols[13])),
            rel_diff=parse(Float64, strip(cols[14])),
            converged=lowercase(strip(cols[15])) == "true",
        ))
    end
    return rows
end

function main()
    rows = load_rows(BASELINE_PATH)
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)

    rtol = 5e-2
    atol = 1e-7

    total = 0
    pass = 0

    for row in rows
        total += 1
        T_fm = row.T / Constants_PNJL.ħc_MeV_fm
        mu_fm = row.mu / Constants_PNJL.ħc_MeV_fm
        eB_fm2 = row.eB / (Constants_PNJL.ħc_MeV_fm^2)
        mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)

        conf = default_magnetic_config(eB_fm2=eB_fm2)
        comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
        rho = calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
        nd = calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
        conv = magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=3e-2)

        checks = [
            isapprox(comp.omega, row.omega; rtol=rtol, atol=atol),
            isapprox(-comp.omega, row.pressure; rtol=rtol, atol=atol),
            isapprox(rho[1], row.rho_u; rtol=rtol, atol=atol),
            isapprox(rho[2], row.rho_d; rtol=rtol, atol=atol),
            isapprox(rho[3], row.rho_s; rtol=rtol, atol=atol),
            isapprox(nd.baryon, row.baryon; rtol=rtol, atol=atol),
            isapprox(comp.G_B, row.G_B; rtol=rtol, atol=atol),
            comp.n_max == row.n_max,
            conv.converged == row.converged,
            isapprox(conv.rel_diff, row.rel_diff; rtol=2e-1, atol=1e-8),
        ]

        ok = all(checks)
        pass += ok ? 1 : 0
        @printf("[magnetic-nightly] T=%.1f mu=%.1f eB=%.1f -> %s\n", row.T, row.mu, row.eB, ok ? "pass" : "FAIL")
    end

    println("magnetic nightly regression: $(pass)/$(total) passed")
    pass == total || error("magnetic nightly regression failed")
end

main()
