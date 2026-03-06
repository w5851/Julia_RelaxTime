using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MAGNETIC_SCOPE = lowercase(get(ENV, "REGRESSION_MAGNETIC_SCOPE", "smoke"))
const BASELINE_FILE = MAGNETIC_SCOPE == "nightly" ?
    "baseline_pnjl_magnetic_fixedpoints_nightly_v1.csv" :
    "baseline_pnjl_magnetic_fixedpoints_smoke_v1.csv"
const MAGNETIC_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", BASELINE_FILE)

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 15 || error("invalid baseline row: $line")
        push!(rows, (
            scope=strip(cols[1]),
            T=parse(Float64, strip(cols[2])),
            mu=parse(Float64, strip(cols[3])),
            eB=parse(Float64, strip(cols[4])),
            xi=parse(Float64, strip(cols[5])),
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

@testset "PNJL magnetic fixedpoint regression ($(MAGNETIC_SCOPE))" begin
    rows = _load_rows(MAGNETIC_BASELINE_PATH)
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)

    rtol = MAGNETIC_SCOPE == "nightly" ? 5e-2 : 8e-2
    atol = MAGNETIC_SCOPE == "nightly" ? 1e-7 : 1e-6

    for row in rows
        T_fm = row.T / Constants_PNJL.ħc_MeV_fm
        mu_fm = row.mu / Constants_PNJL.ħc_MeV_fm
        eB_fm2 = row.eB / (Constants_PNJL.ħc_MeV_fm^2)
        mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)

        conf = PNJL.default_magnetic_config(eB_fm2=eB_fm2)
        comp = PNJL.calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
        rho = PNJL.calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
        nd = PNJL.calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
        conv = PNJL.magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=3e-2)

        @test isapprox(comp.omega, row.omega; rtol=rtol, atol=atol)
        @test isapprox(-comp.omega, row.pressure; rtol=rtol, atol=atol)
        @test isapprox(rho[1], row.rho_u; rtol=rtol, atol=atol)
        @test isapprox(rho[2], row.rho_d; rtol=rtol, atol=atol)
        @test isapprox(rho[3], row.rho_s; rtol=rtol, atol=atol)
        @test isapprox(nd.baryon, row.baryon; rtol=rtol, atol=atol)
        @test isapprox(comp.G_B, row.G_B; rtol=rtol, atol=atol)
        @test comp.n_max == row.n_max
        @test conv.converged == row.converged
        @test isapprox(conv.rel_diff, row.rel_diff; rtol=2e-1, atol=1e-8)
    end
end