using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const NJL2_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "njl", "baseline_njl2_gap_fixedpoints_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    rows = NamedTuple[]
    for line in readlines(path)[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        push!(rows, (
            label=strip(cols[1]), T_MeV=parse(Float64, strip(cols[2])), mu_MeV=parse(Float64, strip(cols[3])),
            phi_u=parse(Float64, strip(cols[4])), phi_d=parse(Float64, strip(cols[5])), phi_s=parse(Float64, strip(cols[6])),
            Phi=parse(Float64, strip(cols[7])), PhiBar=parse(Float64, strip(cols[8])), omega=parse(Float64, strip(cols[9])),
        ))
    end
    return rows
end

@testset "NJL2 gap fixedpoint regression" begin
    rows = _load_rows(NJL2_BASELINE_PATH)
    model = Models.create_model(:NJL2)
    for row in rows
        T_fm = row.T_MeV / 197.327
        mu_fm = row.mu_MeV / 197.327
        st = Models.solve_gap(model, T_fm, mu_fm; p_num=24, t_num=6, xi=0.0, residual_norm_max=1e-6)
        x = Models.state_vector(st)
        omega = Models.omega(model, st, T_fm, mu_fm; p_num=24, t_num=6, xi=0.0)
        @test isapprox(x[1], row.phi_u; rtol=1e-6, atol=1e-10)
        @test isapprox(x[2], row.phi_d; rtol=1e-6, atol=1e-10)
        @test isapprox(x[3], row.phi_s; rtol=1e-12, atol=1e-12)
        @test isapprox(x[4], row.Phi; rtol=1e-12, atol=1e-12)
        @test isapprox(x[5], row.PhiBar; rtol=1e-12, atol=1e-12)
        @test isapprox(omega, row.omega; rtol=1e-6, atol=1e-10)
    end
end