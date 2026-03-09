using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GAP_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "rpnjl", "baseline_rpnjl_gap_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 10 || error("invalid baseline row: $line")
        push!(rows, (
            label=strip(cols[1]),
            T_MeV=parse(Float64, strip(cols[2])),
            mu_MeV=parse(Float64, strip(cols[3])),
            phi_u=parse(Float64, strip(cols[4])),
            phi_d=parse(Float64, strip(cols[5])),
            phi_s=parse(Float64, strip(cols[6])),
            Phi=parse(Float64, strip(cols[7])),
            PhiBar=parse(Float64, strip(cols[8])),
            omega=parse(Float64, strip(cols[9])),
            max_abs_residual=parse(Float64, strip(cols[10])),
        ))
    end

    isempty(rows) && error("baseline CSV has no data rows: $path")
    return rows
end

@testset "RPNJL gap fixedpoint regression" begin
    rows = _load_baseline(GAP_BASELINE_PATH)
    model = Models.create_model(:RPNJL; profile="default", physics_profile="default")
    p_num = 12
    t_num = 4
    xi = 0.0
    residual_norm_max = 1e-5
    rtol = 1e-6
    atol = 1e-10

    @test isapprox(model.ext.kappa, 0.1; rtol=0.0, atol=1e-12)

    for row in rows
        @testset "$(row.label)" begin
            T_fm = row.T_MeV / Constants_PNJL.ħc_MeV_fm
            mu_fm = row.mu_MeV / Constants_PNJL.ħc_MeV_fm
            st = Models.solve_gap(model, T_fm, mu_fm;
                p_num=p_num,
                t_num=t_num,
                xi=xi,
                residual_norm_max=residual_norm_max,
            )
            x = Models.state_vector(st)
            omega = Models.omega(model, st, T_fm, mu_fm; p_num=p_num, t_num=t_num, xi=xi)
            residual = Models.gap_residual(model, st, T_fm, mu_fm; p_num=p_num, t_num=t_num, xi=xi)

            @test isapprox(x[1], row.phi_u; rtol=rtol, atol=atol)
            @test isapprox(x[2], row.phi_d; rtol=rtol, atol=atol)
            @test isapprox(x[3], row.phi_s; rtol=rtol, atol=atol)
            @test isapprox(x[4], row.Phi; rtol=rtol, atol=atol)
            @test isapprox(x[5], row.PhiBar; rtol=rtol, atol=atol)
            @test isapprox(omega, row.omega; rtol=rtol, atol=atol)
            @test maximum(abs.(residual)) <= max(1e-6, row.max_abs_residual * 10)
        end
    end
end