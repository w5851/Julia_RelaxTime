using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_ENTRY = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :create_model) && isdefined(Main.Models, :omega))
    Base.include(Main, _MODELS_ENTRY)
end

@testset "Models bridge smoke: PNJL vs RPNJL(degenerate)" begin
    pnjl = Main.Models.create_model(:PNJL)
    rpnjl_deg = Main.Models.create_model(:RPNJL; use_rpnjl_extensions=false)
    rpnjl_ext = Main.Models.create_model(:RPNJL; use_rpnjl_extensions=true)

    solver = Main.Models.NLsolveGapSolver(
        method=:trust_region,
        jacobian=:finite,
        xtol=1e-10,
        ftol=1e-10,
    )

    points = (
        (T=0.15, mu=0.00),
        (T=0.18, mu=0.00),
        (T=0.20, mu=0.00),
        (T=0.18, mu=0.10),
        (T=0.22, mu=0.00),
        (T=0.24, mu=0.00),
        (T=0.18, mu=0.15),
        (T=0.20, mu=0.20),
    )

    omega_ext_baseline = Dict(
        (0.15, 0.00) => -21.95531477423188,
        (0.18, 0.00) => -21.955196741429955,
        (0.20, 0.00) => -21.95507724709505,
        (0.18, 0.10) => -21.955196741429944,
        (0.22, 0.00) => -21.954915987128526,
        (0.24, 0.00) => -21.954704205654604,
        (0.18, 0.15) => -21.955196741429948,
        (0.20, 0.20) => -21.95507724709505,
    )
    baseline_rtol = 5e-6
    baseline_atol = 5e-7

    for pt in points
        T_fm = pt.T
        mu = pt.mu

        x_pnjl = Main.Models.solve_gap(
            pnjl,
            T_fm,
            mu;
            solver_backend=:models,
            solver=solver,
            p_num=8,
            t_num=4,
            xi=0.0,
            residual_norm_max=1e-4,
        )

        x_rpnjl_deg = Main.Models.solve_gap(
            rpnjl_deg,
            T_fm,
            mu;
            solver=solver,
            p_num=8,
            t_num=4,
            xi=0.0,
            residual_norm_max=1e-4,
        )

        x_rpnjl_ext = Main.Models.solve_gap(
            rpnjl_ext,
            T_fm,
            mu;
            solver=solver,
            p_num=8,
            t_num=4,
            xi=0.0,
            residual_norm_max=1e-4,
        )

        ω_pnjl = Main.Models.omega(pnjl, x_pnjl, T_fm, mu; p_num=8, t_num=4, xi=0.0)
        ω_rpnjl_deg = Main.Models.omega(rpnjl_deg, x_rpnjl_deg, T_fm, mu; p_num=8, t_num=4, xi=0.0)
        ω_rpnjl_ext = Main.Models.omega(rpnjl_ext, x_rpnjl_ext, T_fm, mu; p_num=8, t_num=4, xi=0.0)

        @test isfinite(ω_pnjl)
        @test isfinite(ω_rpnjl_deg)
        @test isfinite(ω_rpnjl_ext)

        @test isapprox(ω_rpnjl_deg, ω_pnjl; rtol=5e-2, atol=1e-6)

        ext_delta = abs(ω_rpnjl_ext - ω_rpnjl_deg)
        @test ext_delta > 1e-10

        ω_ref = omega_ext_baseline[(T_fm, mu)]
        @test isapprox(ω_rpnjl_ext, ω_ref; rtol=baseline_rtol, atol=baseline_atol)
    end
end
