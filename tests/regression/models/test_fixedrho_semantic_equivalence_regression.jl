using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

to_fm_inv(x_mev::Real) = Float64(x_mev) / HBARC_MEV_FM

@testset "fixedrho semantic equivalence regression" begin
    points = [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=130.0, rho_target=1.0),
    ]

    for point in points
        T_fm = to_fm_inv(point.T_MeV)
        mode = P.FixedRho(point.rho_target)
        seed = copy(P.HADRON_SEED_8)

        legacy = P.solve(mode, T_fm;
            seed_strategy=P.DefaultSeed(seed, seed, :hadron),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        model = Models.create_model(:PNJL)
        models = Models.solve_constraint(
            model,
            Models.FixedRho(point.rho_target),
            T_fm;
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test legacy.converged
        @test models.converged
        @test isapprox(models.rho_norm, legacy.rho_norm; rtol=2e-6, atol=1e-8)
        @test isapprox(models.pressure, legacy.pressure; rtol=1e-6, atol=1e-8)
        @test isapprox(models.entropy, legacy.entropy; rtol=1e-6, atol=1e-8)
        @test isapprox(models.energy, legacy.energy; rtol=1e-6, atol=1e-8)
    end
end
