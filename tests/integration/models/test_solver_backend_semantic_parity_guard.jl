using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()
const HBARC_MEV_FM = 197.327

to_fm_inv(x_mev) = Float64(x_mev) / HBARC_MEV_FM

function _solve_fixedrho_backend(backend::Symbol, T_MeV::Real, rho_target::Real; p_num::Int=8, t_num::Int=4)
    T_fm = to_fm_inv(T_MeV)
    mode = P.FixedRho(Float64(rho_target))
    seed = copy(P.HADRON_SEED_8)

    if backend === :models
        model = Models.create_model(:PNJL)
        return Models.solve_constraint(
            model,
            Models.FixedRho(Float64(rho_target)),
            T_fm;
            seed_guess=seed,
            p_num=p_num,
            t_num=t_num,
            residual_norm_max=1e-6,
        )
    elseif backend === :legacy
        return P.solve(mode, T_fm; seed_strategy=P.DefaultSeed(seed, seed, :hadron), p_num=p_num, t_num=t_num, residual_norm_max=1e-6)
    elseif backend === :auto
        return P.TrhoScan._solve_with_models(mode, T_fm;
            xi=0.0,
            model_kind=:PNJL,
            seed_strategy=P.DefaultSeed(seed, seed, :hadron),
            p_num=p_num,
            t_num=t_num,
            residual_norm_max=1e-6,
        )
    end

    throw(ArgumentError("unsupported backend $(backend)"))
end

@testset "solver backend semantic parity guard" begin
    cases = [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=130.0, rho_target=1.0),
    ]

    for case in cases
        legacy = _solve_fixedrho_backend(:legacy, case.T_MeV, case.rho_target)
        models = _solve_fixedrho_backend(:models, case.T_MeV, case.rho_target)
        auto = _solve_fixedrho_backend(:auto, case.T_MeV, case.rho_target)

        @test legacy.converged
        @test models.converged
        @test auto.converged

        @test abs(models.rho_norm - case.rho_target) <= 1e-6
        @test abs(auto.rho_norm - case.rho_target) <= 1e-6

        @test isapprox(auto.pressure, models.pressure; rtol=1e-9, atol=1e-10)
        @test isapprox(auto.entropy, models.entropy; rtol=1e-9, atol=1e-10)
        @test isapprox(auto.energy, models.energy; rtol=1e-9, atol=1e-10)
        @test isapprox(auto.rho_norm, models.rho_norm; rtol=1e-9, atol=1e-10)

        @test isapprox(models.pressure, legacy.pressure; rtol=1e-6, atol=1e-8)
        @test isapprox(models.entropy, legacy.entropy; rtol=1e-6, atol=1e-8)
        @test isapprox(models.energy, legacy.energy; rtol=1e-6, atol=1e-8)
        @test abs(models.rho_norm - legacy.rho_norm) <= 1e-6
    end
end

@testset "FixedRho ProblemSpec converges with default seed pool" begin
    model = Models.create_model(:PNJL)
    T_fm = to_fm_inv(110.0)
    mode = Models.FixedRho(0.6)
    seed = copy(P.HADRON_SEED_8)

    raw = Models.solve_constraint(
        model,
        mode,
        T_fm;
        seed_guess=seed,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
    )

    @test raw.converged
    @test abs(raw.rho_norm - 0.6) <= 1e-6
end
