using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

to_fm_inv(x_mev::Real) = Float64(x_mev) / HBARC_MEV_FM

@testset "models fixedrho precision guard regression" begin
    model = Models.create_model(:PNJL)
    cases = [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=130.0, rho_target=1.0),
    ]

    for case in cases
        T_fm = to_fm_inv(case.T_MeV)
        result = Models.solve_constraint(
            model,
            Models.FixedRho(case.rho_target),
            T_fm;
            seed_guess=PNJL.HADRON_SEED_8,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test result.converged
        @test isfinite(result.residual_norm)
        @test result.residual_norm <= 1e-6
        @test isfinite(result.rho_norm)
        @test abs(result.rho_norm - case.rho_target) <= 1e-6
    end
end
