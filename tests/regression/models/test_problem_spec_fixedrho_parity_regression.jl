using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

to_fm_inv(x_mev::Real) = Float64(x_mev) / HBARC_MEV_FM

@testset "problem spec fixedrho parity regression" begin
    model = Models.create_model(:PNJL)
    points = [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=130.0, rho_target=1.0),
    ]

    for point in points
        mode = Models.FixedRho(point.rho_target)
        spec = Models.build_problem_spec(mode)
        T_fm = to_fm_inv(point.T_MeV)
        seed = Models.pnjl_module().HADRON_SEED_8

        direct = Models.solve_constraint(
            model,
            mode,
            T_fm;
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        via_spec = Models.solve_constraint(
            model,
            mode,
            T_fm;
            problem_spec=spec,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test direct.converged == via_spec.converged
        @test isapprox(direct.rho_norm, via_spec.rho_norm; rtol=1e-10, atol=1e-12)
        @test isapprox(direct.pressure, via_spec.pressure; rtol=1e-10, atol=1e-12)
        @test isapprox(direct.entropy, via_spec.entropy; rtol=1e-10, atol=1e-12)
        @test isapprox(direct.energy, via_spec.energy; rtol=1e-10, atol=1e-12)
        @test isapprox(direct.residual_norm, via_spec.residual_norm; rtol=1e-10, atol=1e-12)
    end
end
