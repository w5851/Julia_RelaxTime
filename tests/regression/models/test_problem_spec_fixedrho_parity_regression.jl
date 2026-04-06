using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

to_fm_inv(x_mev::Real) = Float64(x_mev) / HBARC_MEV_FM

@testset "problem spec fixedrho continuity-like seeds with joint only" begin
    model = Models.create_model(:PNJL)
    base_seed = copy(Models.pnjl_module().HADRON_SEED_8)
    points = [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=100.0, rho_target=0.4),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=120.0, rho_target=0.8),
        (T_MeV=130.0, rho_target=1.0),
    ]

    seed = copy(base_seed)
    for point in points
        mode = Models.FixedRho(point.rho_target)
        spec = Models.build_problem_spec(mode)
        T_fm = to_fm_inv(point.T_MeV)

        result = Models.solve_constraint(
            model,
            mode,
            T_fm;
            problem_spec=spec,
            fixedrho_joint_solve=true,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            nlsolve_method=:trust_region,
            trust_region_fallback=true,
            fallback_method=:trust_region,
            residual_norm_max=1e-6,
        )

        @test result.fixedrho_joint_solve_requested
        @test result.fixedrho_joint_solve_active
        @test result.converged
        @test result.residual_norm <= 1e-6

        seed = copy(result.solution)
    end
end

@testset "problem spec fixedrho parity regression" begin
    model = Models.create_model(:PNJL)
    points = [
        (T_MeV=90.0, rho_target=0.2),
        (T_MeV=110.0, rho_target=0.6),
        (T_MeV=130.0, rho_target=1.0),
    ]

    rolling_seed = copy(Models.pnjl_module().HADRON_SEED_8)

    for point in points
        mode = Models.FixedRho(point.rho_target)
        spec = Models.build_problem_spec(mode)
        T_fm = to_fm_inv(point.T_MeV)
        seed = copy(rolling_seed)

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

        via_spec_joint = Models.solve_constraint(
            model,
            mode,
            T_fm;
            problem_spec=spec,
            fixedrho_joint_solve=true,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test direct.converged == via_spec.converged
        @test isapprox(direct.rho_norm, via_spec.rho_norm; rtol=1e-6, atol=1e-8)
        @test isapprox(direct.pressure, via_spec.pressure; rtol=1e-6, atol=1e-8)
        @test isapprox(direct.entropy, via_spec.entropy; rtol=1e-6, atol=1e-8)
        @test isapprox(direct.energy, via_spec.energy; rtol=1e-6, atol=1e-8)
        @test isapprox(direct.residual_norm, via_spec.residual_norm; rtol=1e-6, atol=1e-10)

        @test via_spec_joint.fixedrho_joint_solve_requested
        @test via_spec_joint.fixedrho_joint_solve_active
        @test !via_spec_joint.fixedrho_joint_fallback
        @test via_spec_joint.selection_reason in (
            :pressure_max_under_constraints,
            :residual_min_under_constraints,
            :no_candidate_passed_constraints,
        )
        @test isfinite(via_spec_joint.residual_norm)
        @test via_spec_joint.residual_norm <= 1e-3

        if via_spec_joint.converged
            rolling_seed = copy(via_spec_joint.solution)
        end
    end
end

@testset "problem spec fixedrho joint parity gates" begin
    model = Models.create_model(:PNJL)
    T_fm = to_fm_inv(100.0)
    mode = Models.FixedRho(0.2)
    spec = Models.build_problem_spec(mode)
    seed = copy(Models.pnjl_module().HADRON_SEED_8)

    via_joint_no_fallback = Models.solve_constraint(
        model,
        mode,
        T_fm;
        problem_spec=spec,
        fixedrho_joint_solve=true,
        seed_guess=seed,
        seed_candidates=(seed,),
        p_num=8,
        t_num=4,
        nlsolve_method=:trust_region,
        trust_region_fallback=false,
        residual_norm_max=1e-6,
    )

    @test via_joint_no_fallback.converged
    @test via_joint_no_fallback.fixedrho_joint_solve_active
    @test haskey(via_joint_no_fallback, :fixedrho_joint_selected_method)
    @test via_joint_no_fallback.fixedrho_joint_selected_method == :trust_region
    @test haskey(via_joint_no_fallback, :fixedrho_joint_fallback_used)
    @test !via_joint_no_fallback.fixedrho_joint_fallback_used
    @test isapprox(via_joint_no_fallback.pressure, 21.62219502138967; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_no_fallback.rho_norm, 0.20000000005601337; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_no_fallback.entropy, 0.1988982719799636; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_no_fallback.energy, -21.36959370985176; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_no_fallback.residual_norm, 6.385724668330235e-11; rtol=1e-4, atol=1e-12)

    via_joint_with_fallback = Models.solve_constraint(
        model,
        mode,
        T_fm;
        problem_spec=spec,
        fixedrho_joint_solve=true,
        seed_guess=seed,
        seed_candidates=(seed,),
        p_num=8,
        t_num=4,
        nlsolve_method=:newton,
        trust_region_fallback=true,
        fallback_method=:trust_region,
        residual_norm_max=1e-6,
    )

    @test via_joint_with_fallback.converged
    @test via_joint_with_fallback.fixedrho_joint_solve_active
    @test haskey(via_joint_with_fallback, :fixedrho_joint_selected_method)
    @test via_joint_with_fallback.fixedrho_joint_selected_method in (:newton, :trust_region)
    @test haskey(via_joint_with_fallback, :fixedrho_joint_fallback_used)
    @test (via_joint_with_fallback.fixedrho_joint_selected_method == :trust_region) == via_joint_with_fallback.fixedrho_joint_fallback_used
    @test isapprox(via_joint_with_fallback.pressure, 21.62219502138967; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_with_fallback.rho_norm, 0.20000000005601337; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_with_fallback.entropy, 0.1988982719799636; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_with_fallback.energy, -21.36959370985176; rtol=1e-6, atol=1e-8)
    @test isapprox(via_joint_with_fallback.residual_norm, 6.385724668330235e-11; rtol=1e-4, atol=1e-12)
end
