using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

to_fm_inv(x_mev::Real) = Float64(x_mev) / HBARC_MEV_FM

function _run_point(model, mode, T_mev::Real; mu_mev::Real=0.0)
    T_fm = to_fm_inv(T_mev)
    mu_fm = to_fm_inv(mu_mev)
    theta = mode isa Models.FixedMu ? [T_fm, mu_fm] : [T_fm]

    spec = Models.build_problem_spec(mode)
    seed = Models.get_seed(Models.DefaultSeed(), theta, mode)

    kwargs = Dict{Symbol,Any}(
        :problem_spec => spec,
        :seed_guess => seed,
        :p_num => 8,
        :t_num => 4,
        :residual_norm_max => 1e-6,
    )
    if mode isa Models.FixedMu
        kwargs[Symbol("μ_fm")] = mu_fm
    end

    raw = Models.solve_constraint(model, mode, T_fm; pairs(kwargs)...)
    return (
        converged = Bool(raw.converged),
        residual_norm = Float64(raw.residual_norm),
        selection_reason = Symbol(get(raw, :selection_reason, :none)),
        seed_index = Int(get(raw, :seed_index, -1)),
        failed_constraints = Symbol.(get(raw, :failed_constraints, Symbol[])),
        quality_tag = Symbol(get(raw, :governed_selected_quality, get(raw, :fixedrho_joint_selected_quality, :none))),
        fallback_used = Bool(get(raw, :governed_fallback_used, get(raw, :fixedrho_joint_fallback_used, false))),
    )
end

@testset "solver attempt-engine convergence regression" begin
    model = Models.create_model(:PNJL)

    # Baseline snapshot fields for convergence-chain stability.
    expected = [
        (
            name = :fixedmu_100_250,
            mode = Models.FixedMu(),
            T_mev = 100.0,
            mu_mev = 250.0,
            converged = true,
            residual_upper = 1e-6,
            selection_reason = :pressure_max_under_constraints,
            seed_index = -1,
            failed_constraints = Symbol[],
        ),
        (
            name = :fixedrho_110_0p2,
            mode = Models.FixedRho(0.2),
            T_mev = 110.0,
            mu_mev = 0.0,
            converged = true,
            residual_upper = 1e-6,
            selection_reason = :pressure_max_under_constraints,
            seed_index = -1,
            failed_constraints = Symbol[],
        ),
        (
            name = :fixedasym_110_0p05,
            mode = Models.FixedAsymmetricRho(0.05, 1.0, 0.0),
            T_mev = 110.0,
            mu_mev = 0.0,
            converged = false,
            residual_upper = Inf,
            selection_reason = :no_candidate_passed_constraints,
            seed_index = -1,
            failed_constraints = Symbol[:solver_failed],
        ),
        (
            name = :fixedsigma_120_10,
            mode = Models.FixedSigma(10.0),
            T_mev = 120.0,
            mu_mev = 0.0,
            converged = false,
            residual_upper = Inf,
            selection_reason = :no_candidate_passed_constraints,
            seed_index = -1,
            failed_constraints = Symbol[:solver_failed],
        ),
    ]

    for cfg in expected
        got = _run_point(model, cfg.mode, cfg.T_mev; mu_mev=cfg.mu_mev)
        @test got.converged == cfg.converged
        @test got.selection_reason == cfg.selection_reason
        @test got.seed_index == cfg.seed_index
        @test got.failed_constraints == cfg.failed_constraints
        if cfg.selection_reason == :no_candidate_passed_constraints
            @test got.quality_tag in (:bad, :degraded, :none)
        elseif cfg.selection_reason == :pressure_max_under_constraints
            @test got.quality_tag in (:good, :fallback, :none)
        end
        if got.fallback_used
            @test got.quality_tag in (:fallback, :none)
        end
        if isfinite(cfg.residual_upper)
            @test isfinite(got.residual_norm)
            @test got.residual_norm <= cfg.residual_upper
        else
            @test !isfinite(got.residual_norm)
        end
    end
end
