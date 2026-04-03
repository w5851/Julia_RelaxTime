"""
    solve_constraint(model, mode, T_fm; kwargs...)

统一约束求解入口。根据 `mode` 的类型分发到 models 域约束求解内核实现。
"""
const SolverResult = ImplicitSolver.SolverResult

@inline function _solve_with_problem_spec_default(
    model::AbstractQCDModel,
    mode::ConstraintMode,
    T_fm::Real,
    legacy_solver::Function,
    kwargs,
)
    use_problem_spec = get(kwargs, :use_problem_spec, true)
    use_problem_spec isa Bool || throw(ArgumentError("use_problem_spec must be Bool, got $(typeof(use_problem_spec))"))
    allow_legacy_path = get(kwargs, :allow_legacy_path, false)
    allow_legacy_path isa Bool || throw(ArgumentError("allow_legacy_path must be Bool, got $(typeof(allow_legacy_path))"))
    warn_on_legacy_path = get(kwargs, :warn_on_legacy_path, true)
    warn_on_legacy_path isa Bool || throw(ArgumentError("warn_on_legacy_path must be Bool, got $(typeof(warn_on_legacy_path))"))

    if use_problem_spec
        haskey(kwargs, :allow_legacy_path) && throw(ArgumentError("allow_legacy_path is only valid when use_problem_spec=false"))
        haskey(kwargs, :warn_on_legacy_path) && throw(ArgumentError("warn_on_legacy_path is only valid when use_problem_spec=false"))
    else
        haskey(kwargs, :problem_spec) && throw(ArgumentError("problem_spec cannot be combined with use_problem_spec=false"))
    end

    forwarded = (; (k => v for (k, v) in pairs(kwargs) if k != :problem_spec && k != :use_problem_spec && k != :warn_on_legacy_path && k != :allow_legacy_path)...)
    if !use_problem_spec
        if !allow_legacy_path
            throw(ArgumentError("use_problem_spec=false is transitional and now gated; pass allow_legacy_path=true only for explicit compatibility or test scenarios"))
        end
        if warn_on_legacy_path
            @warn "use_problem_spec=false selects transitional legacy solver path; this compatibility route is scheduled for removal"
        end
        return legacy_solver(forwarded)
    end

    spec = get(kwargs, :problem_spec, nothing)
    if spec === nothing
        spec = build_problem_spec(mode)
    else
        spec isa ProblemSpec || throw(ArgumentError("problem_spec must be ProblemSpec or nothing, got $(typeof(spec))"))
    end
    return spec.forward_solve(model, T_fm; forwarded...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; μ_fm::Real, kwargs...)
    return _solve_constraint_fixedmu(model, T_fm, μ_fm; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedRho, T_fm::Real;
    problem_spec::Union{Nothing, ProblemSpec}=nothing,
    kwargs...)
    merged = (; kwargs..., problem_spec=problem_spec)
    return _solve_with_problem_spec_default(
        model,
        mode,
        T_fm,
        local_kwargs -> _solve_constraint_fixedrho(model, T_fm, mode.rho_target; local_kwargs...),
        merged,
    )
end

function solve_constraint(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; kwargs...)
    return _solve_with_problem_spec_default(
        model,
        mode,
        T_fm,
        local_kwargs -> _solve_constraint_fixedentropy(model, T_fm, mode.s_target; local_kwargs...),
        kwargs,
    )
end

function solve_constraint(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; kwargs...)
    return _solve_with_problem_spec_default(
        model,
        mode,
        T_fm,
        local_kwargs -> _solve_constraint_fixedsigma(model, T_fm, mode.sigma_target; local_kwargs...),
        kwargs,
    )
end

function solve_constraint(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; kwargs...)
    return _solve_with_problem_spec_default(
        model,
        mode,
        T_fm,
        local_kwargs -> _solve_constraint_fixedasymrho(model, T_fm, mode.rho_target, mode.ud_ratio_target, mode.s_target; local_kwargs...),
        kwargs,
    )
end

function solve(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    return ImplicitSolver.solve(mode, T_fm, μ_fm; kwargs...)
end

function solve(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real; kwargs...)
    return ImplicitSolver.solve(mode, T_fm; kwargs...)
end

function solve_multi(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    return ImplicitSolver.solve_multi(mode, T_fm, μ_fm; kwargs...)
end

function solve_multi(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
    return ImplicitSolver.solve_multi(mode, T_fm; kwargs...)
end

function solve(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve(model, mode, T_fm, μ_fm; kwargs...)
end

function solve(mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve(model, mode, T_fm; kwargs...)
end

function solve_multi(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    model = create_model(:PNJL)
    if haskey(kwargs, :seeds)
        seeds = kwargs[:seeds]
        candidates = [DefaultSeed(Float64.(seed), Float64.(seed), :hadron) for seed in seeds]
        selector = haskey(kwargs, :selector) ? kwargs[:selector] : SeedStrategies.default_omega_selector
        forward_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k != :seeds && k != :selector)...)
        return solve_multi(model, mode, T_fm, μ_fm;
            seed_strategy=MultiSeed(candidates, selector),
            forward_kwargs...)
    end
    return solve_multi(model, mode, T_fm, μ_fm; kwargs...)
end

function solve_multi(mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve_multi(model, mode, T_fm; kwargs...)
end

@inline create_implicit_solver(; kwargs...) = ImplicitSolver.create_implicit_solver(; kwargs...)
@inline solve_weighted_block_fallback(args...; kwargs...) = ImplicitSolver.solve_weighted_block_fallback(args...; kwargs...)
@inline is_physical_solution(x_state::AbstractVector{<:Real}, masses::AbstractVector{<:Real}; kwargs...) =
    ImplicitSolver.default_is_physical_solution(x_state, masses; kwargs...)
@inline solve_with_derivatives(T_fm::Real, μ_fm::Real; kwargs...) = ImplicitSolver.solve_with_derivatives(T_fm, μ_fm; kwargs...)
@inline solve_with_root_diagnostics(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...) = ImplicitSolver.solve_with_root_diagnostics(mode, T_fm, μ_fm; kwargs...)
@inline solve_with_root_diagnostics(mode::FixedRho, T_fm::Real; kwargs...) = ImplicitSolver.solve_with_root_diagnostics(mode, T_fm; kwargs...)
@inline solve_with_root_diagnostics(mode::FixedAsymmetricRho, T_fm::Real; kwargs...) = ImplicitSolver.solve_with_root_diagnostics(mode, T_fm; kwargs...)
@inline solve_with_root_diagnostics(mode::FixedEntropy, T_fm::Real; kwargs...) = ImplicitSolver.solve_with_root_diagnostics(mode, T_fm; kwargs...)
@inline solve_with_root_diagnostics(mode::FixedSigma, T_fm::Real; kwargs...) = ImplicitSolver.solve_with_root_diagnostics(mode, T_fm; kwargs...)

@inline function solve_with_derivatives(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    _ = model, mode
    return ImplicitSolver.solve_with_derivatives(T_fm, μ_fm; kwargs...)
end
