"""
    solve_constraint(model, mode, T_fm; kwargs...)

统一约束求解入口。根据 `mode` 的类型分发到 models 域内核 `solve_fixed*` 实现。
"""
const SolverResult = ImplicitSolver.SolverResult

@inline function solver_migration_map()
    return [
        (
            old_entry="Models.solve_fixedmu_constraint",
            new_entry="Models.solve_constraint(model, FixedMu(), T; μ_fm=...)",
            status=:hard_deprecated,
            removal_wave=:D,
            removal_threshold="hard-deprecate-first; remove only after all external call sites migrate and Wave-D parity/smoke/regression checks stay green",
        ),
        (
            old_entry="Models.solve_fixedrho_constraint",
            new_entry="Models.solve_constraint(model, FixedRho(...), T)",
            status=:hard_deprecated,
            removal_wave=:D,
            removal_threshold="hard-deprecate-first; remove only after all external call sites migrate and Wave-D parity/smoke/regression checks stay green",
        ),
        (
            old_entry="Models.solve_fixedasymrho_constraint",
            new_entry="Models.solve_constraint(model, FixedAsymmetricRho(...), T)",
            status=:hard_deprecated,
            removal_wave=:D,
            removal_threshold="hard-deprecate-first; remove only after all external call sites migrate and Wave-D parity/smoke/regression checks stay green",
        ),
        (
            old_entry="Models.solve_fixedentropy_constraint",
            new_entry="Models.solve_constraint(model, FixedEntropy(...), T)",
            status=:hard_deprecated,
            removal_wave=:D,
            removal_threshold="hard-deprecate-first; remove only after all external call sites migrate and Wave-D parity/smoke/regression checks stay green",
        ),
        (
            old_entry="Models.solve_fixedsigma_constraint",
            new_entry="Models.solve_constraint(model, FixedSigma(...), T)",
            status=:hard_deprecated,
            removal_wave=:D,
            removal_threshold="hard-deprecate-first; remove only after all external call sites migrate and Wave-D parity/smoke/regression checks stay green",
        ),
    ]
end

function solver_migration_status(old_entry::AbstractString)
    for entry in solver_migration_map()
        if entry.old_entry == old_entry
            return (
                old_entry=entry.old_entry,
                unified_entry=entry.new_entry,
                status=entry.status,
                removal_wave=entry.removal_wave,
                removal_threshold=entry.removal_threshold,
                route=:compat_shim,
                deprecation_ready=false,
            )
        end
    end
    throw(ArgumentError("unknown solver compatibility entry: $(old_entry)"))
end

function solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; μ_fm::Real, kwargs...)
    return solve_fixedmu_constraint(model, T_fm, μ_fm; __allow_hard_deprecated_internal=true, kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; kwargs...)
    return solve_fixedrho_constraint(model, T_fm, mode.rho_target; __allow_hard_deprecated_internal=true, kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; kwargs...)
    return solve_fixedentropy_constraint(model, T_fm, mode.s_target; __allow_hard_deprecated_internal=true, kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; kwargs...)
    return solve_fixedsigma_constraint(model, T_fm, mode.sigma_target; __allow_hard_deprecated_internal=true, kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; kwargs...)
    return solve_fixedasymrho_constraint(model, T_fm, mode.rho_target, mode.ud_ratio_target, mode.s_target; __allow_hard_deprecated_internal=true, kwargs...)
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
