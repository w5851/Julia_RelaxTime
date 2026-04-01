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
            status=:active,
            removal_wave=:C,
            removal_threshold="all call sites route via solve_constraint",
        ),
        (
            old_entry="Models.solve_fixedrho_constraint",
            new_entry="Models.solve_constraint(model, FixedRho(...), T)",
            status=:active,
            removal_wave=:C,
            removal_threshold="all call sites route via solve_constraint",
        ),
        (
            old_entry="Models.solve_fixedasymrho_constraint",
            new_entry="Models.solve_constraint(model, FixedAsymmetricRho(...), T)",
            status=:active,
            removal_wave=:C,
            removal_threshold="all call sites route via solve_constraint",
        ),
        (
            old_entry="Models.solve_fixedentropy_constraint",
            new_entry="Models.solve_constraint(model, FixedEntropy(...), T)",
            status=:active,
            removal_wave=:C,
            removal_threshold="all call sites route via solve_constraint",
        ),
        (
            old_entry="Models.solve_fixedsigma_constraint",
            new_entry="Models.solve_constraint(model, FixedSigma(...), T)",
            status=:active,
            removal_wave=:C,
            removal_threshold="all call sites route via solve_constraint",
        ),
    ]
end

function solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; μ_fm::Real, kwargs...)
    # COMPAT: unified primary entrypoint; legacy fixedmu helper retained for migration window.
    return solve_fixedmu_constraint(model, T_fm, μ_fm; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; kwargs...)
    # COMPAT: unified primary entrypoint; legacy fixedrho helper retained for migration window.
    return solve_fixedrho_constraint(model, T_fm, mode.rho_target; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; kwargs...)
    # COMPAT: unified primary entrypoint; legacy fixedentropy helper retained for migration window.
    return solve_fixedentropy_constraint(model, T_fm, mode.s_target; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; kwargs...)
    # COMPAT: unified primary entrypoint; legacy fixedsigma helper retained for migration window.
    return solve_fixedsigma_constraint(model, T_fm, mode.sigma_target; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; kwargs...)
    # COMPAT: unified primary entrypoint; legacy fixedasymrho helper retained for migration window.
    return solve_fixedasymrho_constraint(model, T_fm, mode.rho_target, mode.ud_ratio_target, mode.s_target; kwargs...)
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
