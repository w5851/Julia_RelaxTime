"""
    solve_constraint(model, mode, T_fm; kwargs...)

统一约束求解入口。根据 `mode` 的类型分发到 models 域内核 `solve_fixed*` 实现。
"""
const SolverResult = ImplicitSolver.SolverResult

function solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; μ_fm::Real, kwargs...)
    return solve_fixedmu_constraint(model, T_fm, μ_fm; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; kwargs...)
    return solve_fixedrho_constraint(model, T_fm, mode.rho_target; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; kwargs...)
    return solve_fixedentropy_constraint(model, T_fm, mode.s_target; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; kwargs...)
    return solve_fixedsigma_constraint(model, T_fm, mode.sigma_target; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; kwargs...)
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

@inline create_implicit_solver(; kwargs...) = ImplicitSolver.create_implicit_solver(; kwargs...)
@inline solve_with_derivatives(T_fm::Real, μ_fm::Real; kwargs...) = ImplicitSolver.solve_with_derivatives(T_fm, μ_fm; kwargs...)
