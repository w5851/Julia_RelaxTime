module ScanResultFinalize

using StaticArrays

using ..Models: SolverResult, SolverConvergence, solver_convergence_view

export finalize_solver_result
export promote_near_converged
export is_success
export refine_near_converged
export scan_acceptance_diagnostics

"""finalize_solver_result(result, T_fm, xi; solver_backend, p_num, t_num, model_kind=:PNJL)

统一处理 scan 单点求解后的输出（单路径 models 实现）。
"""
function finalize_solver_result(result::SolverResult, T_fm, xi;
    solver_backend::Symbol,
    p_num::Int,
    t_num::Int,
    model_kind::Symbol=:PNJL,
)
    return result
end

"""promote_near_converged(result; acceptable_residual) -> (result, message)

Scan 辅助：当结果未标记为收敛但残差足够小（近似收敛）时，强制把 `converged=true`。
用于保持 scan 输出的“成功/失败”口径稳定，但不改变数值内容。
"""
function promote_near_converged(result; acceptable_residual::Real)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end

    residual = result.residual_norm
    if !isfinite(residual) || residual > acceptable_residual
        return result, ""
    end

    promoted = SolverResult(
        result.mode,
        true,
        copy(result.solution),
        result.x_state,
        result.mu_vec,
        result.omega,
        result.pressure,
        result.rho_norm,
        result.entropy,
        result.energy,
        result.masses,
        result.iterations,
        residual,
        result.xi;
        convergence=SolverConvergence(; solver_convergence_view(result)..., scan_promoted=true),
    )
    msg = "force-marked converged (residual $(residual))"
    return promoted, msg
end

"""is_success(result; acceptable_residual) -> Bool

Scan 判定：
- `result === nothing` -> false
- `result.converged` -> true
- 否则：`isfinite(residual_norm) && residual_norm <= acceptable_residual`
"""
function is_success(result; acceptable_residual::Real)
    result === nothing && return false
    result.converged && return true
    residual = result.residual_norm
    return isfinite(residual) && residual <= acceptable_residual
end

"""Explain an existing scan decision without changing its acceptance policy."""
function scan_acceptance_diagnostics(result; acceptable_residual::Real,
    accepted::Bool=is_success(result; acceptable_residual), physical_valid=nothing)
    c = solver_convergence_view(result)
    physical = c.physical_valid === nothing ? physical_valid : c.physical_valid
    scan_physical = physical_valid === nothing ? physical : physical_valid
    reason = if result === nothing
        :no_result
    elseif !accepted
        scan_physical === false ? :physical_invalid :
            !isfinite(result.residual_norm) ? :nonfinite_residual :
            result.residual_norm > acceptable_residual ? :residual_rejected : :scan_rejected
    elseif c.scan_promoted || !result.converged
        :near_residual
    elseif c.solver_converged === true && c.residual_ok === true && physical === true
        :solver_converged
    else
        :legacy_converged
    end
    return (; solver_converged=c.solver_converged, solver_residual_ok=c.residual_ok,
        physical_valid=physical, solver_residual_norm_max=c.residual_norm_max,
        convergence_source=c.source, scan_accepted=accepted, acceptance_reason=reason)
end

"""refine_near_converged(result; acceptable_residual, solve_again) -> (result, message)

Scan 辅助：当结果“接近收敛”（残差小于阈值）但未标记 converged 时，
用当前解作为初值再求解一次，尝试把点拉到真正收敛。

- `solve_again(seed_solution)` 需返回 `(refined_result, message)`，签名与 scan 内 `_solve_point` 一致。
"""
function refine_near_converged(result;
    acceptable_residual::Real,
    solve_again::Function,
)
    result === nothing && return nothing, ""
    result.converged && return result, ""

    residual = result.residual_norm
    if !isfinite(residual) || residual > acceptable_residual
        return result, ""
    end

    refined, msg = solve_again(result.solution)
    if refined !== nothing && refined.converged
        return refined, "refined from near-converged seed"
    end
    return result, msg
end

end # module ScanResultFinalize
