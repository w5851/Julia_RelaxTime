module ScanResultFinalize

using StaticArrays

using ..Models: SolverResult

export finalize_solver_result
export promote_near_converged
export is_success
export refine_near_converged

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
        result.xi,
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
