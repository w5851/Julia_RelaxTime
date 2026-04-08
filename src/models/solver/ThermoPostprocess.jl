"""ThermoPostprocess

统一约束求解后处理：
- 从 solver solution 计算热力学量；
- 计算统一 residual_norm；
- 组装标准候选结果字段。
"""

@inline function _thermo_quantities_finite(thermo)::Bool
    return isfinite(thermo.omega) && isfinite(thermo.pressure) &&
           isfinite(thermo.rho_norm) && isfinite(thermo.entropy) && isfinite(thermo.energy)
end

function compute_thermo_from_solution(
    model::AbstractQCDModel,
    solution::AbstractVector,
    T_fm::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
    rho0_scale::Union{Nothing, Real}=rho0,
    state_n::Int=5,
    mu_n::Int=3,
)
    x_state, mu_vec = _unpack_solution(solution; state_n=state_n, mu_n=mu_n)
    thermo = _compute_mode_thermo_quantities(
        model,
        x_state,
        T_fm,
        mu_vec;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        rho0_scale=rho0_scale,
    )
    return (; x_state=x_state, mu_vec=mu_vec, thermo...)
end

function compute_residual_norm_from_solution(
    model::AbstractQCDModel,
    solution::AbstractVector,
    T_fm::Real,
    components...;
    xi::Real,
    p_num::Int,
    t_num::Int,
    state_n::Int=5,
    mu_n::Int=3,
    residual_fn::Union{Nothing, Function}=nothing,
)
    if residual_fn !== nothing
        residual_vec = zeros(Float64, length(solution))
        residual_fn(residual_vec, solution)
        return sqrt(sum(abs2, residual_vec))
    end

    x_state, mu_vec = _unpack_solution(solution; state_n=state_n, mu_n=mu_n)
    return _compose_mode_residual_norm(
        model,
        x_state,
        mu_vec,
        T_fm,
        components...;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )
end

function build_solver_candidate(
    solution::AbstractVector,
    thermo,
    residual_norm::Real;
    converged::Bool,
    iterations::Integer,
    residual_norm_max::Real,
    extra_fields::NamedTuple=NamedTuple(),
)
    candidate = (
        converged=Bool(converged),
        solution=Float64.(solution),
        x_state=thermo.x_state,
        mu_vec=thermo.mu_vec,
        omega=Float64(thermo.omega),
        pressure=Float64(thermo.pressure),
        rho_norm=Float64(thermo.rho_norm),
        entropy=Float64(thermo.entropy),
        energy=Float64(thermo.energy),
        masses=thermo.masses,
        iterations=Int(iterations),
        residual_norm=Float64(residual_norm),
        residual_norm_max=Float64(residual_norm_max),
    )
    return merge(candidate, extra_fields)
end
