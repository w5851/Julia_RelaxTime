"""
    Conditions

PNJL 求解器条件函数定义。

## 核心思想
将能隙方程的条件定义与求解器分离，支持多种约束模式。

## 主要函数
- `gap_conditions`: 核心能隙方程（5维，所有模式共用）
- `build_conditions`: 根据模式构建完整条件函数
- `build_residual!`: 构建 NLsolve 兼容的残差函数
"""
module Conditions

using StaticArrays
using ForwardDiff
import ..Models

# 从 Models 域导入
using ..Models: ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma
using ..Models: ModelStateSchema, state_dim, state_var_dim, mu_var_dim, schema_for_model, state_view, mu_view
const cached_nodes = Models.cached_nodes
using Main.Constants_PNJL: ρ0_inv_fm3
const ρ0 = ρ0_inv_fm3

using ..Models: AbstractQCDModel, AbstractPNJLModel, PNJLModel, PNJLMagneticModel, RPNJLModel
using ..Models: model_pressure, model_rho, model_thermo, calculate_mass_vec
using ..Models: AbstractNJLModel, NJL2Model, gap_residual

export gap_conditions, gap_core_residual!, build_conditions, build_residual!
export GapParams
export explicit_residual, explicit_residual!

@inline _model_kind_symbol(::AbstractPNJLModel) = :PNJL
@inline _model_kind_symbol(::PNJLModel) = :PNJL
@inline _model_kind_symbol(::PNJLMagneticModel) = :PNJL
@inline _model_kind_symbol(::RPNJLModel) = :RPNJL
@inline _model_kind_symbol(::AbstractNJLModel) = :NJL
@inline _model_kind_symbol(::NJL2Model) = :NJL2
@inline _model_kind_symbol(::AbstractQCDModel) = :PNJL

@inline function _get_model(model_kind::Symbol)
    if model_kind === :PNJL || model_kind === :RPNJL || model_kind === :NJL || model_kind === :NJL2
        return Models.get_cached_model(model_kind)
    end
    error("Unsupported model kind in Conditions: $(model_kind)")
end

# ============================================================================
# 参数结构
# ============================================================================

"""
    GapParams

能隙方程求解所需的参数集合。

# 字段
- `T_fm`: 温度 (fm⁻¹)
- `thermal_nodes`: 积分节点
- `xi`: 各向异性参数
"""
struct GapParams{TT, TN, TX}
    T_fm::TT
    thermal_nodes::TN
    xi::TX

    p_num::Int
    t_num::Int
    model_kind::Symbol
end

@inline function GapParams(
    T_fm::TT,
    thermal_nodes::TN,
    xi::TX;
    p_num::Int=size(thermal_nodes[1], 1),
    t_num::Int=size(thermal_nodes[1], 2),
    model_kind::Symbol=:PNJL,
) where {TT, TN, TX}
    return GapParams(T_fm, thermal_nodes, xi, p_num, t_num, model_kind)
end

@inline function _rho_vec(x_state, mu_vec, T_fm, params::GapParams)
    return model_rho(
        _get_model(params.model_kind),
        x_state,
        mu_vec,
        T_fm;
        p_num=params.p_num,
        t_num=params.t_num,
        xi=params.xi,
    )
end

@inline function _thermo_tuple(x_state, mu_vec, T_fm, params::GapParams)
    return model_thermo(
        _get_model(params.model_kind),
        x_state,
        mu_vec,
        T_fm;
        p_num=params.p_num,
        t_num=params.t_num,
        xi=params.xi,
    )
end

@inline function _mass_vec(x_state, params::GapParams)
    phi = SVector{3}(x_state[1], x_state[2], x_state[3])
    return calculate_mass_vec(_get_model(params.model_kind), phi)
end

@inline function _safe_density_ratio(num, den; eps::Float64=1e-12)
    abs(den) >= eps && return num / den
    den_safe = den >= 0 ? eps : -eps
    return num / den_safe
end

# ============================================================================
# 核心能隙条件（所有模式共用）
# ============================================================================

"""
    gap_conditions(x_state, mu_vec, params) -> SVector{5}

计算 5 个能隙方程的残差：∂Ω/∂φ_u, ∂Ω/∂φ_d, ∂Ω/∂φ_s, ∂Ω/∂Φ, ∂Ω/∂Φ̄

这是所有求解模式的核心，通过对压强求梯度实现。

# 参数
- `x_state`: 状态向量 [φ_u, φ_d, φ_s, Φ, Φ̄]
- `mu_vec`: 化学势向量 [μ_u, μ_d, μ_s]
- `params`: GapParams 参数结构

# 返回
- SVector{5}: 能隙方程残差（平衡态时为零）
"""
function gap_conditions(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, params::GapParams) where {TF, TM}
    model = _get_model(params.model_kind)
    pressure_fn = y -> begin
        eltp = typeof(y[1])
        y_s = SVector{5, eltp}(Tuple(y))
        model_pressure(
            model,
            y_s,
            mu_vec,
            params.T_fm;
            p_num=params.p_num,
            t_num=params.t_num,
            xi=params.xi,
        )
    end
    grad = ForwardDiff.gradient(pressure_fn, x_state)
    grad_type = typeof(grad[1])
    return SVector{5, grad_type}(Tuple(grad))
end

@inline function _gap_conditions_with_model(
    model::AbstractQCDModel,
    x_state::SVector{5, TF},
    mu_vec::AbstractVector{TM},
    params::GapParams,
) where {TF, TM}
    pressure_fn = y -> begin
        eltp = typeof(y[1])
        y_s = SVector{5, eltp}(Tuple(y))
        model_pressure(
            model,
            y_s,
            mu_vec,
            params.T_fm;
            p_num=params.p_num,
            t_num=params.t_num,
            xi=params.xi,
        )
    end
    grad = ForwardDiff.gradient(pressure_fn, x_state)
    grad_type = typeof(grad[1])
    return SVector{5, grad_type}(Tuple(grad))
end

@inline function gap_core_residual!(F::AbstractVector, x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, params::GapParams) where {TF, TM}
    length(F) == 5 || throw(ArgumentError("gap_core_residual! expects output length 5, got $(length(F))"))
    core = gap_conditions(x_state, mu_vec, params)
    @inbounds for i in 1:5
        F[i] = core[i]
    end
    return F
end

@inline function gap_core_residual!(
    F::AbstractVector,
    model::AbstractQCDModel,
    x_state::AbstractVector,
    mu_vec::AbstractVector,
    params::GapParams,
)
    length(mu_vec) == 3 || throw(ArgumentError("gap_core_residual! expects mu_vec length 3, got $(length(mu_vec))"))

    nx = length(x_state)
    length(F) == nx || throw(ArgumentError("gap_core_residual! expects output length $nx, got $(length(F))"))

    if nx == 5
        Tx = typeof(x_state[1])
        x5 = SVector{5, Tx}(x_state[1], x_state[2], x_state[3], x_state[4], x_state[5])
        core = _gap_conditions_with_model(model, x5, mu_vec, params)
        @inbounds for i in 1:5
            F[i] = core[i]
        end
        return F
    elseif nx == 2 || nx == 3
        residual = gap_residual(
            model,
            x_state,
            params.T_fm,
            mu_vec;
            p_num=params.p_num,
            t_num=params.t_num,
            xi=params.xi,
        )
        @inbounds for i in 1:nx
            F[i] = residual[i]
        end
        return F
    end

    throw(ArgumentError("gap_core_residual! only supports state_dim in (2,3,5), got $nx"))
end

@inline function gap_core_residual!(F::AbstractVector, x_state::AbstractVector, mu_vec::AbstractVector, params::GapParams)
    return gap_core_residual!(F, _get_model(params.model_kind), x_state, mu_vec, params)
end

# ============================================================================
# 模式特定条件构建
# ============================================================================

"""
    build_conditions(mode::FixedMu, params::GapParams) -> Function

构建固定化学势模式的条件函数。

返回函数签名：(θ, x) -> residual
- θ = [T, μ]（参数）
- x = [φ_u, φ_d, φ_s, Φ, Φ̄]（状态变量）
"""
function build_conditions(::FixedMu, params::GapParams)
    schema = schema_for_model(params.model_kind)
    return (θ, x) -> begin
        T_fm = θ[1]
        μ_fm = θ[2]
        length(x) == state_dim(schema) || throw(ArgumentError("FixedMu expects x length $(state_dim(schema)), got $(length(x))"))
        mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
        x_state = SVector{5}(Tuple(state_view(schema, x)))
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)
        return Vector(gap_conditions(x_state, mu_vec, local_params))
    end
end

"""
    build_conditions(mode::FixedRho, params::GapParams) -> Function

构建固定密度模式的条件函数。

返回函数签名：(θ, x) -> residual
- θ = [T]（参数）
- x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（状态变量）
"""
function build_conditions(mode::FixedRho, params::GapParams)
    schema = schema_for_model(params.model_kind)
    return build_conditions(mode, params, schema; mu_dim=mu_var_dim(mode))
end

@inline function _extract_state_mu(schema::ModelStateSchema, x::AbstractVector; mu_dim::Int=3)
    state_slice = state_view(schema, x)
    mu_slice = mu_view(schema, x; mu_dim=mu_dim)
    return state_slice, mu_slice
end

@inline function _validate_schema_mode_dims(mode::ConstraintMode, schema::ModelStateSchema, mu_dim::Int)
    expected_state_dim = state_var_dim(mode)
    schema_state_dim = state_dim(schema)
    schema_state_dim == expected_state_dim || throw(ArgumentError("schema state_dim mismatch for $(typeof(mode)): expected $expected_state_dim, got $schema_state_dim"))

    expected_mu_dim = mu_var_dim(mode)
    mu_dim == expected_mu_dim || throw(ArgumentError("schema mu_dim mismatch for $(typeof(mode)): expected $expected_mu_dim, got $mu_dim"))
    return nothing
end

@inline function _gap_conditions_dynamic(
    mode::ConstraintMode,
    schema::ModelStateSchema,
    x_state::AbstractVector,
    mu_vec::AbstractVector,
    params::GapParams;
    mu_dim::Int,
)
    _validate_schema_mode_dims(mode, schema, mu_dim)
    state_n = length(x_state)
    mu_n = length(mu_vec)
    if (state_n == 2 || state_n == 3 || state_n == 5) && mu_n == 3
        Tout = promote_type(eltype(x_state), eltype(mu_vec), typeof(params.T_fm))
        out = Vector{Tout}(undef, state_n)
        gap_core_residual!(out, x_state, mu_vec, params)
        return out
    end
    throw(ArgumentError("schema/mode dimensions currently support state_dim in (2,3,5) and mu_dim=3, got state_dim=$state_n, mu_dim=$mu_n"))
end

"""
    build_conditions(mode::FixedRho, params::GapParams, schema::ModelStateSchema; mu_dim=3) -> Function

基于 schema 的条件函数构建路径（迁移期接口）。
"""
function build_conditions(mode::FixedRho, params::GapParams, schema::ModelStateSchema; mu_dim::Int=3)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state, mu_vec = _extract_state_mu(schema, x; mu_dim=mu_dim)
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)

        gap = _gap_conditions_dynamic(mode, schema, x_state, mu_vec, local_params; mu_dim=mu_dim)
        mu_eq1 = mu_vec[1] - mu_vec[2]
        mu_eq2 = mu_vec[2] - mu_vec[3]
        rho_vec = _rho_vec(x_state, mu_vec, T_fm, local_params)
        rho_constraint = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target

        return [gap..., mu_eq1, mu_eq2, rho_constraint]
    end
end

"""
    build_conditions(mode::FixedAsymmetricRho, params::GapParams) -> Function

构建固定非对称约束模式的条件函数。

返回函数签名：(θ, x) -> residual
- θ = [T]（参数）
- x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（状态变量）
"""
function build_conditions(mode::FixedAsymmetricRho, params::GapParams)
    schema = schema_for_model(params.model_kind)
    return build_conditions(mode, params, schema; mu_dim=mu_var_dim(mode))
end

function build_conditions(mode::FixedAsymmetricRho, params::GapParams, schema::ModelStateSchema; mu_dim::Int=3)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state, mu_vec = _extract_state_mu(schema, x; mu_dim=mu_dim)
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)

        gap = _gap_conditions_dynamic(mode, schema, x_state, mu_vec, local_params; mu_dim=mu_dim)

        rho_vec = _rho_vec(x_state, mu_vec, T_fm, local_params)
        rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]

        nB_constraint = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target
        ud_ratio_constraint = _safe_density_ratio(rho_u, rho_d) - mode.ud_ratio_target
        s_constraint = rho_s - mode.s_target

        return [gap..., nB_constraint, ud_ratio_constraint, s_constraint]
    end
end

"""
    build_conditions(mode::FixedEntropy, params::GapParams) -> Function

构建固定熵密度模式的条件函数。
"""
function build_conditions(mode::FixedEntropy, params::GapParams)
    schema = schema_for_model(params.model_kind)
    return build_conditions(mode, params, schema; mu_dim=mu_var_dim(mode))
end

function build_conditions(mode::FixedEntropy, params::GapParams, schema::ModelStateSchema; mu_dim::Int=3)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state, mu_vec = _extract_state_mu(schema, x; mu_dim=mu_dim)
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)

        gap = _gap_conditions_dynamic(mode, schema, x_state, mu_vec, local_params; mu_dim=mu_dim)
        mu_eq1 = mu_vec[1] - mu_vec[2]
        mu_eq2 = mu_vec[2] - mu_vec[3]

        _, _, entropy, _ = _thermo_tuple(x_state, mu_vec, T_fm, local_params)
        s_constraint = entropy - mode.s_target

        return [gap..., mu_eq1, mu_eq2, s_constraint]
    end
end

"""
    build_conditions(mode::FixedSigma, params::GapParams) -> Function

构建固定比熵模式的条件函数。
"""
function build_conditions(mode::FixedSigma, params::GapParams)
    schema = schema_for_model(params.model_kind)
    return build_conditions(mode, params, schema; mu_dim=mu_var_dim(mode))
end

function build_conditions(mode::FixedSigma, params::GapParams, schema::ModelStateSchema; mu_dim::Int=3)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state, mu_vec = _extract_state_mu(schema, x; mu_dim=mu_dim)
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)

        gap = _gap_conditions_dynamic(mode, schema, x_state, mu_vec, local_params; mu_dim=mu_dim)
        mu_eq1 = mu_vec[1] - mu_vec[2]
        mu_eq2 = mu_vec[2] - mu_vec[3]

        _, rho_norm, entropy, _ = _thermo_tuple(x_state, mu_vec, T_fm, local_params)
        n_B = rho_norm * ρ0
        sigma = n_B > 1e-12 ? entropy / n_B : 0.0
        sigma_constraint = sigma - mode.sigma_target

        return [gap..., mu_eq1, mu_eq2, sigma_constraint]
    end
end

# ============================================================================
# NLsolve 兼容的残差函数构建
# ============================================================================

"""
    build_residual!(mode::FixedMu, mu_vec, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定化学势模式）。

返回函数签名：(F, x) -> nothing（原地修改 F）
"""
function build_residual!(::FixedMu, mu_vec::SVector{3}, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x))
        gap_core_residual!(F, x_state, mu_vec, params)
        return nothing
    end
end

"""
    build_residual!(mode::FixedRho, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定密度模式）。
"""
function build_residual!(mode::FixedRho, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])
        
        gap_core_residual!(@view(F[1:5]), x_state, mu_state, params)
        
        # 化学势相等
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
        
        # 密度约束
        rho = _rho_vec(x_state, mu_state, params.T_fm, params)
        F[8] = sum(rho) / (3.0 * ρ0) - mode.rho_target
        
        return nothing
    end
end

"""
    build_residual!(mode::FixedAsymmetricRho, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定非对称约束模式）。
"""
function build_residual!(mode::FixedAsymmetricRho, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])

        gap_core_residual!(@view(F[1:5]), x_state, mu_state, params)

        rho = _rho_vec(x_state, mu_state, params.T_fm, params)
        rho_u, rho_d, rho_s = rho[1], rho[2], rho[3]

        F[6] = sum(rho) / (3.0 * ρ0) - mode.rho_target
        F[7] = _safe_density_ratio(rho_u, rho_d) - mode.ud_ratio_target
        F[8] = rho_s - mode.s_target

        return nothing
    end
end

"""
    build_residual!(mode::FixedEntropy, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定熵密度模式）。
"""
function build_residual!(mode::FixedEntropy, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])
        
        gap_core_residual!(@view(F[1:5]), x_state, mu_state, params)
        
        # 化学势相等
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
        
        # 熵密度约束
        _, _, entropy, _ = _thermo_tuple(x_state, mu_state, params.T_fm, params)
        F[8] = entropy - mode.s_target
        
        return nothing
    end
end

"""
    build_residual!(mode::FixedSigma, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定比熵模式）。
"""
function build_residual!(mode::FixedSigma, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])
        
        gap_core_residual!(@view(F[1:5]), x_state, mu_state, params)
        
        # 化学势相等
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
        
        # 比熵约束
        _, rho_norm, entropy, _ = _thermo_tuple(x_state, mu_state, params.T_fm, params)
        n_B = rho_norm * ρ0
        sigma = n_B > 1e-12 ? entropy / n_B : 0.0
        F[8] = sigma - mode.sigma_target
        
        return nothing
    end
end

function build_residual!(model::AbstractQCDModel, mode::FixedMu, T_fm::Real, μ_fm::Real;
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
)
    thermal_nodes = cached_nodes(p_num, t_num; p_max_inv_fm=Models.thermal_p_max_inv_fm(model))
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num,
        t_num=t_num,
        model_kind=_model_kind_symbol(model),
    )
    mu_vec = SVector{3}(Float64(μ_fm), Float64(μ_fm), Float64(μ_fm))
    return build_residual!(mode, mu_vec, params)
end

function build_residual!(model::AbstractQCDModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real;
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
)
    thermal_nodes = cached_nodes(p_num, t_num; p_max_inv_fm=Models.thermal_p_max_inv_fm(model))
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num,
        t_num=t_num,
        model_kind=_model_kind_symbol(model),
    )
    return build_residual!(mode, params)
end

@inline function _local_gap_params(T_fm, params::GapParams)
    return GapParams(T_fm, params.thermal_nodes, params.xi;
        p_num=params.p_num,
        t_num=params.t_num,
        model_kind=params.model_kind,
    )
end

@inline function explicit_residual(mode::FixedMu, x::AbstractVector, θ::AbstractVector, params::GapParams)
    out = zeros(promote_type(eltype(x), eltype(θ)), state_dim(mode))
    explicit_residual!(out, x, θ, params, mode)
    return out
end

@inline function explicit_residual(mode::FixedRho, x::AbstractVector, θ::AbstractVector, params::GapParams)
    out = zeros(promote_type(eltype(x), eltype(θ)), state_dim(mode))
    explicit_residual!(out, x, θ, params, mode)
    return out
end

@inline function explicit_residual(mode::FixedAsymmetricRho, x::AbstractVector, θ::AbstractVector, params::GapParams)
    out = zeros(promote_type(eltype(x), eltype(θ)), state_dim(mode))
    explicit_residual!(out, x, θ, params, mode)
    return out
end

@inline function explicit_residual(mode::FixedEntropy, x::AbstractVector, θ::AbstractVector, params::GapParams)
    out = zeros(promote_type(eltype(x), eltype(θ)), state_dim(mode))
    explicit_residual!(out, x, θ, params, mode)
    return out
end

@inline function explicit_residual(mode::FixedSigma, x::AbstractVector, θ::AbstractVector, params::GapParams)
    out = zeros(promote_type(eltype(x), eltype(θ)), state_dim(mode))
    explicit_residual!(out, x, θ, params, mode)
    return out
end

@inline function explicit_residual!(F::AbstractVector, x::AbstractVector, θ::AbstractVector, params::GapParams, mode::FixedMu)
    T_fm = θ[1]
    μ_fm = θ[2]
    x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    gap = gap_conditions(x_state, mu_vec, _local_gap_params(T_fm, params))

    @inbounds for i in 1:5
        F[i] = gap[i]
    end
    return nothing
end

@inline function explicit_residual!(F::AbstractVector, x::AbstractVector, θ::AbstractVector, params::GapParams, mode::FixedRho)
    T_fm = θ[1]
    x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
    mu_vec = SVector{3}(x[6], x[7], x[8])
    local_params = _local_gap_params(T_fm, params)

    gap = gap_conditions(x_state, mu_vec, local_params)
    @inbounds for i in 1:5
        F[i] = gap[i]
    end
    @inbounds begin
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
    end
    rho_vec = _rho_vec(x_state, mu_vec, T_fm, local_params)
    F[8] = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target
    return nothing
end

@inline function explicit_residual!(F::AbstractVector, x::AbstractVector, θ::AbstractVector, params::GapParams, mode::FixedAsymmetricRho)
    T_fm = θ[1]
    x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
    mu_vec = SVector{3}(x[6], x[7], x[8])
    local_params = _local_gap_params(T_fm, params)

    gap = gap_conditions(x_state, mu_vec, local_params)
    @inbounds for i in 1:5
        F[i] = gap[i]
    end
    rho_vec = _rho_vec(x_state, mu_vec, T_fm, local_params)
    rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]
    F[6] = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target
    F[7] = _safe_density_ratio(rho_u, rho_d) - mode.ud_ratio_target
    F[8] = rho_s - mode.s_target
    return nothing
end

@inline function explicit_residual!(F::AbstractVector, x::AbstractVector, θ::AbstractVector, params::GapParams, mode::FixedEntropy)
    T_fm = θ[1]
    x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
    mu_vec = SVector{3}(x[6], x[7], x[8])
    local_params = _local_gap_params(T_fm, params)

    gap = gap_conditions(x_state, mu_vec, local_params)
    @inbounds for i in 1:5
        F[i] = gap[i]
    end
    @inbounds begin
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
    end
    _, _, entropy, _ = _thermo_tuple(x_state, mu_vec, T_fm, local_params)
    F[8] = entropy - mode.s_target
    return nothing
end

@inline function explicit_residual!(F::AbstractVector, x::AbstractVector, θ::AbstractVector, params::GapParams, mode::FixedSigma)
    T_fm = θ[1]
    x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
    mu_vec = SVector{3}(x[6], x[7], x[8])
    local_params = _local_gap_params(T_fm, params)

    gap = gap_conditions(x_state, mu_vec, local_params)
    @inbounds for i in 1:5
        F[i] = gap[i]
    end
    @inbounds begin
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
    end
    _, rho_norm, entropy, _ = _thermo_tuple(x_state, mu_vec, T_fm, local_params)
    n_B = rho_norm * ρ0
    sigma = n_B > 1e-12 ? entropy / n_B : 0.0
    F[8] = sigma - mode.sigma_target
    return nothing
end

@inline function explicit_residual!(F::AbstractVector, x::AbstractVector, θ::AbstractVector, params::GapParams, mode::ConstraintMode)
    throw(ArgumentError("unsupported mode for explicit_residual!: $(typeof(mode))"))
    return nothing
end

end # module Conditions
