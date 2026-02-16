"""state.jl

新架构下的平均场状态类型。

动机：
- 明确 `x_state` 每个分量的物理含义，减少索引错误
- 方便后续 PNJL/rPNJL 扩展（例如引入更多状态分量）
- 对 ForwardDiff 等自动微分更友好（稳定的字段访问）

约定：
- `phi` = (φu, φd, φs)
- `Phi`/`PhiBar` 为 Polyakov loop 变量；NJL 可取 1（或任意值但会被忽略）
"""

using StaticArrays

export MeanFieldState
export meanfield_state, state_vector
export normalize_mu_vec

struct MeanFieldState{T}
    phi::SVector{3, T}
    Phi::T
    PhiBar::T
end

function MeanFieldState(phi::SVector{3, T}; Phi::T=one(T), PhiBar::T=one(T)) where {T}
    return MeanFieldState{T}(phi, Phi, PhiBar)
end

function MeanFieldState(x_state::AbstractVector{T}) where {T}
    n = length(x_state)
    if n >= 5
        phi = SVector{3, T}(x_state[1], x_state[2], x_state[3])
        return MeanFieldState{T}(phi, x_state[4], x_state[5])
    elseif n == 3
        phi = SVector{3, T}(x_state[1], x_state[2], x_state[3])
        return MeanFieldState{T}(phi, one(T), one(T))
    end
    throw(ArgumentError("x_state length must be 3 or >= 5, got $n"))
end

function MeanFieldState(x_state::NamedTuple)
    if haskey(x_state, :φ)
        phi = x_state.φ
    elseif haskey(x_state, :phi)
        phi = x_state.phi
    else
        throw(ArgumentError("x_state NamedTuple must contain :φ or :phi"))
    end

    Phi = haskey(x_state, :Φ) ? x_state.Φ : (haskey(x_state, :Phi) ? x_state.Phi : 1.0)
    PhiBar = haskey(x_state, :Φbar) ? x_state.Φbar : (haskey(x_state, :PhiBar) ? x_state.PhiBar : (haskey(x_state, :Phibar) ? x_state.Phibar : 1.0))

    T = promote_type(eltype(phi), typeof(Phi), typeof(PhiBar))
    return MeanFieldState(SVector{3, T}(phi[1], phi[2], phi[3]); Phi=T(Phi), PhiBar=T(PhiBar))
end

"""meanfield_state(x_state) -> MeanFieldState

将多种输入形式（MeanFieldState / AbstractVector / NamedTuple）规范化为 MeanFieldState。
"""
@inline meanfield_state(x_state) = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)

"""state_vector(x_state) -> SVector{5}

将平均场状态转换为 5 维向量表示：`(φu, φd, φs, Φ, Φbar)`。
"""
@inline function state_vector(st::MeanFieldState{T}) where {T}
    return SVector{5, T}(st.phi[1], st.phi[2], st.phi[3], st.Phi, st.PhiBar)
end

@inline state_vector(x_state) = state_vector(meanfield_state(x_state))

"""normalize_mu_vec(mu_vec) -> SVector{3}

化学势输入契约：内部统一为三味 `SVector{3}`。

- `mu_vec::Real`：按对称情形扩展为 (μ, μ, μ)
- `mu_vec::AbstractVector`：必须长度为 3
"""
@inline function normalize_mu_vec(mu::Real)
    # Keep the numeric type (may be Dual) to support AD/implicit differentiation.
    μ = mu
    return SVector{3, typeof(μ)}(μ, μ, μ)
end

@inline normalize_mu_vec(mu_vec::SVector{3}) = mu_vec

@inline function normalize_mu_vec(mu_vec::AbstractVector)
    length(mu_vec) == 3 || throw(ArgumentError("mu_vec must have length 3, got $(length(mu_vec))"))
    μ1, μ2, μ3 = mu_vec[1], mu_vec[2], mu_vec[3]
    T = promote_type(typeof(μ1), typeof(μ2), typeof(μ3))
    return SVector{3, T}(T(μ1), T(μ2), T(μ3))
end
