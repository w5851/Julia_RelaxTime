"""
IntervalQuadratureStrategies.jl

这是一份“无 module”的共享 include 文件：它不会创建新的命名空间，
而是把定义注入到 `include(...)` 发生处的模块中。

设计目标：
- 复用 OneLoopIntegralsAniso 里的通用求积节点生成/变换策略（聚簇/幂次/DE/HYBRID）。
- 不触碰物理相关部分（integrand、找根、区间构建等），避免循环依赖。

何时“需要”封装为 module：
- 当多个模块/脚本频繁 include 导致命名污染、重定义、或加载顺序难以管理时，建议升级为 `module IntervalQuadratureStrategies`。
- 若当前只在少数模块内 include，且你希望最小改动、避免额外命名空间层级，那么维持 include 更省心。

依赖假设（由包含方负责提供）：
- 已经存在 `GaussLegendre` 模块（例如先 include("GaussLegendre.jl")）。
"""

# ============================================================================
# 数值积分策略配置
# ============================================================================

"""积分策略枚举"""
@enum IntegrationStrategy begin
    STRATEGY_QUADGK       # 原始 QuadGK 自适应积分
    STRATEGY_INTERVAL_GL  # 区间分割 + 标准 GL
    STRATEGY_CLUSTER_GL   # 区间分割 + 聚簇 GL (tanh 对称)
    STRATEGY_HYBRID       # 混合策略：根据奇点位置自适应选择变换
end

"""默认积分策略"""
const DEFAULT_STRATEGY = STRATEGY_HYBRID

"""聚簇 GL / 变换默认参数"""
const DEFAULT_CLUSTER_BETA = 8.0
const DEFAULT_CLUSTER_N = 32
const DEFAULT_POWER_ALPHA = 0.35
const DEFAULT_DE_H = 0.15

"""诊断输出结构"""
struct IntegrationDiagnostics
    strategy::IntegrationStrategy
    n_roots::Int
    roots::Vector{Float64}
    n_intervals::Int
    intervals::Vector{Tuple{Float64,Float64}}
    real_part::Float64
    imag_part::Float64
    elapsed_ms::Float64
end

# ============================================================================
# 聚簇 / 变换节点生成（含默认参数的预计算加速）
# ============================================================================

# 标准节点/权重（用于所有变换）
const _STD_32_NODES, _STD_32_WEIGHTS = GaussLegendre.gausslegendre(32)
const _STD_16_NODES, _STD_16_WEIGHTS = GaussLegendre.gausslegendre(16)

# power_left / power_right 预计算（alpha=DEFAULT_POWER_ALPHA）
const _POWER_LEFT_INV_ALPHA = 1.0 / DEFAULT_POWER_ALPHA
const _POWER_LEFT_U_HALF_32 = (_STD_32_NODES .+ 1) ./ 2
const _POWER_LEFT_T_32 = _POWER_LEFT_U_HALF_32 .^ _POWER_LEFT_INV_ALPHA
const _POWER_LEFT_T_PRIME_32 = _POWER_LEFT_INV_ALPHA .* _POWER_LEFT_U_HALF_32 .^ (_POWER_LEFT_INV_ALPHA - 1) ./ 2

const _POWER_LEFT_U_HALF_16 = (_STD_16_NODES .+ 1) ./ 2
const _POWER_LEFT_T_16 = _POWER_LEFT_U_HALF_16 .^ _POWER_LEFT_INV_ALPHA
const _POWER_LEFT_T_PRIME_16 = _POWER_LEFT_INV_ALPHA .* _POWER_LEFT_U_HALF_16 .^ (_POWER_LEFT_INV_ALPHA - 1) ./ 2

const _POWER_RIGHT_T_32 = 1 .- (1 .- _POWER_LEFT_U_HALF_32) .^ _POWER_LEFT_INV_ALPHA
const _POWER_RIGHT_T_PRIME_32 = _POWER_LEFT_INV_ALPHA .* (1 .- _POWER_LEFT_U_HALF_32) .^ (_POWER_LEFT_INV_ALPHA - 1) ./ 2

const _POWER_RIGHT_T_16 = 1 .- (1 .- _POWER_LEFT_U_HALF_16) .^ _POWER_LEFT_INV_ALPHA
const _POWER_RIGHT_T_PRIME_16 = _POWER_LEFT_INV_ALPHA .* (1 .- _POWER_LEFT_U_HALF_16) .^ (_POWER_LEFT_INV_ALPHA - 1) ./ 2

# tanh 聚簇预计算（beta=DEFAULT_CLUSTER_BETA, n=32）
const _TANH_BETA = tanh(DEFAULT_CLUSTER_BETA)
const _TANH_PHI = tanh.(DEFAULT_CLUSTER_BETA .* _STD_32_NODES) ./ _TANH_BETA
const _TANH_PHI_PRIME = DEFAULT_CLUSTER_BETA .* (1 .- tanh.(DEFAULT_CLUSTER_BETA .* _STD_32_NODES).^2) ./ _TANH_BETA

# DE (tanh-sinh) 预计算（h=DEFAULT_DE_H, n=32 -> 33 points）
const _DE_KS = collect(-16:16)
const _DE_TS = _DE_KS .* DEFAULT_DE_H
const _DE_PHI = tanh.(π/2 .* sinh.(_DE_TS))
const _DE_PHI_PRIME = DEFAULT_DE_H .* (π/2) .* cosh.(_DE_TS) .* (sech.(π/2 .* sinh.(_DE_TS))).^2

"""生成聚簇 GL 节点（tanh 映射，默认参数走预计算）"""
function clustered_gl_nodes(a::Float64, b::Float64, n::Int; beta::Float64=DEFAULT_CLUSTER_BETA)
    if n == 32 && beta == DEFAULT_CLUSTER_BETA
        half = (b - a) / 2
        center = (a + b) / 2
        xs = center .+ half .* _TANH_PHI
        wx = _STD_32_WEIGHTS .* half .* _TANH_PHI_PRIME
        return xs, wx
    end

    us, ws = GaussLegendre.gauleg(-1.0, 1.0, n)
    tanh_beta = tanh(beta)
    phi = tanh.(beta .* us) ./ tanh_beta
    phi_prime = beta .* (1 .- tanh.(beta .* us).^2) ./ tanh_beta
    half = (b - a) / 2
    center = (a + b) / 2
    xs = center .+ half .* phi
    wx = ws .* half .* phi_prime
    return xs, wx
end

"""生成单侧幂次聚簇节点（聚簇于左端 a；默认参数走预计算）"""
function power_left_nodes(a::Float64, b::Float64, n::Int; alpha::Float64=DEFAULT_POWER_ALPHA)
    len = b - a
    if alpha == DEFAULT_POWER_ALPHA
        if n == 32
            xs = a .+ len .* _POWER_LEFT_T_32
            wx = _STD_32_WEIGHTS .* len .* _POWER_LEFT_T_PRIME_32
            return xs, wx
        elseif n == 16
            xs = a .+ len .* _POWER_LEFT_T_16
            wx = _STD_16_WEIGHTS .* len .* _POWER_LEFT_T_PRIME_16
            return xs, wx
        end
    end

    us, ws = GaussLegendre.gauleg(-1.0, 1.0, n)
    inv_alpha = 1.0 / alpha
    u_half = (us .+ 1) ./ 2
    t = u_half .^ inv_alpha
    t_prime = inv_alpha .* u_half .^ (inv_alpha - 1) ./ 2
    xs = a .+ len .* t
    wx = ws .* len .* t_prime
    return xs, wx
end

"""生成单侧幂次聚簇节点（聚簇于右端 b；默认参数走预计算）"""
function power_right_nodes(a::Float64, b::Float64, n::Int; alpha::Float64=DEFAULT_POWER_ALPHA)
    len = b - a
    if alpha == DEFAULT_POWER_ALPHA
        if n == 32
            xs = a .+ len .* _POWER_RIGHT_T_32
            wx = _STD_32_WEIGHTS .* len .* _POWER_RIGHT_T_PRIME_32
            return xs, wx
        elseif n == 16
            xs = a .+ len .* _POWER_RIGHT_T_16
            wx = _STD_16_WEIGHTS .* len .* _POWER_RIGHT_T_PRIME_16
            return xs, wx
        end
    end

    us, ws = GaussLegendre.gauleg(-1.0, 1.0, n)
    inv_alpha = 1.0 / alpha
    u_half = (us .+ 1) ./ 2
    t = 1 .- (1 .- u_half) .^ inv_alpha
    t_prime = inv_alpha .* (1 .- u_half) .^ (inv_alpha - 1) ./ 2
    xs = a .+ len .* t
    wx = ws .* len .* t_prime
    return xs, wx
end

"""生成 DE (tanh-sinh) 变换节点（默认参数走预计算）"""
function de_nodes(a::Float64, b::Float64, n::Int; h::Float64=DEFAULT_DE_H)
    if n == 32 && h == DEFAULT_DE_H
        half = (b - a) / 2
        center = (a + b) / 2
        xs = center .+ half .* _DE_PHI
        wx = half .* _DE_PHI_PRIME
        return xs, wx
    end

    ks = collect(-n÷2:n÷2)
    ts = ks .* h
    half = (b - a) / 2
    center = (a + b) / 2
    xs = center .+ half .* tanh.(π/2 .* sinh.(ts))
    wx = h .* half .* (π/2) .* cosh.(ts) .* (sech.(π/2 .* sinh.(ts))).^2
    return xs, wx
end

"""奇点位置枚举"""
@enum SingularityPosition SING_NONE SING_LEFT SING_RIGHT SING_BOTH

"""根据奇点位置选择变换"""
function hybrid_nodes(a::Float64, b::Float64, n::Int, sing_pos::SingularityPosition;
    alpha::Float64=DEFAULT_POWER_ALPHA, beta::Float64=DEFAULT_CLUSTER_BETA, h::Float64=DEFAULT_DE_H)
    if sing_pos == SING_LEFT
        return power_left_nodes(a, b, n; alpha=alpha)
    elseif sing_pos == SING_RIGHT
        return power_right_nodes(a, b, n; alpha=alpha)
    elseif sing_pos == SING_BOTH
        return de_nodes(a, b, n; h=h)
    else
        return power_left_nodes(a, b, n; alpha=alpha)
    end
end

# ============================================================================
# 高性能：无分配的 hybrid 区间求积
#
# 设计：直接复用本文件里预计算的映射数组（_POWER_* / _DE_*），
# 在循环内计算节点与权重，避免为每个区间分配 xs/wx 向量。
#
# 说明：
# - SING_LEFT  -> power_left
# - SING_RIGHT -> power_right
# - SING_BOTH  -> DE (tanh-sinh)
# - SING_NONE  -> power_left（经验上对端点更稳）
#
# 若传入非默认参数或非常用 n，会自动回退到生成节点的实现。
# ============================================================================

@inline function _integrate_power_left_32(f::F, a::Float64, b::Float64) where {F}
    len = b - a
    total = 0.0
    @inbounds @simd for i in eachindex(_STD_32_WEIGHTS)
        x = a + len * _POWER_LEFT_T_32[i]
        w = _STD_32_WEIGHTS[i] * len * _POWER_LEFT_T_PRIME_32[i]
        v = f(x)
        if isfinite(v)
            total += w * v
        end
    end
    return total
end

@inline function _integrate_power_left_16(f::F, a::Float64, b::Float64) where {F}
    len = b - a
    total = 0.0
    @inbounds @simd for i in eachindex(_STD_16_WEIGHTS)
        x = a + len * _POWER_LEFT_T_16[i]
        w = _STD_16_WEIGHTS[i] * len * _POWER_LEFT_T_PRIME_16[i]
        v = f(x)
        if isfinite(v)
            total += w * v
        end
    end
    return total
end

@inline function _integrate_power_right_32(f::F, a::Float64, b::Float64) where {F}
    len = b - a
    total = 0.0
    @inbounds @simd for i in eachindex(_STD_32_WEIGHTS)
        x = a + len * _POWER_RIGHT_T_32[i]
        w = _STD_32_WEIGHTS[i] * len * _POWER_RIGHT_T_PRIME_32[i]
        v = f(x)
        if isfinite(v)
            total += w * v
        end
    end
    return total
end

@inline function _integrate_power_right_16(f::F, a::Float64, b::Float64) where {F}
    len = b - a
    total = 0.0
    @inbounds @simd for i in eachindex(_STD_16_WEIGHTS)
        x = a + len * _POWER_RIGHT_T_16[i]
        w = _STD_16_WEIGHTS[i] * len * _POWER_RIGHT_T_PRIME_16[i]
        v = f(x)
        if isfinite(v)
            total += w * v
        end
    end
    return total
end

@inline function _integrate_de_default(f::F, a::Float64, b::Float64) where {F}
    half = (b - a) / 2
    center = (a + b) / 2
    total = 0.0
    @inbounds @simd for i in eachindex(_DE_PHI)
        x = center + half * _DE_PHI[i]
        w = half * _DE_PHI_PRIME[i]
        v = f(x)
        if isfinite(v)
            total += w * v
        end
    end
    return total
end

"""在区间 [a,b] 上用 hybrid 规则积分（高性能路径，尽量无分配）。"""
@inline function integrate_hybrid_interval(f::F, a::Float64, b::Float64, sing_pos::SingularityPosition;
    n::Int=32, alpha::Float64=DEFAULT_POWER_ALPHA, h::Float64=DEFAULT_DE_H) where {F}

    if !(isfinite(a) && isfinite(b)) || b <= a
        return 0.0
    end

    # 常用默认参数走无分配快路径
    if alpha == DEFAULT_POWER_ALPHA && h == DEFAULT_DE_H
        if sing_pos == SING_BOTH
            return _integrate_de_default(f, a, b)
        elseif sing_pos == SING_RIGHT
            return (n == 16) ? _integrate_power_right_16(f, a, b) : _integrate_power_right_32(f, a, b)
        else
            # SING_LEFT / SING_NONE
            return (n == 16) ? _integrate_power_left_16(f, a, b) : _integrate_power_left_32(f, a, b)
        end
    end

    # 回退：生成节点再积分（仍然正确，但会分配）
    xs, wx = hybrid_nodes(a, b, n, sing_pos; alpha=alpha, h=h)
    total = 0.0
    @inbounds for i in eachindex(xs)
        v = f(xs[i])
        if isfinite(v)
            total += wx[i] * v
        end
    end
    return total
end
