"""
# GaussLegendre.jl

高斯-勒让德积分节点生成模块

## 单位约定
- 积分区间 (a, b): fm⁻¹ 单位（与被积函数的自变量单位一致）
- 节点位置 xᵢ: fm⁻¹ 单位
- 权重 wᵢ: fm⁻¹ 单位（使得 ∫f(x)dx ≈ Σwᵢf(xᵢ) 的单位匹配）

## 依赖
本模块依赖 FastGaussQuadrature.jl 包来计算高斯-勒让德节点和权重。
"""
module GaussLegendre

using FastGaussQuadrature
using Base.Threads: ReentrantLock

export gauleg, gausslegendre
export standard_nodes_weights
export affine_transform, transform_standard8

const _STANDARD_NODES_WEIGHTS_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const _STANDARD_NODES_WEIGHTS_LOCK = ReentrantLock()

"""
    standard_nodes_weights(n::Int) -> (nodes_std, weights_std)

返回标准区间 [-1, 1] 上的 n 个高斯-勒让德节点与权重，并在模块内缓存。

注意：返回的向量应视为只读。
"""
function standard_nodes_weights(n::Int)
    lock(_STANDARD_NODES_WEIGHTS_LOCK)
    try
        return get!(_STANDARD_NODES_WEIGHTS_CACHE, n) do
            FastGaussQuadrature.gausslegendre(n)
        end
    finally
        unlock(_STANDARD_NODES_WEIGHTS_LOCK)
    end
end

"""
    gauleg(a::Float64, b::Float64, n::Int) -> (nodes, weights)

生成区间 (a, b) 上的 n 个高斯-勒让德积分节点和对应权重。

# 参数
- `a::Float64`: 积分区间左端点（fm⁻¹ 单位）
- `b::Float64`: 积分区间右端点（fm⁻¹ 单位）
- `n::Int`: 节点数量（必须 > 0）

# 返回值
- `nodes::Vector{Float64}`: n 个积分节点，位于区间 (a, b) 内（fm⁻¹ 单位）
- `weights::Vector{Float64}`: 对应的积分权重（fm⁻¹ 单位）

# 数学原理
高斯-勒让德积分公式：
```math
∫ₐᵇ f(x)dx ≈ Σᵢ₌₁ⁿ wᵢ f(xᵢ)
```

标准区间 [-1, 1] 到任意区间 [a, b] 的变换：
```math
x = ((b-a)t + (b+a))/2,  t ∈ [-1, 1]
w' = w × (b-a)/2
```

# 示例
```julia
# 计算区间 [0.0, 1.0] 上的 5 个节点
nodes, weights = gauleg(0.0, 1.0, 5)

# 使用节点进行数值积分
f(x) = x^2
result = sum(weights .* f.(nodes))
# 应该接近 ∫₀¹ x² dx = 1/3 ≈ 0.333...
```

# 注意事项
- 节点数 n 必须大于 0
- 要求 a < b
- 对于光滑函数，节点数越多精度越高（指数收敛）
- 典型应用中 n = 10-100 已足够精确

# 异常
- `ArgumentError`: 如果 n ≤ 0 或 a ≥ b

# 参考
- FastGaussQuadrature.jl: https://github.com/JuliaApproximation/FastGaussQuadrature.jl
- Golub & Welsch (1969): "Calculation of Gauss Quadrature Rules"
"""
function gauleg(a::Float64, b::Float64, n::Int)
    # 参数检查
    if n <= 0
        throw(ArgumentError("节点数 n 必须大于 0，当前值: $n"))
    end

    if a >= b
        throw(ArgumentError("区间端点必须满足 a < b，当前值: a=$a, b=$b"))
    end

    # 获取标准区间 [-1,1] 上的节点和权重（带缓存）
    nodes_std, weights_std = standard_nodes_weights(n)

    # 变换到区间 [a, b]
    # x = ((b-a)*t + (b+a))/2, t ∈ [-1, 1]
    scale = (b - a) / 2.0
    shift = (b + a) / 2.0
    nodes = similar(nodes_std)
    weights = similar(weights_std)
    @inbounds @simd for i in 1:n
        nodes[i] = scale * nodes_std[i] + shift
        weights[i] = scale * weights_std[i]
    end

    return nodes, weights
end

const DEFAULT_momentum_POINTS = 64
const DEFAULT_theta_POINTS = 8

# 预计算一个用于标准区间 [-1,1] 的固定 8 节点（及权重），方便在需要固定节点数的场景中直接复用。
# 这样上层代码可以 `include("GaussLegendre.jl")` 后直接使用 `DEFAULT_STANDARD_8_NODES/WEIGHTS`，
# 对区间 [a,b] 只需做仿射变换，避免重复求解标准节点。
const DEFAULT_STANDARD_8_NODES, DEFAULT_STANDARD_8_WEIGHTS = gausslegendre(8)
# 预计算 16 节点与 32 节点版本，供需要更高精度的固定节点情形使用
const DEFAULT_STANDARD_16_NODES, DEFAULT_STANDARD_16_WEIGHTS = gausslegendre(16)
const DEFAULT_STANDARD_32_NODES, DEFAULT_STANDARD_32_WEIGHTS = gausslegendre(32)

"""对标准区间节点/权重做仿射变换到区间 [a, b]

返回 (nodes, weights)
"""
function affine_transform(nodes_std::Vector{Float64}, weights_std::Vector{Float64}, a::Float64, b::Float64)
    nodes = @. ((b - a) * nodes_std + (b + a)) / 2.0
    weights = @. weights_std * (b - a) / 2.0
    return nodes, weights
end

"""把预计算的标准 8 节点映射到区间 [a,b] 并返回 (nodes, weights)"""
function transform_standard8(a::Float64, b::Float64)
    return affine_transform(DEFAULT_STANDARD_8_NODES, DEFAULT_STANDARD_8_WEIGHTS, a, b)
end

"""把预计算的标准 16 节点映射到区间 [a,b] 并返回 (nodes, weights)"""
function transform_standard16(a::Float64, b::Float64)
    return affine_transform(DEFAULT_STANDARD_16_NODES, DEFAULT_STANDARD_16_WEIGHTS, a, b)
end

"""把预计算的标准 32 节点映射到区间 [a,b] 并返回 (nodes, weights)"""
function transform_standard32(a::Float64, b::Float64)
    return affine_transform(DEFAULT_STANDARD_32_NODES, DEFAULT_STANDARD_32_WEIGHTS, a, b)
end

"""默认的角度积分节点和权重 (cosθ ∈ [-1, 1]，用于一般角度积分)"""
const DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS = gauleg(-1.0, 1.0, DEFAULT_theta_POINTS)

"""默认的半区间角度积分节点和权重 (cosθ ∈ [0, 1]，用于对称被积函数的高精度积分)

当被积函数关于 cosθ = 0 对称时（即 f(cosθ) = f(-cosθ)），可以利用对称性：
∫₋₁¹ f(cosθ) d(cosθ) = 2 ∫₀¹ f(cosθ) d(cosθ)

使用相同节点数时，半区间积分的精度更高，因为节点更密集地分布在 [0, 1] 区间。
"""
const DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS = gauleg(0.0, 1.0, DEFAULT_theta_POINTS)

"""默认的动量积分节点和权重 (p ∈ [0, 10] fm⁻¹，用于热力学积分)"""
const DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS = gauleg(0.0, 10.0, DEFAULT_momentum_POINTS)

export DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS
export DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS
export DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
export DEFAULT_STANDARD_8_NODES, DEFAULT_STANDARD_8_WEIGHTS
export DEFAULT_STANDARD_16_NODES, DEFAULT_STANDARD_16_WEIGHTS
export DEFAULT_STANDARD_32_NODES, DEFAULT_STANDARD_32_WEIGHTS
export transform_standard16, transform_standard32

end # module GaussLegendre
