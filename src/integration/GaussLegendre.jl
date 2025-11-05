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

include("../../src/Constants_PNJL.jl")
using .Constants_PNJL: Λ_inv_fm
using FastGaussQuadrature

export gauleg

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
    
    # 获取标准区间 [-1, 1] 上的节点和权重
    nodes_std, weights_std = gausslegendre(n)
    
    # 变换到区间 [a, b]
    # x = ((b-a)*t + (b+a))/2, t ∈ [-1, 1]
    nodes = @. ((b - a) * nodes_std + (b + a)) / 2.0
    
    # 权重需要乘以雅可比行列式 (b-a)/2
    weights = @. weights_std * (b - a) / 2.0
    
    return nodes, weights
end

const DEFAULT_momentum_POINTS = 64
const DEFAULT_theta_POINTS = 32
function build_default_nodes_weights()
    cosθ_NODES, cosθ_WEIGHTS = gauleg(-1.0, 1.0, DEFAULT_theta_POINTS)
    momentum_NODES_toINF, momentum_WEIGHTS_toINF = gauleg(0.0, 20.0, DEFAULT_momentum_POINTS)
    momentum_NODES_toPNJLΛ, momentum_WEIGHTS_toPNJLΛ = gauleg(0.0, Λ_inv_fm, DEFAULT_momentum_POINTS)
    return cosθ_NODES, cosθ_WEIGHTS,
        momentum_NODES_toINF, momentum_WEIGHTS_toINF,
        momentum_NODES_toPNJLΛ, momentum_WEIGHTS_toPNJLΛ
end
export build_default_nodes_weights
end # module GaussLegendre
