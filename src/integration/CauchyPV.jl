"""
# CauchyPV.jl

柯西主值积分节点生成模块

## 单位约定
- 积分区间 (a, b): fm⁻¹ 单位
- 奇点位置 x₀: fm⁻¹ 单位
- 节点位置 xᵢ: fm⁻¹ 单位
- 权重 wᵢ: fm⁻¹ 单位

## 原理
柯西主值积分用于计算含有奇点的积分：
```math
P.V. ∫ₐᵇ f(x)/(x-x₀) dx = lim_{ε→0} [∫ₐ^(x₀-ε) f(x)/(x-x₀)dx + ∫_(x₀+ε)^b f(x)/(x-x₀)dx]
```

本模块使用复化梯形法则，在奇点附近避开一个小区间（0.001h），
将积分区间分为两部分：[a, x₀-δ] 和 [x₀+δ, b]。
"""
module CauchyPV

export cpvi

"""
    cpvi(a::Float64, b::Float64, x0::Float64, n::Int) -> (nodes, weights, n_actual)

生成柯西主值积分的节点和权重，避开奇点 x₀。

使用复化梯形法则，在区间 [a, b] 上生成节点，在奇点 x₀ 附近保持一个小间隙（0.001h）。

# 参数
- `a::Float64`: 积分区间左端点（fm⁻¹ 单位）
- `b::Float64`: 积分区间右端点（fm⁻¹ 单位）
- `x0::Float64`: 奇点位置（fm⁻¹ 单位），必须满足 a < x₀ < b
- `n::Int`: 期望的节点数量（实际节点数可能略有不同）

# 返回值
- `nodes::Vector{Float64}`: 积分节点（fm⁻¹ 单位）
- `weights::Vector{Float64}`: 对应的积分权重（fm⁻¹ 单位）
- `n_actual::Int`: 实际生成的节点数

# 算法原理

## 步长计算
```math
h = (b - a) / (n - 1)
```

## 区间分割
- 左区间：[a, x₀ - 0.001h]，包含 i 个节点
- 右区间：[x₀ + 0.001h, b]，包含 j 个节点
- 总节点数：n_actual = i + j

## 权重分配
使用复化梯形法则：
- 内部节点权重：w = h
- 端点权重：w = h/2
- 奇点附近的端点（索引 i 和 i+1）权重：w = h/2

# 示例

**示例 1：基本用法**
```julia
using RelaxTime.CauchyPV

# 在区间 [0.0, 2.0] 上生成节点，避开奇点 x₀ = 1.0
nodes, weights, n_actual = cpvi(0.0, 2.0, 1.0, 100)
println("实际节点数: ", n_actual)
```

**示例 2：计算柯西主值积分**
```julia
# 计算 P.V. ∫₀² 1/(x-1) dx = 0 （奇点在 x=1）
a, b, x0 = 0.0, 2.0, 1.0
nodes, weights, n_actual = cpvi(a, b, x0, 200)

# 被积函数（不包含奇异因子 1/(x-x₀)）
f(x) = 1.0

# 数值积分（包含奇异因子）
result = sum(weights[i] * f(nodes[i]) / (nodes[i] - x0) for i in 1:n_actual)
println("积分结果: ", result)  # 应该接近 0
```

**示例 3：物理应用 - 色散关系**
```julia
# 计算色散积分（自然单位制）
T = 0.1      # 温度，fm⁻¹
μ = 0.3      # 化学势，fm⁻¹
k0 = 0.5     # 外动量，fm⁻¹
Λf = 3.05    # 截断参数，fm⁻¹

# 在 [0, Λf] 上积分，k0 是奇点
nodes, weights, n_actual = cpvi(0.0, Λf, k0, 100)

# 被积函数：谱函数
function spectral(k)
    ω = sqrt(k^2 + 0.3^2)  # 粒子能量
    return 1.0 / (exp((ω - μ) / T) + 1.0)
end

# 色散积分
dispersion = sum(weights[i] * spectral(nodes[i]) / (nodes[i] - k0) 
                 for i in 1:n_actual)
println("色散积分: ", dispersion, " fm")
```

# 注意事项

1. **奇点位置**：必须保证 a < x₀ < b，否则抛出异常
2. **节点数调整**：实际节点数 n_actual 可能与输入 n 略有不同
3. **步长选择**：步长 h = (b-a)/(n-1)，较小的步长提高精度
4. **奇点间隙**：在 x₀ 附近保持 0.001h 的间隙，避免数值不稳定
5. **对称性**：对于对称的被积函数和居中的奇点，积分结果应接近零

# 误差来源

- 梯形法则的截断误差：O(h²)
- 奇点附近的数值误差
- 有限间隙（0.001h）引入的误差

# 异常
- `ArgumentError`: 如果 n ≤ 1、a ≥ b 或 x₀ 不在 (a, b) 内

# 参考
- Numerical Recipes in Fortran 90
- Cauchy Principal Value Integration Methods
"""
function cpvi(a::Float64, b::Float64, x0::Float64, n::Int)
    # 参数检查
    if n <= 1
        throw(ArgumentError("节点数 n 必须大于 1，当前值: $n"))
    end
    
    if a >= b
        throw(ArgumentError("区间端点必须满足 a < b，当前值: a=$a, b=$b"))
    end
    
    if x0 <= a || x0 >= b
        throw(ArgumentError("奇点 x₀ 必须在区间 (a, b) 内，当前值: x₀=$x0, a=$a, b=$b"))
    end
    
    # 计算步长
    h = (b - a) / ((n - 1) * 1.0)
    
    # 计算左区间 [a, x₀ - 0.001h] 的节点数
    i = floor(Int, ((x0 - 0.001 * h) - a) / h) + 1
    
    # 计算右区间 [x₀ + 0.001h, b] 的节点数
    j = floor(Int, (b - (x0 + 0.001 * h)) / h) + 1
    
    # 实际节点总数
    n_actual = i + j
    
    # 初始化节点和权重数组
    nodes = zeros(Float64, n_actual)
    weights = fill(h, n_actual)
    
    # 设置权重（梯形法则）
    weights[1] = h / 2.0              # 左端点
    if i > 0
        weights[i] = h / 2.0          # 左区间右端点
    end
    if j > 0 && i + 1 <= n_actual
        weights[i + 1] = h / 2.0      # 右区间左端点
    end
    weights[n_actual] = h / 2.0       # 右端点
    
    # 生成节点位置
    # 左区间：从 x₀ - 0.001h 向左，步长为 -h
    if i > 0
        for k in 1:i
            nodes[i - k + 1] = (x0 - 0.001 * h) - (k - 1) * h
        end
    end
    
    # 右区间：从 x₀ + 0.001h 向右，步长为 h
    if j > 0
        for k in 1:j
            nodes[i + k] = (x0 + 0.001 * h) + (k - 1) * h
        end
    end
    
    return nodes, weights, n_actual
end

end # module CauchyPV
