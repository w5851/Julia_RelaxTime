# GaussLegendre 模块 API 文档

## 模块概述

`GaussLegendre` 模块提供高斯-勒让德积分节点和权重的生成功能，用于高精度数值积分计算。

## 依赖

- `FastGaussQuadrature.jl` - 用于计算标准高斯-勒让德节点

## 单位约定

遵循项目统一的自然单位制 (ℏ = c = 1)：

- **积分区间** (a, b): fm⁻¹ 单位
- **节点位置** xᵢ: fm⁻¹ 单位  
- **积分权重** wᵢ: fm⁻¹ 单位

---

## API 参考

### `gauleg(a, b, n)`

生成区间 (a, b) 上的 n 个高斯-勒让德积分节点和对应权重。

#### 函数签名

```julia
gauleg(a::Float64, b::Float64, n::Int) -> (nodes::Vector{Float64}, weights::Vector{Float64})
```

#### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `a` | `Float64` | 积分区间左端点 | fm⁻¹ |
| `b` | `Float64` | 积分区间右端点 | fm⁻¹ |
| `n` | `Int` | 节点数量（必须 > 0） | 无量纲 |

#### 返回值

返回一个元组 `(nodes, weights)`：

- `nodes::Vector{Float64}`: 长度为 n 的向量，包含 n 个积分节点，位于区间 (a, b) 内
- `weights::Vector{Float64}`: 长度为 n 的向量，包含对应的积分权重

#### 数学原理

高斯-勒让德积分公式在区间 [a, b] 上的形式：

$$
\int_a^b f(x) \, dx \approx \sum_{i=1}^{n} w_i f(x_i)
$$

其中：
- $x_i$ 是勒让德多项式 $P_n(x)$ 在标准区间 [-1, 1] 上的零点，经过线性变换映射到 [a, b]
- $w_i$ 是对应的积分权重

**坐标变换：**

从标准区间 [-1, 1] 到任意区间 [a, b]：

$$
x = \frac{(b-a)t + (b+a)}{2}, \quad t \in [-1, 1]
$$

$$
w' = w \times \frac{b-a}{2}
$$

#### 精度特性

- 对于 n 个节点的高斯-勒让德积分，可以精确积分至多 2n-1 次的多项式
- 对于光滑函数，积分误差随节点数呈**指数衰减**
- 典型应用中，n = 10-100 个节点已足够精确

#### 使用示例

**示例 1：基本用法**

```julia
using RelaxTime
using RelaxTime.GaussLegendre

# 生成区间 [0.0, 1.0] 上的 5 个节点
nodes, weights = gauleg(0.0, 1.0, 5)

println("节点: ", nodes)
println("权重: ", weights)
```

**示例 2：计算简单积分**

```julia
# 计算 ∫₀¹ x² dx = 1/3
f(x) = x^2
nodes, weights = gauleg(0.0, 1.0, 5)
result = sum(weights .* f.(nodes))
println("数值结果: ", result)
println("解析解: ", 1/3)
println("相对误差: ", abs(result - 1/3) / (1/3))
```

**示例 3：物理应用（自然单位制）**

```julia
# 计算动量积分 ∫₀^Λf f(k) dk，Λf ≈ 3.05 fm⁻¹
Λf = 3.05  # fm⁻¹
n = 50     # 使用 50 个节点

# 生成节点
k_nodes, k_weights = gauleg(0.0, Λf, n)

# 定义被积函数（例如：费米-狄拉克分布）
T = 0.1    # 温度，fm⁻¹ 单位
μ = 0.3    # 化学势，fm⁻¹ 单位
fermi_dirac(k) = 1.0 / (exp((k - μ) / T) + 1.0)

# 计算积分
integral = sum(k_weights .* fermi_dirac.(k_nodes))
println("积分结果: ", integral, " fm⁻¹")
```

**示例 4：验证节点数对精度的影响**

```julia
# 计算 ∫₀^π sin(x) dx = 2
f(x) = sin(x)
exact = 2.0

for n in [5, 10, 20, 50]
    nodes, weights = gauleg(0.0, π, n)
    result = sum(weights .* f.(nodes))
    error = abs(result - exact)
    println("n=$n: 结果=$result, 误差=$error")
end
```

#### 异常处理

函数在以下情况会抛出 `ArgumentError`：

1. **节点数无效**: `n ≤ 0`
   ```julia
   gauleg(0.0, 1.0, 0)  # ArgumentError: 节点数 n 必须大于 0
   ```

2. **区间无效**: `a ≥ b`
   ```julia
   gauleg(1.0, 0.0, 10)  # ArgumentError: 区间端点必须满足 a < b
   ```

#### 性能建议

1. **节点数选择**：
   - 光滑函数：n = 10-20 通常足够
   - 振荡函数：增加到 n = 50-100
   - 奇异点附近：考虑分段积分

2. **重复使用节点**：
   - 如果多次在同一区间积分，预先计算并存储节点和权重

3. **向量化计算**：
   - 利用 Julia 的广播语法 `.` 进行向量化操作
   - 示例：`result = sum(weights .* f.(nodes))`

#### 数值稳定性

- 使用 `FastGaussQuadrature.jl` 提供的高精度算法
- 节点计算基于 Golub-Welsch 算法，数值稳定
- 对于极大的 n (>1000)，可能出现数值精度问题

#### 参考文献

1. **Golub, G. H., & Welsch, J. H. (1969)**  
   "Calculation of Gauss Quadrature Rules"  
   *Mathematics of Computation*, 23(106), 221-230

2. **FastGaussQuadrature.jl**  
   https://github.com/JuliaApproximation/FastGaussQuadrature.jl

3. **Press, W. H., et al. (2007)**  
   *Numerical Recipes: The Art of Scientific Computing* (3rd ed.)  
   Cambridge University Press

---

## 预定义常量

模块提供了预定义的积分节点和权重，可以直接使用：

### 角度积分节点

#### `DEFAULT_COSΘ_NODES`, `DEFAULT_COSΘ_WEIGHTS`

**全区间角度积分** (cosθ ∈ [-1, 1])

```julia
using .GaussLegendre: DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS

# 用于一般的角度积分
result = sum(DEFAULT_COSΘ_WEIGHTS .* f.(DEFAULT_COSΘ_NODES))
```

- **节点数**: 32
- **积分区间**: [-1, 1]
- **用途**: 通用的角度积分

#### `DEFAULT_COSΘ_HALF_NODES`, `DEFAULT_COSΘ_HALF_WEIGHTS`

**半区间角度积分** (cosθ ∈ [0, 1])

```julia
using .GaussLegendre: DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS

# 对于对称函数 f(x) = f(-x)，利用对称性:
# ∫₋₁¹ f(cosθ) d(cosθ) = 2 ∫₀¹ f(cosθ) d(cosθ)
result = 2.0 * sum(DEFAULT_COSΘ_HALF_WEIGHTS .* f.(DEFAULT_COSΘ_HALF_NODES))
```

- **节点数**: 32
- **积分区间**: [0, 1]
- **用途**: 对称被积函数的高精度积分
- **优势**: 相同节点数下，精度比全区间方法更高（节点在半区间内更密集）

**适用场景**:
- ✓ 偶函数: f(x) = f(-x)
- ✓ 各向同性分布
- ✓ 偶数阶勒让德多项式
- ✗ 奇函数或非对称函数

**精度对比示例**:

```julia
# 测试函数: f(x) = x²
exact = 2.0/3.0

# 方法1: 全区间 (32节点)
result_full = sum(DEFAULT_COSΘ_WEIGHTS .* (x -> x^2).(DEFAULT_COSΘ_NODES))
error_full = abs(result_full - exact) / exact * 100  # ~3e-14 %

# 方法2: 半区间对称 (32节点)
result_half = 2.0 * sum(DEFAULT_COSΘ_HALF_WEIGHTS .* (x -> x^2).(DEFAULT_COSΘ_HALF_NODES))
error_half = abs(result_half - exact) / exact * 100  # ~0 % (机器精度)

# 精度提升: 可达数倍到数十倍
```

### 动量积分节点

#### `DEFAULT_MOMENTUM_NODES`, `DEFAULT_MOMENTUM_WEIGHTS`

**热力学动量积分** (p ∈ [0, 10] fm⁻¹)

```julia
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS

# 计算热力学量的动量积分
result = sum(DEFAULT_MOMENTUM_WEIGHTS .* f.(DEFAULT_MOMENTUM_NODES))
```

- **节点数**: 64
- **积分区间**: [0, 10] fm⁻¹
- **用途**: 热力学积分（费米-狄拉克分布在 p ~ 10 fm⁻¹ 时已充分衰减）

#### PNJL Λ-截断节点

Λ 取决于 PNJL 配置文件，目前默认的 `[0, Λ]` 节点仅在 `PNJL.AnisoGapSolver` 内部使用。若需要同样的节点分布，可通过下述方式自行构造：

```julia
include("src/integration/GaussLegendre.jl")
include("src/Constants_PNJL.jl")
using .GaussLegendre: gauleg
using .Constants_PNJL: Λ_inv_fm

nodes_Λ, weights_Λ = gauleg(0.0, Λ_inv_fm, 64)
```

这样可以完全自定义节点数，同时显式绑定当前使用的 PNJL 参数集。

---

## 模块导入

```julia
# 方式 1：通过主模块导入
using RelaxTime
using RelaxTime.GaussLegendre

# 方式 2：直接包含文件
include("src/integration/GaussLegendre.jl")
using .GaussLegendre

# 方式 3：只导入需要的常量
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS,
                       DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS
```

---

## 常见问题

**Q: 如何选择合适的节点数 n？**

A: 从较小的 n 开始（如 10），逐步增加并检查结果收敛性。如果相邻两次结果变化小于所需精度，则当前 n 已足够。

**Q: 为什么权重单位是 fm⁻¹？**

A: 权重的单位等于积分区间的单位，使得 Σwᵢf(xᵢ) 与 ∫f(x)dx 单位一致。

**Q: 可以用于无穷区间积分吗？**

A: 高斯-勒让德积分适用于有限区间。对于无穷区间，建议：
- 使用高斯-拉盖尔积分（半无穷区间）
- 使用高斯-埃尔米特积分（全无穷区间）
- 或进行变量替换将无穷区间映射到有限区间

**Q: 如何处理被积函数的奇点？**

A: 
- 如果奇点在积分区间内，考虑分段积分避开奇点
- 如果奇点在端点，考虑使用特殊的积分方法（如高斯-雅可比积分）
- 对于可积奇点，可以解析处理奇异部分

---

## 更新日志

### v0.1.0 (2025-10-14)
- 初始版本
- 实现 `gauleg` 函数
- 支持任意有限区间 [a, b]
- 完整的单位约定和文档
