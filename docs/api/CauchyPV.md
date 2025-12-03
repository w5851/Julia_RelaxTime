# CauchyPV 模块 API 文档

## 模块概述

`CauchyPV` 模块提供柯西主值积分（Cauchy Principal Value Integration）的节点和权重生成功能。用于计算包含奇点的积分，通过在奇点附近避开一个小区间，使用复化梯形法则进行数值积分。

## 依赖

无外部依赖，仅使用 Julia 标准库。

## 单位约定

遵循项目统一的自然单位制 (ℏ = c = 1)：

- **积分区间** (a, b): fm⁻¹ 单位
- **奇点位置** x₀: fm⁻¹ 单位
- **节点位置** xᵢ: fm⁻¹ 单位  
- **积分权重** wᵢ: fm⁻¹ 单位

---

## 数学背景

### 柯西主值积分

对于含有简单极点的积分，柯西主值定义为：

$$
P.V. \int_a^b \frac{f(x)}{x - x_0} dx = \lim_{\varepsilon \to 0} \left[ \int_a^{x_0-\varepsilon} \frac{f(x)}{x - x_0} dx + \int_{x_0+\varepsilon}^b \frac{f(x)}{x - x_0} dx \right]
$$

其中 $x_0 \in (a, b)$ 是奇点。

### 数值实现

本模块使用复化梯形法则，将积分区间分为两部分：
- 左区间：$[a, x_0 - \delta]$
- 右区间：$[x_0 + \delta, b]$

其中 $\delta = 0.001h$，$h$ 是基本步长。

---

## API 参考

### `cpvi(a, b, x0, n)`

生成柯西主值积分的节点和权重。

#### 函数签名

```julia
cpvi(a::Float64, b::Float64, x0::Float64, n::Int) -> (nodes::Vector{Float64}, weights::Vector{Float64}, n_actual::Int)
```

#### 参数

| 参数 | 类型 | 说明 | 单位 | 约束 |
|------|------|------|------|------|
| `a` | `Float64` | 积分区间左端点 | fm⁻¹ | a < b |
| `b` | `Float64` | 积分区间右端点 | fm⁻¹ | a < b |
| `x0` | `Float64` | 奇点位置 | fm⁻¹ | a < x₀ < b |
| `n` | `Int` | 期望节点数 | 无量纲 | n > 1 |

#### 返回值

返回一个元组 `(nodes, weights, n_actual)`：

- `nodes::Vector{Float64}`: 长度为 n_actual 的向量，包含积分节点
- `weights::Vector{Float64}`: 长度为 n_actual 的向量，包含对应的积分权重
- `n_actual::Int`: 实际生成的节点数（可能与输入 n 略有不同）

#### 算法详解

**1. 步长计算**

$$
h = \frac{b - a}{n - 1}
$$

**2. 区间划分**

- 左区间节点数：$i = \lfloor \frac{(x_0 - 0.001h) - a}{h} \rfloor + 1$
- 右区间节点数：$j = \lfloor \frac{b - (x_0 + 0.001h)}{h} \rfloor + 1$
- 实际总节点数：$n_{actual} = i + j$

**3. 节点生成**

左区间节点（从奇点向左）：
$$
x_k = (x_0 - 0.001h) - (k-1)h, \quad k = 1, 2, \ldots, i
$$

右区间节点（从奇点向右）：
$$
x_{i+k} = (x_0 + 0.001h) + (k-1)h, \quad k = 1, 2, \ldots, j
$$

**4. 权重分配（复化梯形法则）**

- 内部节点：$w_k = h$
- 左端点：$w_1 = h/2$
- 左区间右端点：$w_i = h/2$
- 右区间左端点：$w_{i+1} = h/2$
- 右端点：$w_{n_{actual}} = h/2$

#### 使用示例

**示例 1：基本用法**

```julia
using RelaxTime
using RelaxTime.CauchyPV

# 在区间 [0.0, 2.0] 上生成节点，奇点在 x₀ = 1.0
nodes, weights, n_actual = cpvi(0.0, 2.0, 1.0, 100)

println("期望节点数: 100")
println("实际节点数: ", n_actual)
println("节点范围: [", minimum(nodes), ", ", maximum(nodes), "]")
println("权重和: ", sum(weights))
```

**示例 2：验证对称性**

```julia
# 对于对称函数和居中奇点，主值积分应为零
a, b, x0 = -1.0, 1.0, 0.0
nodes, weights, n_actual = cpvi(a, b, x0, 200)

# f(x) = x（奇函数）
f(x) = x
result = sum(weights[i] * f(nodes[i]) / nodes[i] for i in 1:n_actual)
println("奇函数积分（应接近0）: ", result)

# f(x) = x²（偶函数）
g(x) = x^2
result2 = sum(weights[i] * g(nodes[i]) / nodes[i] for i in 1:n_actual)
println("偶函数积分: ", result2)
```

**示例 3：计算 Hilbert 变换**

```julia
# Hilbert 变换: H[f](x₀) = (1/π) P.V. ∫_{-∞}^{∞} f(x)/(x₀-x) dx
# 对于有限区间近似

function hilbert_transform(f, x0, a, b, n)
    nodes, weights, n_actual = cpvi(a, b, x0, n)
    integral = sum(weights[i] * f(nodes[i]) / (x0 - nodes[i]) 
                   for i in 1:n_actual)
    return integral / π
end

# 测试函数：f(x) = exp(-x²)
f(x) = exp(-x^2)
x0 = 0.5
result = hilbert_transform(f, x0, -5.0, 5.0, 500)
println("Hilbert 变换在 x₀=$x0: ", result)
```

**示例 4：物理应用 - 色散关系（Kramers-Kronig）**

```julia
# 在量子场论中，实部和虚部通过色散关系联系
# Re[Π(k₀)] = (1/π) P.V. ∫ Im[Π(k)] / (k₀ - k) dk

T = 0.15     # 温度，fm⁻¹
μ = 0.3      # 化学势，fm⁻¹
k0 = 0.5     # 外动量，fm⁻¹
Λf = 3.05    # 截断，fm⁻¹

# 偏振函数的虚部（简化模型）
function im_polarization(k)
    ω = sqrt(k^2 + 0.14^2)  # 准粒子能量（夸克质量 ≈ 0.14 fm⁻¹）
    nF = 1.0 / (exp((ω - μ) / T) + 1.0)  # 费米分布
    return k^2 * nF * (1 - nF)
end

# 通过色散关系计算实部
nodes, weights, n_actual = cpvi(0.0, Λf, k0, 200)
re_polarization = (1/π) * sum(weights[i] * im_polarization(nodes[i]) / (k0 - nodes[i]) 
                               for i in 1:n_actual)

println("Re[Π(k₀=$k0)] = ", re_polarization, " fm²")
```

**示例 5：收敛性测试**

```julia
# 测试节点数对精度的影响
a, b, x0 = 0.0, 2.0, 1.0
f(x) = x^2

function compute_pv_integral(n)
    nodes, weights, n_actual = cpvi(a, b, x0, n)
    return sum(weights[i] * f(nodes[i]) / (nodes[i] - x0) for i in 1:n_actual)
end

println("节点数收敛性测试：")
for n in [50, 100, 200, 500, 1000]
    result = compute_pv_integral(n)
    println("n = $n: 结果 = $result")
end
```

**示例 6：与解析解比较**

```julia
# 对于某些简单情况，可以验证数值结果
# P.V. ∫₋₁¹ 1/(x-0) dx = 0 （奇函数）

a, b, x0 = -1.0, 1.0, 0.0
nodes, weights, n_actual = cpvi(a, b, x0, 500)

f(x) = 1.0
result = sum(weights[i] * f(nodes[i]) / (nodes[i] - x0) for i in 1:n_actual)

println("数值结果: ", result)
println("解析解: 0.0")
println("绝对误差: ", abs(result))
```

#### 精度特性

1. **梯形法则精度**：局部误差 O(h³)，全局误差 O(h²)
2. **节点密度**：h = (b-a)/(n-1)，增加 n 可提高精度
3. **奇点处理**：通过 0.001h 的间隙避免数值不稳定
4. **对称检验**：对于对称问题，可验证奇函数积分 ≈ 0

#### 注意事项

1. **奇点位置验证**
   - 必须保证 a < x₀ < b
   - 奇点过于接近端点可能导致精度下降

2. **节点数调整**
   - 实际节点数 n_actual 通常接近但不完全等于输入 n
   - 这是由于在奇点附近的离散化调整

3. **数值稳定性**
   - 避免在靠近奇点的节点处直接求值
   - 0.001h 间隙是经验值，对大多数问题有效

4. **被积函数光滑性**
   - 除了 1/(x-x₀) 奇点外，f(x) 应该光滑
   - 额外的奇点或不连续性需要特殊处理

5. **物理单位**
   - 确保所有输入量的单位一致（fm⁻¹）
   - 积分结果的单位取决于被积函数

#### 误差分析

**主要误差来源：**

1. **截断误差**：梯形法则 O(h²)
2. **奇点间隙**：0.001h 引入的 O(h) 误差
3. **舍入误差**：浮点运算累积误差

**改进方法：**

- 增加节点数 n
- 使用高阶积分方法（Simpson 法则等）
- 自适应网格细化

#### 异常处理

函数在以下情况会抛出 `ArgumentError`：

1. **节点数无效**: `n ≤ 1`
   ```julia
   cpvi(0.0, 1.0, 0.5, 1)  # ArgumentError
   ```

2. **区间无效**: `a ≥ b`
   ```julia
   cpvi(1.0, 0.0, 0.5, 100)  # ArgumentError
   ```

3. **奇点位置无效**: `x₀ ≤ a` 或 `x₀ ≥ b`
   ```julia
   cpvi(0.0, 1.0, 1.5, 100)  # ArgumentError: 奇点在区间外
   cpvi(0.0, 1.0, 0.0, 100)  # ArgumentError: 奇点在端点
   ```

#### 性能建议

1. **节点数选择**
   - 光滑函数：n = 100-500 通常足够
   - 快速变化函数：增加到 n = 1000-5000
   - 先用小 n 测试，观察收敛性

2. **重复计算优化**
   - 对于固定的 a, b, x₀, n，预计算节点和权重
   - 存储并重用于不同的被积函数

3. **向量化**
   - 利用 Julia 的向量化操作加速计算
   - 使用生成器表达式节省内存

#### 与 GaussLegendre 的比较

| 特性 | CauchyPV | GaussLegendre |
|------|----------|---------------|
| 适用场景 | 含奇点积分 | 光滑函数积分 |
| 精度 | O(h²) | 指数收敛 |
| 节点分布 | 等距 + 奇点避让 | 非等距最优 |
| 计算成本 | 低 | 中等 |
| 实现复杂度 | 简单 | 中等 |

---

## 理论背景

### 柯西主值的物理意义

在物理学中，柯西主值积分常见于：

1. **色散关系（Dispersion Relations）**
   - Kramers-Kronig 关系
   - 光学响应函数
   - 偏振张量的实部和虚部关系

2. **格林函数（Green's Functions）**
   - 传播子的极点贡献
   - 因果性条件

3. **多体理论（Many-Body Theory）**
   - 自能的计算
   - 准粒子谱函数

4. **散射理论（Scattering Theory）**
   - 光学定理
   - 相移分析

### 数学性质

1. **线性性**：$P.V. \int [af + bg] = a \cdot P.V.\int f + b \cdot P.V.\int g$

2. **对称性**：对于对称区间 $[-L, L]$ 和 $x_0 = 0$，
   奇函数的主值积分为零

3. **Plemelj 公式**：
   $$
   \lim_{\varepsilon \to 0^+} \int \frac{f(x)}{x - x_0 \pm i\varepsilon} dx = P.V. \int \frac{f(x)}{x - x_0} dx \mp i\pi f(x_0)
   $$

---

## 常见问题

**Q: 为什么实际节点数与输入不完全一致？**

A: 由于需要在奇点附近调整节点位置，实际节点数根据步长和奇点位置动态计算，通常与输入 n 接近但不完全相同。

**Q: 0.001h 的间隙是否可以调整？**

A: 可以修改源代码中的系数。较小的间隙提高精度但可能降低数值稳定性；较大的间隙更稳定但精度下降。0.001 是经验平衡值。

**Q: 如何判断积分是否收敛？**

A: 逐步增加节点数，如果结果变化小于所需精度，则认为已收敛。也可以用奇函数验证（结果应为零）。

**Q: 可以处理多个奇点吗？**

A: 当前实现只支持单一奇点。多个奇点需要分段积分或使用更复杂的方法。

**Q: 与解析主值积分的关系？**

A: 本方法是数值近似。对于可以解析计算的情况（如多项式除以一次因子），解析方法更精确。

---

## 参考文献

1. **Press, W. H., et al. (2007)**  
   *Numerical Recipes: The Art of Scientific Computing* (3rd ed.)  
   Cambridge University Press - Chapter on Singular Integrals

2. **Davis, P. J., & Rabinowitz, P. (1984)**  
   *Methods of Numerical Integration* (2nd ed.)  
   Academic Press - Chapter on Cauchy Principal Values

3. **Fetter, A. L., & Walecka, J. D. (2003)**  
   *Quantum Theory of Many-Particle Systems*  
   Dover Publications - Dispersion Relations

4. **Mahan, G. D. (2000)**  
   *Many-Particle Physics* (3rd ed.)  
   Springer - Green's Functions and Self-Energy

---

## 更新日志

### v0.1.0 (2025-10-14)
- 初始版本
- 实现 `cpvi` 函数
- 基于复化梯形法则
- 完整的单位约定和文档
- 从 Fortran 代码移植并优化
