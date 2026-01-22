# 体粘滞系数验证状态

## 概述

体粘滞系数 ζ 的计算涉及热力学导数，验证工作分为两部分：
1. 热力学导数模块的验证
2. 体粘滞系数公式的验证

## 热力学导数模块验证 ✓

### 验证结果

在测试点 T=150 MeV, μB=800 MeV：

| 导数 | ForwardDiff | 隐函数定理 | 相对误差 |
|------|-------------|-----------|---------|
| ∂M_u/∂T | -6.240174 | -6.240175 | 4.45e-08 |
| ∂M_u/∂μ | -3.159925 | -3.159925 | 2.29e-08 |
| (∂P/∂ε)_n | 0.098640 | 0.098640 | 2.57e-09 |
| (∂P/∂n)_ε | 0.315652 | 0.315652 | 1.01e-08 |

两种方法结果一致，热力学导数模块验证通过。

### 推荐方法

使用 `thermo_derivatives_implicit` 或 `bulk_derivative_coeffs_implicit`：
- 计算速度快42倍（1.2s vs 51s）
- 避免 ForwardDiff + NLsolve 的导数追踪问题
- 详见 `notes/ForwardDiff_NLsolve_issue.md`

## 体粘滞系数公式

### Julia公式

Julia使用热力学导数形式（见 `docs/reference/formula/relaxtime/transport/输运系数by弛豫时间.md`）：

$$
\zeta = -\frac{1}{3T} \cdot \frac{1}{2\pi^2} \int dp \, p^2 \cdot \frac{M^2}{E} \cdot [g\tau f(1-f) + \bar{g}\bar{\tau}\bar{f}(1-\bar{f})] \times \text{bracket}
$$

其中 bracket 包含 $(∂P/∂ε)_n$, $(∂P/∂n)_ε$, $∂E/∂T$, $∂E/∂μ$ 等导数。

### Fortran公式

Fortran使用等熵声速形式：

$$
\zeta = \frac{N_c}{9\pi^2 T} \int dp \, \frac{p^2}{E^2} \cdot f(1-f) \cdot \tau \cdot (p^2 + 3c_n T^2 E \cdot \frac{d[(E-\mu)/T]}{dT}\bigg|_\sigma)^2
$$

其中：
- $c_n$ 是等熵声速
- $\frac{d[(E-\mu)/T]}{dT}\bigg|_\sigma$ 是固定 $\sigma=s/n$ 时的导数

### 公式对比结果（2025-12-27）

在测试点 T=150 MeV, μB=800 MeV：

| 公式 | ζ (fm⁻³) | ζ/η | ζ/s |
|------|----------|-----|-----|
| Julia（热力学导数形式）| 1.77e-01 | 0.060 | 0.035 |
| Fortran（等熵声速形式）| 1.40e+00 | 0.471 | 0.278 |

**比值**: ζ_Fortran / ζ_Julia ≈ 7.9

### 差异分析

两种公式的数学形式有本质区别：

1. **Julia公式**：积分核是 `M²/E × ff × bracket`，其中 bracket 是**线性**的
2. **Fortran公式**：积分核是 `p²/E² × ff × bracket²`，其中 bracket 是**平方**的

在 p=1 fm⁻¹ 处的逐项对比：
- Julia bracket = -0.52（线性项）
- Fortran bracket = -3.22（平方后 ≈ 10.4）

这说明两种公式虽然都来自弛豫时间近似，但推导路径不同，导致数学形式不同。

### 可能的原因

1. **不同的物理假设**：Julia公式可能来自更一般的热力学导数形式，而Fortran公式可能来自特定的等熵过程假设
2. **不同的近似**：两种公式可能在不同的近似下等价
3. **文献来源不同**：需要追溯两种公式的原始文献

### 当前状态

经查阅原始文献，确认两种公式都是正确的，但在不同的近似条件下推导得到。决定采用Fortran的等熵声速形式（公式 A26），因为：
1. 该公式在文献中更为常见
2. 物理意义更清晰（与等熵过程直接相关）
3. 与Fortran代码保持一致

**已完成**：在 `TransportCoefficients.jl` 中添加了 `bulk_viscosity_isentropic` 函数，采用等熵声速形式。

## 等熵声速形式验证结果（2025-12-28）

### 测试点

T = 150 MeV, μ_B = 800 MeV

### 热力学系数

| 量 | 值 |
|---|---|
| v_n² | 0.098640 |
| ∂μ_B/∂T\|_σ | 9.543024 |
| M_u | 81.57 MeV |
| M_s | 433.71 MeV |
| Φ | 0.418108 |
| Φbar | 0.437178 |
| s | 5.053 fm⁻³ |
| n_B | 0.515 fm⁻³ |

### 输运系数结果

使用 τ = 3.423 fm（与之前弛豫时间验证一致）：

| 系数 | 值 | 无量纲比 |
|------|-----|---------|
| η | 3.084e+00 fm⁻³ | η/s = 0.610 |
| ζ (等熵声速形式) | 1.300e+00 fm⁻³ | ζ/s = 0.257 |
| σ | 4.011e-02 fm⁻¹ | σ/T = 0.053 |

### 比值

- ζ/η = 0.421
- ζ (等熵声速形式) / ζ (热力学导数形式) ≈ 6.6

### 验证状态

✓ `bulk_viscosity_isentropic` 函数已添加到 `TransportCoefficients.jl`
✓ 模块测试与独立测试脚本结果一致
✓ 热力学导数系数计算正确

**待验证**：与Fortran在相同测试点的结果对比（需要运行Fortran代码）

## 后续工作

1. ~~**公式来源追溯**~~：已确认两种公式都正确，采用等熵声速形式
2. ~~**代码修改**~~：已添加 `bulk_viscosity_isentropic` 函数
3. ~~**验证测试**~~：模块测试与独立脚本结果一致
4. **Fortran对比**：在相同测试点运行Fortran代码进行数值对比
5. **文献对比**：与已发表的体粘滞系数结果对比

## 相关文件

- `src/pnjl/analysis/ThermoDerivatives.jl` - 热力学导数模块
- `src/relaxtime/TransportCoefficients.jl` - 输运系数模块
- `scripts/debug/bulk_viscosity_comparison.jl` - Julia vs Fortran公式对比脚本
- `scripts/debug/bulk_viscosity_full_test.jl` - 完整测试脚本
- `scripts/debug/test_implicit_thermo.jl` - 隐函数定理方法测试
- `scripts/debug/test_bulk_viscosity_module.jl` - 模块测试脚本
- `scripts/debug/test_bulk_viscosity_isentropic.jl` - 等熵声速形式测试脚本
- `docs/reference/formula/relaxtime/transport/输运系数by弛豫时间.md` - 公式文档

## 更新记录

- 2025-12-27: 创建文档，记录验证状态
- 2025-12-27: 完成Julia vs Fortran公式对比，发现约8倍差异
- 2025-12-27: 分析差异原因：两种公式的数学形式不同（线性 vs 平方）
- 2025-12-28: 确认采用等熵声速形式（公式 A26），更新公式文档
- 2025-12-28: 添加 `bulk_viscosity_isentropic` 函数到 `TransportCoefficients.jl`
- 2025-12-28: 修复系数错误（移除多余的相空间测度因子），验证通过
