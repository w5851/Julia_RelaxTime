# RelaxationTime 模块 API

弛豫时间计算，基于平均散射率表达式。

## 导出符号
- `REQUIRED_PROCESSES`：需要的散射过程符号列表
- `compute_average_rates(quark_params, thermo_params, K_coeffs; ...)`
- `relaxation_rates(densities, rates)`
- `relaxation_times(quark_params, thermo_params, K_coeffs; densities, ...)`

## 数学定义
$$\tau_i^{-1} = \sum_j \rho_j \; \bar{w}_{ij}$$

## 默认行为

当不传入 `p_grid`/`p_w` 时：
- 动量积分范围：`[0, Λ]`（Λ ≈ 3.05 fm⁻¹）
- σ(s) 缓存范围：基于 Λ 截断的 w0cdf 设计
- 数密度积分：半无穷 `[0, ∞)`（物理要求，由调用方提供）

这确保了动量积分范围与 σ(s) 缓存范围的一致性。

## 参数说明

- `quark_params`：夸克参数 `(m, μ, A)`
- `thermo_params`：热力学参数 `(T, Φ, Φbar, ξ)`
- `K_coeffs`：有效耦合系数
- `densities`：数密度（**必须使用半无穷积分计算**）
- `p_grid/p_w`：可选，自定义动量积分节点
- `sigma_cutoff`：σ(s) 有效范围的动量截断（默认 Λ）


## 典型用法

```julia
using .RelaxationTime
using .GaussLegendre: gauleg
using .Constants_PNJL: Λ_inv_fm

# 方式1：使用默认行为（推荐）
result = relaxation_times(quark_params, thermo_params, K_coeffs;
    densities=densities)

# 方式2：有限截断积分（与 C++ 一致）
p_grid, p_w = gauleg(0.0, 15.0, 64)
result = relaxation_times(quark_params, thermo_params, K_coeffs;
    densities=densities,
    p_grid=p_grid, p_w=p_w,
    sigma_cutoff=Λ_inv_fm)

println(result.tau)      # 各味的 τ
println(result.tau_inv)  # 各味的 τ^-1
```

## 与 C++ 实现对比

| 参数 | Julia 默认 | C++ 实现 | 匹配方式 |
|------|-----------|----------|----------|
| 动量积分范围 | [0, Λ] | [0, 15] fm⁻¹ | 传入 `p_grid=gauleg(0,15,n)` |
| σ(s) 缓存范围 | 基于 Λ 截断 | 基于 Λ 截断 | 已一致 |
| 数密度积分 | [0, ∞) | [0, ∞) | 已一致 |

使用 `p_grid=gauleg(0,15,n)` + `sigma_cutoff=Λ_inv_fm` 后，与 C++ 结果误差 < 0.3%。

## 物理说明

- **数密度积分范围**：必须使用半无穷积分 [0, ∞)，这是物理要求
- **σ(s) 有效范围**：由 PNJL 模型的动量截断 Λ 决定
- **重夸克注意**：当夸克质量 m 接近 Λ 时（如低温下的 s 夸克），相空间被限制
