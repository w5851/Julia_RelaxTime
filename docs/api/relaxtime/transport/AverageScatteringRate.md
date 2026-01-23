# AverageScatteringRate.jl

平均散射率计算模块（支持各向异性）。

## 概述

- 入口：`average_scattering_rate(process, quark_params, thermo_params, K_coeffs; ...)`
- 依赖：`TotalCrossSection.total_cross_section`（通过 **w0cdf 取点 + PCHIP 插值**预计算成 σ(s) 表）
- 默认节点：动量 20 点、角度 4 点、方位角 8 点
- 各向异性：`ξ=0` 时退化为各向同性；`ξ≠0` 时使用 `quark_distribution_aniso`

## 动量积分模式

### 默认行为（推荐）

通过 `RelaxationTime.relaxation_times` 调用时：
- 动量积分范围：`[0, Λ]`（Λ ≈ 3.05 fm⁻¹）
- σ(s) 缓存范围：基于 Λ 截断的 w0cdf 设计
- 数密度积分：半无穷 `[0, ∞)`（物理要求）

### 有限截断模式（与 C++ 一致）

传入 `p_grid`/`p_w` 和 `sigma_cutoff` 参数：
```julia
p_grid, p_w = gauleg(0.0, 15.0, 64)
ω = average_scattering_rate(:uu_to_uu, quark_params, thermo_params, K_coeffs;
    p_grid=p_grid, p_w=p_w,
    sigma_cutoff=Λ_inv_fm)
```

## 核心公式

详见 `docs/reference/formula/relaxtime/scattering/AverageScatteringRate_FromCrossSection.md`。

平均散射率：
$$\omega_{ij} = \frac{d_q^2}{32\pi^5\rho_i\rho_j} \int d^3p_i\,d^3p_j\, f_i f_j \, v_{\text{rel}} \, \sigma_{ij}(s)$$

数密度（**必须**使用半无穷积分）：
$$\rho_i = \frac{d_q}{2\pi^2} \int_0^\infty dp\,p^2 \int_0^1 d\cos\theta\; f_i(p,\cos\theta)$$


## 主要 API

- `CrossSectionCache(process::Symbol)`: 创建空截面缓存
- `precompute_cross_section!(cache, s_grid, ...)`: 预计算 σ(s) 表
- `design_w0cdf_s_grid(process, ...; p_cutoff=nothing)`: 设计 σ(s) 采样网格
  - `p_cutoff=nothing`: 半无穷积分权重设计
  - `p_cutoff=Λ_inv_fm`: 有限截断权重设计（**推荐**）
- `build_w0cdf_pchip_cache(process, ...; p_cutoff=nothing)`: 构建 σ(s) 缓存
- `average_scattering_rate(process, ...; sigma_cutoff=nothing)`: 计算平均散射率
- `number_density(flavor, ...)`: 计算数密度（始终使用半无穷积分）

## 使用示例

```julia
using .AverageScatteringRate
using .Constants_PNJL: Λ_inv_fm

# 方式1：使用默认行为（Λ 截断，推荐）
ω = average_scattering_rate(:uu_to_uu, quark_params, thermo_params, K_coeffs;
    sigma_cutoff=Λ_inv_fm)

# 方式2：有限截断积分（与 C++ 一致）
p_grid, p_w = gauleg(0.0, 15.0, 64)
ω = average_scattering_rate(:uu_to_uu, quark_params, thermo_params, K_coeffs;
    p_grid=p_grid, p_w=p_w,
    sigma_cutoff=Λ_inv_fm)
```

## 与 C++ 实现对比

| 参数 | Julia 默认 | C++ 实现 | 匹配方式 |
|------|-----------|----------|----------|
| 动量积分范围 | [0, Λ] | [0, 15] fm⁻¹ | 传入 `p_grid=gauleg(0,15,n)` |
| σ(s) 有效范围 | 基于 Λ 截断 | 基于 Λ 截断 | 已一致 |
| 数密度积分 | [0, ∞) | [0, ∞) | 已一致 |

使用 `p_grid=gauleg(0,15,n)` + `sigma_cutoff=Λ_inv_fm` 后，与 C++ 结果误差 < 0.3%。

## 注意事项

- 截断 Λ 取自 `Constants_PNJL.Λ_inv_fm`（约 3.05 fm⁻¹）
- 数密度积分**必须**使用半无穷积分，这是物理要求
- 近阈值处若 ρ_i 或 ρ_j 为 0，返回 0 以避免数值不稳定
- 当质心能量 s 超出 σ(s) 缓存覆盖区间时，对应贡献按 0 处理（不会把 s 钳制到边界）
