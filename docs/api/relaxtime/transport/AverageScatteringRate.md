# AverageScatteringRate.jl

平均散射率计算模块（支持各向异性）。

## 概述

- 入口：`average_scattering_rate(process, quark_params, thermo_params, K_coeffs; ...)`
- 依赖：`TotalCrossSection.total_cross_section`（通过 **w0cdf 取点 + PCHIP 插值**预计算成 σ(s) 表）
- 默认节点：动量 20 点、角度 4 点、方位角 8 点
- 各向异性：`ξ=0` 时退化为各向同性；`ξ≠0` 时使用 `quark_distribution_aniso`

## 动量积分模式

### 默认行为（推荐）

通过 `RelaxationTime.relaxation_times` 调用时，默认截断、数密度积分与 σ(s) 缓存口径应遵循 `RelaxationTime.md` 的统一说明；本页只补充平均散射率这一单过程积分层的细节。

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

## 截面缓存指纹（fingerprint）

`CrossSectionCache` 会在通过 `precompute_cross_section!` 或 `build_w0cdf_pchip_cache` 生成 σ(s) 表时写入 `fingerprint` 元数据。该元数据记录过程名、夸克质量/化学势/A、热力学参数、`K_coeffs`、σ 积分点数、阈值渐近正则化参数，以及 σ 网格摘要（长度、端点和哈希）。w0cdf 构建路径还记录 `N`、设计节点数、`p_cutoff` 与 `scale`，其中 `N` 与设计节点数用于溯源，不要求等于 `average_scattering_rate` 的当前默认值。

当调用方把已填充的 `cs_cache` 传给 `average_scattering_rate` 时，默认会校验当前参数上下文与缓存指纹。若过程、物理参数、`K_coeffs`、`n_sigma_points`、有效 `threshold_subtraction` 相关参数，或 w0cdf 的 `p_cutoff`/`scale` 口径不一致，会抛出 `ArgumentError`，避免静默复用不同参数点生成的 σ(s) 表。自定义 `N` 或设计节点数生成的 w0cdf 缓存可以复用，但缓存自身的 s-grid 摘要必须与已存表一致。

兼容策略：

- `CrossSectionCache(process, s_vals, sigma_vals)` 和手工 `insert_sigma!` 填充的旧缓存默认没有 `fingerprint`；调用时允许复用，但会发出 warning。
- `load_cross_section_caches_from_dir` 读取的两列 CSV 旧表同样没有 `fingerprint`；默认允许但 warning。
- 若希望严格拒绝无指纹缓存，传入 `require_cache_fingerprint=true`。
- 多个参数点之间只应复用同一生成指纹下的缓存；跨 `quark_params`、`thermo_params`、`K_coeffs`、`n_sigma_points`、有效 `threshold_subtraction`、`p_cutoff` 或 `scale` 复用时，应重新构建缓存。

## 阈值渐近正则化（threshold_subtraction）

- 新增字段与参数（2026-01-29 提交引入）：
    - `CrossSectionCache` 字段：`asym_enabled::Bool`, `asym_s0::Float64`, `asym_A::Float64`, `asym_requested::Bool`。
    - `precompute_cross_section!(cache, s_grid, quark_params, thermo_params, K_coeffs; n_points, threshold_subtraction=false, asym_window=0.6, asym_fit_min_points=8, asym_extra_points=10)`：新增关键字参数以启用阈值附近的额外采样与拟合，用于检测并减去形如 $A/\sqrt{s-s_0}$ 的渐近项，再把已正则化的值缓存起来。
    - `build_w0cdf_pchip_cache(...; threshold_subtraction=false, asym_window=0.6, asym_fit_min_points=8, asym_extra_points=10)`：构建缓存时会采用与 `average_scattering_rate` 相同的自动阈值减法判定；若初态质量平方和大于末态质量平方和，实际写入 fingerprint 的 `threshold_subtraction` 会被解析为 `true`，并传给 `precompute_cross_section!`。

- 查询与返回语义：
    - 当查询点 `s` 在缓存外且 `asym_enabled` 时，`get_sigma(cache,s)` 会返回仅由解析渐近项给出的近似值 $A/\sqrt{s-s_0}$（若适用）。
    - 当查询点在缓存内且缓存为已正则化的 `sigma_vals = raw - asym`，`get_sigma` 会把插值结果与解析渐近项相加后返回最终的 σ(s)。

- 其他行为：
    - 若在构建缓存时为阈值附近额外采样（`asym_extra_points`），这些点会加入缓存以改善插值质量；默认不开启以保证向后兼容性。


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

与 C++ 的默认口径对齐方式与 `RelaxationTime.md` 相同：

- 通过 `p_grid=gauleg(0,15,n)` 指定有限截断积分
- 通过 `sigma_cutoff=Λ_inv_fm` 固定 σ(s) 有效范围
- 数密度仍保持半无穷积分

如果你关心的是 τ 聚合层面的默认行为与 C++ 总体对齐，请直接看 `RelaxationTime.md`；本页不再重复完整对比表。

## 注意事项

- 截断 Λ 取自 `Constants_PNJL.Λ_inv_fm`（约 3.05 fm⁻¹）
- 数密度积分**必须**使用半无穷积分，这是物理要求
- 近阈值处若 ρ_i 或 ρ_j 为 0，返回 0 以避免数值不稳定
- 当质心能量 s 超出 σ(s) 缓存覆盖区间时，对应贡献按 0 处理（不会把 s 钳制到边界）
