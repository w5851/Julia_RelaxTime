# AverageScatteringRate.jl

平均散射率计算模块（实验性，支持各向异性）。

- 入口：`average_scattering_rate(process, quark_params, thermo_params, K_coeffs; p_nodes=20, angle_nodes=4, phi_nodes=8, cs_cache=nothing, n_sigma_points=32, ...)`
- 依赖：`TotalCrossSection.total_cross_section`（通过 **w0cdf 取点 + PCHIP 插值**预计算成 σ(s) 表），`quark_distribution(_aniso)`/`antiquark_distribution(_aniso)` 作为分布函数。
- 默认节点：动量 20 点、角度 4 点、方位角 8 点；可按需调大节点数，或直接传入预计算的 Gauss-Legendre 节点/权重避免重复生成。
- 各向异性：`ξ=0` 时退化为各向同性；`ξ≠0` 时分布调用 `quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)`，角度取自积分节点。
- 预期用途：为弛豫时间近似提供过程级平均散射率 \(\omega_{ij}\)。

## 核心公式

- 平均散射率（利用对称性后的积分）：

$$
\omega_{ij} = \frac{d_q^2}{4\pi^5\rho_i\rho_j} \int dp_i\,p_i^2 \int dp_j\,p_j^2 \int_0^1 d\cos\theta_i \int_0^1 d\cos\theta_j \int_0^\pi d\phi \, f_i f_j \, v_{\text{rel}} \, \sigma_{ij\to cd}(s, T, \mu_q)
ω = average_scattering_rate(:uu_to_uu, quark_params, thermo_params, K_coeffs;
  cs_cache=cache, p_nodes=6, angle_nodes=2, phi_nodes=4)
- 数密度（各向异性可选）：
# 也可以复用外部预生成的节点/权重（示例）
$$
\rho_i = \frac{d_q}{2\pi^2} \int_0^\infty dp\,p^2 \int_0^1 d\cos\theta\; f_i(p,\cos\theta)
$

**重要说明**：数密度积分应使用半无穷积分范围 $[0, \infty)$，而非有限截断 $[0, \Lambda]$。$

  cs_cache=cache, p_grid=p_grid, p_w=p_w, cos_grid=cos_grid, cos_w=cos_w, phi_grid=phi_grid, phi_w=phi_w)
- 运动学：
  - \(s = m_i^2 + m_j^2 + 2(E_iE_j - p_i p_j \cos\Theta)\)，\(\cos\Theta = \cos\theta_i\cos\theta_j + \sin\theta_i\sin\theta_j \cos\phi\)
  - 相对速度 \(v_{\text{rel}} = \sqrt{(p_i\cdot p_j)^2 - m_i^2 m_j^2}/(E_iE_j)\)

## 主要 API

- `CrossSectionCache(process::Symbol)`: 创建空截面缓存；内部存储 `(s_vals, sigma_vals)`。
- `precompute_cross_section!(cache, s_grid, quark_params, thermo_params, K_coeffs; n_points=32)`: 在给定 `s_grid` 上预计算并填充截面表。
- `build_w0cdf_pchip_cache(process, quark_params, thermo_params, K_coeffs; N=60, design_p_nodes=14, design_angle_nodes=4, design_phi_nodes=8, ...)`: 生成生产默认的 σ(s) 缓存（w0cdf + PCHIP）。
- `average_scattering_rate(process, quark_params, thermo_params, K_coeffs; cs_cache=nothing, ...)`: 计算给定过程的平均散射率；若 `cs_cache === nothing`，会自动构建 w0cdf+PCHIP 的 σ(s) 缓存。

### σ(s) 缓存策略（固定）

- 生产默认：在“ω 的权重 $w_0$”下对 $\sqrt{s}$ 构造加权 CDF，并用分位点选取 $N$ 个采样点（默认 $N=60$）。
- 对这些点逐点预计算 `total_cross_section(process, s, ...)` 得到 σ(s) 表。
- 查询时只做 PCHIP 插值（并对端点做钳制），避免积分过程中触发昂贵的 σ(s) 精确计算。

## 使用示例

```julia
using .AverageScatteringRate

quark_params = (
    m = (u=1.52, d=1.52, s=2.50),
    μ = (u=0.3, d=0.3, s=0.0)
)
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
K_coeffs = (K_σπ=2.0, K_σK=2.0, K_σ=3.0, K_δπ=1.5, K_δK=1.5)

cache = CrossSectionCache(:uu_to_uu)
precompute_cross_section!(cache, range(10.0, 40.0, length=5), quark_params, thermo_params, K_coeffs)

ω = average_scattering_rate(:uu_to_uu, quark_params, thermo_params, K_coeffs;
    cs_cache=cache, p_nodes=32, angle_nodes=4, phi_nodes=4)

# 或者使用生产默认（自动构建 w0cdf+PCHIP cache）：
ω2 = average_scattering_rate(:uu_to_uu, quark_params, thermo_params, K_coeffs)
```

## 注意事项

- 截面计算仍依赖 `TotalCrossSection` 链路（标记“修复中”），用于研究或诊断；生产前需校验。
- 若算力不足，可降低 `p_nodes`、`angle_nodes`、`phi_nodes` 或减少 `n_sigma_points`；或预先将 σ(s) 表存盘再加载。
- 近阈值处若 \(\rho_i\) 或 \(\rho_j\) 为 0，将返回 0 以避免数值不稳定。
- 截断 \(\Lambda\) 取自 `Constants_PNJL.Λ_inv_fm`。
