# AverageScatteringRate.jl

平均散射率计算模块（实验性，支持各向异性）。

- 入口：`average_scattering_rate(process, quark_params, thermo_params, K_coeffs; p_nodes=32, angle_nodes=4, phi_nodes=4, cs_cache=CrossSectionCache(process), n_sigma_points=32)`
- 依赖：`TotalCrossSection.total_cross_section`（可通过 `CrossSectionCache` 预计算+插值），`quark_distribution(_aniso)`/`antiquark_distribution(_aniso)` 作为分布函数。
- 默认节点：动量 32 点、角度 4 点；Gauss-Legendre 节点来自 `integration/GaussLegendre.jl`。
- 各向异性：`ξ=0` 时退化为各向同性；`ξ≠0` 时分布调用 `quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)`，角度取自积分节点。
- 预期用途：为弛豫时间近似提供过程级平均散射率 \(\omega_{ij}\)。

## 核心公式

- 平均散射率（利用对称性后的积分）：

$$
\omega_{ij} = \frac{d_q^2}{4\pi^5\rho_i\rho_j} \int dp_i\,p_i^2 \int dp_j\,p_j^2 \int_0^1 d\cos\theta_i \int_0^1 d\cos\theta_j \int_0^\pi d\phi \, f_i f_j \, v_{\text{rel}} \, \sigma_{ij\to cd}(s, T, \mu_q)
$$

- 数密度（各向异性可选）：
$$
\rho_i = \frac{d_q}{2\pi^2} \int_0^\Lambda dp\,p^2 \int_0^1 d\cos\theta\; f_i(p,\cos\theta)
$$

- 运动学：
  - \(s = m_i^2 + m_j^2 + 2(E_iE_j - p_i p_j \cos\Theta)\)，\(\cos\Theta = \cos\theta_i\cos\theta_j + \sin\theta_i\sin\theta_j \cos\phi\)
  - 相对速度 \(v_{\text{rel}} = \sqrt{(p_i\cdot p_j)^2 - m_i^2 m_j^2}/(E_iE_j)\)

## 主要 API

- `CrossSectionCache(process::Symbol)`: 创建截面缓存；内部存储 `(s_vals, sigma_vals)`，同时持有基于拟合的 σ(s) 近似（见下）。
- `precompute_cross_section!(cache, s_grid, quark_params, thermo_params, K_coeffs; n_points=32)`: 预计算并填充截面表。
- `average_scattering_rate(process, quark_params, thermo_params, K_coeffs; p_nodes=32, angle_nodes=4, phi_nodes=4, cs_cache=CrossSectionCache(process), n_sigma_points=32)`: 计算给定过程的平均散射率。

### σ(s) 拟合与调用流程

- 首次调用 `get_sigma` 会估算积分所需的 s 范围：阈值 `max((m_i+m_j)^2, (m_c+m_d)^2)` 到 `p_max=Λ_inv_fm` 推出的上界并乘以 1.05 余量。
- 在该范围内取约 24 个左端偏置节点（指数偏置，λ=4），逐点调用一次 `total_cross_section` 真实计算并拟合：
  - 初态更重（m_i+m_j≥m_c+m_d）：用 `sqrt` 型基函数 `A/√(s-s0)+B(s-s0)+C(s-s0)^2+D` 最小二乘拟合。
  - 初态更轻：用分段线性插值（外推到边界）。
- 随后所有落在拟合范围内的 σ 查询都直接用拟合值，积分时不再触发真实计算；仅当 s 超出拟合范围时才回退单点精确计算并缓存。
- 拟合成本：每个 process 首次约 24 次 `total_cross_section` 调用；之后常数时间查表。

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
```

## 注意事项

- 截面计算仍依赖 `TotalCrossSection` 链路（标记“修复中”），用于研究或诊断；生产前需校验。
- 若算力不足，可降低 `p_nodes`、`angle_nodes`、`phi_nodes` 或减少 `n_sigma_points`；或预先将截面插值表存盘再加载。
- 近阈值处若 \(\rho_i\) 或 \(\rho_j\) 为 0，将返回 0 以避免数值不稳定。
- 截断 \(\Lambda\) 取自 `Constants_PNJL.Λ_inv_fm`。
