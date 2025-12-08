# ThermoDerivatives (PNJL)

基于自动微分的热力学偏导与体粘滞系数组合导数。求导时会重新求解能隙方程，保证量满足平衡约束。

## 主要接口

- `solve_equilibrium_mu(T_mev, mu_mev; xi=0.0, seed_state, p_num, t_num)`
  - 返回在给定 `(T, μ, ξ)` 下的平衡解：`pressure`, `energy`(能量密度 ε), `rho`(净重子数密度 = Σρ_i/3)，`rho_norm`, `entropy`，以及 `x_state`, `mu_vec` 等元数据。

- `thermo_derivatives(T_mev, mu_mev; xi=0.0, ...)`
  - 返回基础量及一阶导数：`dP_dT`, `dP_dmu`, `dEpsilon_dT`, `dEpsilon_dmu`, `dn_dT`, `dn_dmu`，以及组合导数 `dP_depsilon_n`, `dP_dn_epsilon`（用于体粘滞公式）。

- `dP_dT/dP_dmu`
  - 压强对 `T`/`μ` 的高阶导数，`order` 默认 1。

- `dEpsilon_dT/dEpsilon_dmu`
  - 能量密度 ε 对 `T`/`μ` 的高阶导数（旧版的 `dE_*` 语义）。

- `quasiparticle_energy(T_mev, mu_mev, p_fm; flavor=1, xi=0.0, ...)`
  - 取平衡解质量后返回各向同性色散关系 `E = sqrt(p^2 + m^2)`。

- `dE_dT/dE_dmu (T_mev, mu_mev, p_fm; m, dM_dT|dM_dmu)`
  - 单粒子能量对 `T`/`μ` 的导数，需传入已知的质量 `m` 与其导数 `dM_dT`/`dM_dmu`，内部仅做链式法则 `m/E·dM`，不再重新求解。

- `bulk_derivative_coeffs(T_mev, mu_mev; xi=0.0, ...)`
  - 便捷返回 `(∂P/∂ε)_n` 与 `(∂P/∂n)_ε` 两个组合导数。

## 使用示例

```julia
using .PNJL.ThermoDerivatives

T = 130.0; mu = 320.0; xi = 0.2
# 宏观导数（体粘滞用）
k = thermo_derivatives(T, mu; xi=xi)
println(k.dP_depsilon_n, k.dP_dn_epsilon)

# 单粒子能量及其导数（动量单位：fm⁻¹）
p = 0.5
E = quasiparticle_energy(T, mu, p; xi=xi)
E_T = dE_dT(T, mu, p; xi=xi)
```

## 性能与注意事项
- 求导会多次调用非线性方程求解器，可适当降低 `p_num/t_num` 或复用 `seed_state` 以减少开销。
- `dP_depsilon_n`/`dP_dn_epsilon` 分母退化时返回 `NaN`，需检查物理解或调整区间。
- 单粒子能量使用各向同性色散关系；各向异性只进入分布函数，不影响色散。 
