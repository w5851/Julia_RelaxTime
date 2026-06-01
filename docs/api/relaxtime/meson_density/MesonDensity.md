# MesonDensity 模块 API 文档

## 模块概述

`MesonDensity` 模块提供当前介子数密度主线的最小实现入口，当前覆盖四层最小数值 helper：

- `π/K` 聚合通道默认简并因子
- 玻色分布函数
- 稳定粒子极限数密度
- reduced strict BW 数密度
- Stage2 q 依赖复极点 strict BW 数密度
- `K/π` 比值与温度扫描
- Phase E3 当前最小口径下的 `π/K` 相移双积分 helper
- BU2020 审计所需的 real-axis 分支、相位口径与 Bose-domain policy metadata

后续 BU / BW / 各向异性扩展将在该模块基础上继续演化。

## 单位约定

与项目其他 `relaxtime` 模块一致，内部计算使用自然单位制：

- `mass`, `T`, `μ`, `q`：`fm^-1`
- 返回的数密度：`fm^-3`

## API 参考

### `meson_degeneracy(meson; charge_resolved=false)`

返回当前主线的 `π/K` 简并因子。

- 聚合通道：`d_π = 3`、`d_K = 4`
- 电荷分辨通道：`d = 1`

### `bose_distribution(E, μ, T)`

计算玻色分布函数：

```math
g(E) = \frac{1}{\exp((E-\mu)/T)-1}.
```

约束：

- `T >= 0`
- `E > μ`

### `stable_meson_number_density(mass, T; μ=0.0, degeneracy=1, qmax=nothing, num_q_nodes=256)`

计算稳定粒子极限介子数密度：

```math
n_M = d_M \int_0^\infty \frac{dq\,q^2}{2\pi^2}
\frac{1}{\exp((E_M-\mu_M)/T)-1},
\qquad E_M=\sqrt{q^2+m_M^2}.
```

### `stable_kpi_ratio(m_pi, m_K, T; μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...)`

返回稳定粒子极限的聚合通道 `K/π` 数密度比值。

### `stable_kpi_scan(temperatures; m_pi, m_K, μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...)`

对一组温度执行稳定粒子极限扫描，返回：

- `temperatures`
- `n_pi`
- `n_K`
- `kpi_ratio`

### `strict_bw_meson_number_density(mass, gamma, T; ...)`

计算当前 Stage-1 的 reduced strict BW 单通道介子数密度 helper。

当前采用：

```math
E_M(q)=\sqrt{q^2+m_M^2},
\qquad
\Gamma_M(q)=\Gamma_M,
\qquad
\omega=E_M(q)+\Delta\omega.
```

并计算：

```math
n_M^{BW,red}(T)
= d_M \int_0^\infty \frac{dq\,q^2}{2\pi^2}
\int_0^{\omega_{\max}} \frac{d\Delta\omega}{2\pi}
g(E_M(q)+\Delta\omega)
\frac{\Gamma_M/2}{\Delta\omega^2+\Gamma_M^2/4}.
```

当前返回除 `density` 外，还包含：

- `q_integral_estimate`
- `omega_shell_at_qmax`
- `mode`
- 当前积分配置回显

### `strict_bw_meson_density_summary(pi_mass, pi_gamma, k_mass, k_gamma, T; ...)`

基于 `strict_bw_meson_number_density` 聚合 `π/K` 两个通道，返回：

- `n_pi`
- `n_K`
- `kpi_ratio`
- `pi_density`
- `k_density`

### `strict_bw_qpole_meson_number_density(meson, mass0, gamma0, quark_params, thermo_params; ...)`

计算当前 Stage-2 的 `q` 依赖复极点 strict BW 单通道介子数密度 helper。

当前做法是：

1. 在 `q` 网格上逐点调用介子极点方程求解器；
2. 用上一个 `q` 点的 `(mass, gamma)` 作为下一个 `q` 点的 continuation seed；
3. 将得到的 `E_M(q)` 与 `Gamma_M(q)` 代回 strict BW 双积分核；
4. 内层积分按严格口径直接执行 `\omega \in [0,\omega_{\max}]`，而不是只积 `\omega \ge E_M(q)` 的右半边。

当前返回除 `density` 外，还包含：

- `q_values`
- `E_values`
- `gamma_values`
- `residual_norms`
- `converged_flags`
- `accepted_flags`

### `strict_bw_qpole_density_summary(pi_mass, pi_gamma, k_mass, k_gamma, quark_params, thermo_params; ...)`

基于 `strict_bw_qpole_meson_number_density` 聚合 `π/K` 两个通道，返回：

- `n_pi`
- `n_K`
- `kpi_ratio`
- `pi_density`
- `k_density`

### `phase_shift_meson_number_density(meson, quark_params, thermo_params; scheme=:current, ...)`

当前 Phase E3 最小口径下的单通道相移介子数密度 helper。

当前约束：

- 仅支持 `xi = 0`
- 支持 `:pi` / `:K` 聚合通道以及 `:pi_plus` / `:pi_minus` / `:K_plus` / `:K_minus` 电荷分辨通道
- 积分方案固定为 GL + 硬截断

当前默认参数：

- `qmax = 12`
- `q_nodes = 48`
- `omega_min = 0.05`
- `omega_max = 10`
- `omega_nodes = 48`
- `eta = 1e-6`
- `real_axis_mode = :finite_eta`
- `phase_convention = :arg_propagator`
- `phase_display = :unwrapped`
- `density_policy = :strict_normal_domain`
- `noanom_policy = :none`

当前支持两个正式 `scheme`：

- `:current` -> `:phase_shift_current`
  - 当前生产默认口径
  - 使用 `F(\delta)=\delta`
- `:gbu_reference`（兼容别名 `:gbu`, `:generalized_bu`）-> `:phase_shift_gbu_reference`
  - 当前更严格参考口径
  - 使用 `F(\delta)=\delta-\frac{1}{2}\sin 2\delta`

real-axis 分支：

- `real_axis_mode=:finite_eta`
  - legacy / 默认路径
  - 通过有限虚部展宽调用 finite-width polarization backend
  - 要求 `eta > 0`，不允许把 `eta=0` 静默塞进 finite-width backend
- `real_axis_mode=:pv_b0_eta0`（兼容别名 `:bu2020_pv_eta0`, `:eta0_pv_b0`）
  - BU2020 FIG2 审计所需的独立实轴 principal-value 分支
  - 返回 metadata 中 `eta=0.0`
  - 通过 real-axis polarization 路径进入 `OneLoopIntegrals.B0` 的 PV 实部与解析虚部处理

相位口径：

- `phase_convention=:arg_propagator`：legacy 默认，对传播子取相位
- `phase_convention=:arg_inverse_propagator`：BU2020 诊断口径，对 inverse propagator 取相位并按当前符号约定返回
- `phase_display=:unwrapped`：默认密度核使用 unwrap 后的相移，不强制限制到 `0..pi`
- `phase_display=:fold_0_pi`：显式诊断/FIG3-like 口径，先用 `pi - abs(mod(delta, 2pi) - pi)` 映射到 `0..pi` 再进入密度权重

Bose-domain policy：

- `density_policy=:strict_normal_domain`
  - 默认生产/审计口径
  - 当积分支持触及 `omega <= μ_M` 时不继续积分，返回 `density=NaN` 与 `status=:unsafe_bose_domain`
- `density_policy=:excitation_only_E_gt_mu`
  - 显式诊断延拓，只在 `omega > μ_M` 区域积分
- `density_policy=:x_min_cut`
  - 显式诊断延拓，按 `(omega-μ_M)/T >= bose_x_min` 设定下界

No-anomalous policy：

- `noanom_policy=:none`
  - 默认 full phase-shift 输出，不做 anomalous-mode 扣除
- `noanom_policy=:low_energy_branch_subtraction`
  - reconstructed diagnostic policy
  - 按 temp7 审计口径，在 `K_plus` 的 folded `0..pi` display phase 上识别 positive-energy Landau window 内的低能正相位连通分量，并把该低能分量置零
  - 该口径用于显式输出 no-anomalous 对照，不改变默认生产口径，也不改变上游 FixedMu 分支选择策略

治理约束：

- `current` 保持默认正式生产主线
- `gbu_reference` 作为可重复运行的 stricter reference / analysis branch
- 二者共用同一套 `workflow + continuation + scan` 契约，不允许脚本层平行重组流程
- `pv_b0_eta0` 与 `finite_eta` 必须在 API 和返回 metadata 中保持可区分；`B0` PV 奇点处理不能和 Bose occupation pole 混为一谈

返回值除 `density` 外，还包含：

- `q_integral_estimate`
- `omega_shell_at_qmax`
- `scheme`
- `real_axis_mode`
- `eta`
- `polarization_backend`
- `phase_convention`
- `phase_display`
- `density_policy`
- `noanom_policy`
- `noanom_applied`
- `noanom_removed_component_count`
- `noanom_removed_omega_min`
- `noanom_removed_omega_max`
- `noanom_landau_omega_min`
- `noanom_landau_omega_max`
- `unsafe_bose_count`
- `min_E_minus_mu`
- `bose_x_min`
- `status`
- `message`
- 当前积分配置回显

### `phase_shift_meson_density_summary(quark_params, thermo_params; scheme=:current, ...)`

基于 `phase_shift_meson_number_density` 聚合 `π/K` 两个通道，返回：

- `n_pi`
- `n_K`
- `kpi_ratio`
- `pi_density`
- `k_density`
- `scheme`
- `real_axis_mode`
- `phase_convention`
- `density_policy`
- `noanom_policy`
- `noanom_removed_component_count`
- `unsafe_bose_count`
- `min_E_minus_mu`
- `status`
- `message`

当前它是 workflow 层 Phase-E3 后处理入口所依赖的正式数值 helper。
