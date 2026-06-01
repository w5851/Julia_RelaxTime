# MesonDensityWorkflow

将 meson workflow 的返回值后处理为介子数密度。

本页是 meson density workflow 的领域细节页，重点说明它如何与现有 meson workflow 对接，以及为什么它被约束为后处理层。

如果你只是想判断“应该从哪个 `Models` 入口开始”，优先阅读 `docs/api/models/workflows/MesonDensityWorkflow.md`。

## 入口

### `solve_meson_density_from_meson_point`

```julia
solve_meson_density_from_meson_point(
    meson_point;
    pi_channel=:pi,
    k_channel=:K,
    μ_pi=0.0,
    μ_K=0.0,
    d_pi=3,
    d_K=4,
    qmax_pi=nothing,
    qmax_K=nothing,
    num_q_nodes=256,
)
```

它要求输入对象至少包含：

- `thermo_params.T`
- `thermo_params.ξ`
- `meson_results[:pi].mass`
- `meson_results[:K].mass`

workflow 会直接从这些字段读取温度与介子质量，然后调用 `Main.MesonDensity.stable_meson_number_density` 完成数密度计算。

### `solve_gap_and_meson_density_point`

```julia
solve_gap_and_meson_density_point(T_fm, mu_fm; density_kwargs=(;), kwargs...)
```

它内部先调用 `solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)`，再执行：

```julia
solve_meson_density_from_meson_point(meson_point; density_kwargs...)
```

因此它的物理主链固定为：

```text
gap / equilibrium
  -> meson mass / threshold / gap
  -> meson density
```

而不是：

```text
script 手工取温度与质量
  -> 自行拼数密度
```

### `solve_strict_bw_meson_density_from_meson_point`

```julia
solve_strict_bw_meson_density_from_meson_point(
    meson_point;
    qmax=12.0,
    q_nodes=48,
    omega_max=10.0,
    omega_nodes=48,
    gamma_zero_tol=1e-12,
)
```

它继续消费现有 meson workflow 返回值中的：

- `thermo_params`
- `meson_results[:pi].mass`
- `meson_results[:pi].gamma`
- `meson_results[:K].mass`
- `meson_results[:K].gamma`

然后把 reduced strict BW 数值内核委托给
`Main.MesonDensity.strict_bw_meson_density_summary`。

当前它是 **Stage-1 reduced strict BW** 的正式 workflow helper，物理上保持以下约束：

- 只消费 workflow 当前点给出的 `q=0` 质量与宽度
- 采用 `E(q)=sqrt(q^2+m^2)`
- 采用 `Gamma(q)=Gamma(q=0)`
- 尚未进入 `q` 依赖复极点求解

当前同一入口也支持：

- `stage = :stage1_reduced`
- `stage = :stage2_qpole`

当选择 `stage = :stage2_qpole` 时：

1. `q` 网格上的每个节点都会重解一次介子复极点；
2. 前一个 `q` 点的 `(mass, gamma)` 会被用作下一个 `q` 点的续算 seed；
3. 内层 BW 积分会按严格 `\omega \in [0,\omega_{\max}]` 口径执行；
4. 输出中会带上 `q_values / E_values / gamma_values / residual_norms` 等 pole 诊断量。

### `solve_gap_and_strict_bw_meson_density_point`

```julia
solve_gap_and_strict_bw_meson_density_point(T_fm, mu_fm; density_kwargs=(;), kwargs...)
```

它内部先调用：

```julia
solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
```

再执行：

```julia
solve_strict_bw_meson_density_from_meson_point(meson_point; density_kwargs...)
```

因此它继续保持同一个 workflow 契约：

```text
gap / equilibrium
  -> meson mass / threshold / width
  -> strict BW meson density
```

### `solve_phase_shift_meson_density_from_meson_point`

```julia
solve_phase_shift_meson_density_from_meson_point(
    meson_point;
    scheme=:current,
    qmax=12.0,
    q_nodes=48,
    omega_min=0.05,
    omega_max=10.0,
    omega_nodes=48,
    eta=1e-6,
    real_axis_mode=:finite_eta,
    phase_convention=:arg_propagator,
    density_policy=:strict_normal_domain,
    bose_x_min=0.0,
    noanom_policy=:none,
)
```

它消费现有 meson workflow 返回值中的：

- `quark_params`
- `thermo_params`
- `meson_results`

然后把数值内核委托给 `Main.MesonDensity.phase_shift_meson_density_summary`。

当前它是 **Phase E3 最小 BU 相移双积分** 的正式 workflow helper，物理上仍保持以下约束：

- 仅支持 `xi = 0`
- 支持 `π/K` 聚合通道以及 `pi_plus/pi_minus/K_plus/K_minus` 电荷分辨通道
- 积分方案固定为 GL + 硬截断

当前 `scheme` 治理口径：

- `:current`
  - `F(\delta)=\delta`
  - 默认正式生产主线
- `:gbu_reference`（兼容 `:gbu` / `:generalized_bu`）
  - `F(\delta)=\delta-\frac{1}{2}\sin 2\delta`
  - 可重复运行的 stricter reference 输出链

当前 real-axis / Bose-domain 治理口径：

- `real_axis_mode=:finite_eta`
  - legacy 默认路径，使用有限虚部展宽；`eta` 必须大于 0
- `real_axis_mode=:pv_b0_eta0`
  - BU2020/temp7 FIG2 审计所需的独立 `eta=0 + PV B0` 实轴分支
  - 返回 metadata 中 `eta=0.0` 与 `polarization_backend=:pv_b0_real_axis`
- `phase_convention=:arg_propagator`
  - legacy 默认相位口径
- `phase_convention=:arg_inverse_propagator`
  - BU2020 诊断用 inverse-propagator phase 口径
- `density_policy=:strict_normal_domain`
  - 默认遇到 `omega <= μ_M` 支持时返回 `status=:unsafe_bose_domain` 与 `density=NaN`
- `density_policy=:excitation_only_E_gt_mu` / `:x_min_cut`
  - 只作为显式诊断延拓，不能视作文献明示公式
- `noanom_policy=:low_energy_branch_subtraction`
  - 按 temp7 reconstructed diagnostic 口径删除 `K_plus` low-energy anomalous branch
  - 不改变上游 FixedMu 默认分支选择，也不作为 full phase-shift 默认值

### `solve_gap_and_phase_shift_meson_density_point`

```julia
solve_gap_and_phase_shift_meson_density_point(T_fm, mu_fm; density_kwargs=(;), kwargs...)
```

它内部先调用：

```julia
solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
```

再执行：

```julia
solve_phase_shift_meson_density_from_meson_point(meson_point; density_kwargs...)
```

因此它继续保持同一个 workflow 契约：

```text
gap / equilibrium
  -> meson mass / threshold / gap
  -> phase-shift meson density
```

## 当前最小扫描入口

当前已提供脚本级最小扫描入口：

```text
scripts/relaxtime/run_strict_bw_meson_density_scan.jl
scripts/relaxtime/run_phase_shift_meson_density_scan.jl
```

其中 strict BW 扫描固定通过：

```text
Models.solve_gap_and_strict_bw_meson_density_point
```

并输出：

- `n_pi`
- `n_K`
- `kpi_ratio`
- `gamma_pi`, `gamma_K`
- `qmax`, `q_nodes`
- `omega_max`, `omega_nodes`
- `pi/K` 两个通道的 `q_integral_estimate`
- `pi/K` 两个通道的 `omega_shell_at_qmax`
- `pi/K` 两个通道的 `mode`

phase-shift 扫描固定通过：

```text
Models.solve_gap_and_phase_shift_meson_density_point
```

进入 workflow，并把下列结果输出到 scan CSV：

- `n_pi`
- `n_K`
- `kpi_ratio`
- `scheme`
- `qmax`, `q_nodes`
- `omega_min`, `omega_max`, `omega_nodes`
- `real_axis_mode`, `eta`, `phase_convention`
- `density_policy`, `unsafe_bose_count`, `min_E_minus_mu`, `bose_x_min`, `status`
- `noanom_policy`, `noanom_removed_component_count`, `noanom_landau_omega_min`, `noanom_landau_omega_max`
- `pi/K` 两个通道的 `q_integral_estimate`
- `pi/K` 两个通道的 `omega_shell_at_qmax`

BU2020/temp7 主线审计脚本：

```text
scripts/relaxtime/run_bu2020_meson_density_audit_scan.jl
```

该脚本不复制 temp7 代码，只通过 `Models` workflow 输出一个可审计 CSV/README，覆盖 stable、strict BW Stage1、`phase_shift_current`、`phase_shift_gbu_reference`，并显式记录 charged-channel 化学势、`pv_b0_eta0`、inverse-propagator phase、Bose-domain policy 和 `low_energy_branch_subtraction` no-anomalous 诊断状态。

## 当前设计原则

1. workflow 是唯一计算入口；
2. 脚本只能消费 workflow 返回值；
3. 数密度层不重写 meson 求解链；
4. 稳定粒子与相移双积分都应作为 meson workflow 的后处理层演进；
5. Stage1 reduced strict BW、Stage2 q-pole strict BW、以及未来更完整的 BU 扩展都应在此后处理边界内继续演进。
