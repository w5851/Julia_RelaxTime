# Models Meson Density Workflow

本文档从 `Models` 统一入口视角描述 meson density workflow。底层实现位于 `src/models/workflow_apps/MesonDensityWorkflow.jl`，领域细节页保留在 `docs/api/relaxtime/workflow/MesonDensityWorkflow.md`。

## 首选入口

- `Models.solve_meson_density_from_meson_point`
- `Models.solve_gap_and_meson_density_point`
- `Models.solve_strict_bw_meson_density_from_meson_point`
- `Models.solve_gap_and_strict_bw_meson_density_point`
- `Models.solve_phase_shift_meson_density_from_meson_point`
- `Models.solve_gap_and_phase_shift_meson_density_point`

## `solve_meson_density_from_meson_point`

这是后处理入口。它只消费 `Models.solve_gap_and_meson_point` 的返回值，不重复求平衡态或介子极点。
当上游平衡态来自 `Models.solve(model, FixedAsymmetricRho(...), T_fm)` 等非
FixedMu source 时，应先用 `Models.solve_meson_point_from_equilibrium` 构造同构
`meson_point`，再进入本入口。

它适合：

- 已经拿到 `meson_results`，想在同一个点上继续算 `n_pi`、`n_K`、`K/π`
- 保持“介子数密度挂在介子质量 workflow 之后”的主链结构
- 避免脚本层自行拼接 `thermo_params`、介子质量与数密度公式

## `solve_gap_and_meson_density_point`

这是完整闭环入口：先求平衡态与介子谱，再在其返回值上追加 `meson_density` 字段。

它适合：

- 调用方只有 `(T_fm, mu_fm, xi)`，想直接得到介子质量与稳定粒子极限数密度
- 希望脚本只走 `Models` 统一入口，而不是直接 include workflow 实现文件

## `solve_strict_bw_meson_density_from_meson_point`

这是 strict BW workflow 的后处理入口。它继续只消费
`Models.solve_gap_and_meson_point` 的返回值，但把稳定粒子极限替换为 reduced strict BW 双积分后处理。

它适合：

- 已经拿到 `meson_results`，想在同一点上继续算有限宽度极点近似下的 `n_pi`、`n_K`、`K/π`
- 保持 workflow 唯一入口，不让 strict BW 又退回到脚本层私有实现

当前口径：

- `E(q)=sqrt(q^2+m^2)`
- `Gamma(q)=Gamma(q=0)`
- Stage1 在有限 `omega_min..omega_max` 谱窗口上积分单位 Lorentzian 权重
- `omega_min` 默认 `0.05 fm^-1`，必须高于介子化学势以避开 Bose pole
- 内层实现使用 `theta = atan(2(omega-E(q))/Gamma)` 的等价变量变换，以保证小宽度极限连续回到 stable fallback

当前也支持：

- `stage = :stage1_reduced`
- `stage = :stage2_qpole`

其中 `stage2_qpole` 会在 `q` 网格上逐点重解复极点 `z_p(q)`，并按严格口径对
`\omega \in [0,\omega_{\max}]` 执行 BW 双积分。

## `solve_gap_and_strict_bw_meson_density_point`

这是 reduced strict BW 版本的完整闭环入口：先求平衡态与介子谱，再追加
`strict_bw_meson_density` 字段。

## `solve_phase_shift_meson_density_from_meson_point`

这是当前 Phase E3 最小 BU 相移双积分的后处理入口。它仍然只消费
`Models.solve_gap_and_meson_point` 的返回值，但把稳定粒子极限替换为相移双积分后处理。

它适合：

- 已经拿到 `meson_results`，想在同一点上继续算当前最小口径下的 `n_pi`、`n_K`、`K/π`
- 保持 workflow 唯一入口，不让分析脚本直接拼接 `quark_params` / `thermo_params` / 相移积分

当前默认口径：

- `qmax = 12`
- `q_nodes = 48`
- `omega_max = 10`
- `omega_nodes = 48`
- `eta = 1e-6`
- `real_axis_mode = :finite_eta`
- `phase_convention = :arg_propagator`
- `phase_display = :unwrapped`
- `density_policy = :strict_normal_domain`
- `noanom_policy = :none`

当前同一 workflow 入口支持两套 `scheme`：

- `scheme = :current`
  - 输出 `scheme = :phase_shift_current`
  - 作为默认正式生产主线
- `scheme = :gbu_reference`（兼容 `:gbu` / `:generalized_bu`）
  - 输出 `scheme = :phase_shift_gbu_reference`
  - 作为可重复运行的 stricter reference / analysis branch

BU2020/temp7 审计相关参数：

- `real_axis_mode=:finite_eta`
  - legacy finite-width smoothing 路径
  - 要求 `eta > 0`
- `real_axis_mode=:pv_b0_eta0`
  - 独立的 `eta=0` real-axis principal-value 分支
  - 返回 `eta=0.0` 与 `polarization_backend=:pv_b0_real_axis`
- `phase_convention=:arg_inverse_propagator`
  - 用于 BU2020 FIG2 审计的 inverse-propagator phase 口径
- `phase_display=:fold_0_pi`
  - 显式 FIG3-like / temp7 display 诊断口径；把相移映射到 `0..pi` 后再用于密度权重
- `density_policy=:strict_normal_domain`
  - 默认遇到 `omega <= μ_M` 支持时返回 `status=:unsafe_bose_domain`，不静默 clamp/skip
- `density_policy=:excitation_only_E_gt_mu` / `:x_min_cut`
  - 仅作为显式诊断延拓，不应解释为文献明示事实
  - 在 combined scan 中只作用于 `phase_shift_current` 与 `phase_shift_gbu_reference`；`stable` 与 `strict_bw_stage*` 继续按 strict Bose-domain guard 输出 unsafe/NaN 行
- `noanom_policy=:low_energy_branch_subtraction`
  - 按 temp7 审计通过的 reconstructed diagnostic 口径删除 `K_plus` 低能 anomalous 分支
  - 只影响相移 density kernel 的 no-anomalous 对照，不改变 FixedMu 默认分支选择

## `solve_gap_and_phase_shift_meson_density_point`

这是相移双积分版本的完整闭环入口：先求平衡态与介子谱，再追加
`phase_shift_meson_density` 字段。

## 返回结构

`solve_gap_and_meson_density_point` 返回值保留 `solve_gap_and_meson_point` 的核心字段：

- `equilibrium`
- `quark_params`, `thermo_params`
- `meson_results`

并新增：

- `meson_density`
  - `m_pi`, `m_K`
  - `n_pi`, `n_K`
  - `kpi_ratio`
  - 简并因子、化学势与积分节点配置
- `strict_bw_meson_density`
  - `m_pi`, `m_K`
  - `gamma_pi`, `gamma_K`
  - `n_pi`, `n_K`
  - `kpi_ratio`
  - `qmax`, `q_nodes`, `omega_min`, `omega_max`, `omega_nodes`
- `phase_shift_meson_density`
  - `m_pi`, `m_K`
  - `n_pi`, `n_K`
  - `kpi_ratio`
  - `scheme`
  - `pi_density`, `k_density`
  - `qmax`, `q_nodes`, `omega_min`, `omega_max`, `omega_nodes`, `eta`
  - `real_axis_mode`, `polarization_backend`, `phase_convention`
  - `density_policy`, `unsafe_bose_count`, `min_E_minus_mu`, `bose_x_min`
  - `noanom_policy`, `noanom_removed_component_count`, `noanom_landau_omega_min`, `noanom_landau_omega_max`
  - `status`, `message`

## 当前边界

当前 workflow 已覆盖：

1. 稳定粒子极限
2. strict BW Stage1 reduced
3. strict BW Stage2 q-pole
3. Phase E3 最小相移双积分

其中 strict BW 与 phase-shift 入口当前仍处于严格受限状态：

- reduced strict BW：
  - 只消费 workflow 当前点给出的 `mass/gamma`
  - 对应 `stage = :stage1_reduced`
- q-pole strict BW：
  - 在 `q` 网格上逐点调用介子极点方程
  - 内层 `\omega` 积分按完整 ` [0,\omega_{\max}] ` 口径执行
  - 当前依赖 continuation seed 串行续算
  - 对应 `stage = :stage2_qpole`
- Phase E3 phase-shift：
  - 仅支持 `xi = 0`
  - 支持 `π/K` 聚合通道以及 `pi_plus/pi_minus/K_plus/K_minus` 电荷分辨通道
  - 积分方案固定为 GL + 硬截断
  - `current` 为默认生产口径，`gbu_reference` 为并列参考分支
  - `finite_eta` 为 legacy 默认 real-axis 分支，`pv_b0_eta0` 为 BU2020 审计所需的独立 PV 分支
  - strict Bose-domain policy 默认标记 unsafe；诊断延拓必须显式传入
  - no-anomalous 只作为显式 reconstructed diagnostic policy；默认 full phase shift 不扣除 anomalous 分支

脚本层组合入口：

- `scripts/relaxtime/run_combined_meson_density_scan.jl`
  - Bridge-style 组合：`scan path` × `density regime`
- 当前实现 `--path tmu` 与 smoke-only `--path trho_asymmetric`
- `--path tmu` 默认四口径为 `stable,strict_bw_stage1,phase_shift_current,phase_shift_gbu_reference`
- 支持 `--muq-values` 或 `--mumin/--mumax/--mustep` 生成多个固定 `mu_q` 的 T 扫描；多 `mu_q` 输出会生成 FIG3-like heatmap SVG
- `--path trho_asymmetric` 使用 `FixedAsymmetricRho(rho_target, asym_ud_ratio_target, asym_s_target)` 作为 density-constrained equilibrium source，再通过 `Models.solve_meson_point_from_equilibrium` 后处理；当前只作为 smoke / diagnostic path，不作为正式高精度生产入口
- `trho_asymmetric` 默认按温度分组并在每个温度内做 rho 连续扫描，`--trho-reverse-rho=true` 时从高 rho 到低 rho 传递 equilibrium seed；`--no-trho-reverse-rho` 仅用于顺序敏感诊断
- `--density-policy x_min_cut --bose-x-min <x>` 可把 Bose x 下界传给 BU/GBU phase-shift regimes；该选项不延拓 stable/BW 口径
- `trho_asymmetric` 额外输出 `constraint_mode`、`rho_target`、`rho_norm`、`rho_u_fm3`、`rho_d_fm3`、`rho_s_fm3`、`rho_u_over_rho_d`、`asym_ud_ratio_target`、`asym_s_target`、`constraint_residual_norm`、`mu_u_MeV`、`mu_d_MeV`、`mu_s_MeV`、`muB_MeV`、`muQ_MeV`、`muS_MeV`
- 正式数据默认写入 `data/outputs/results/...`；图像和 `plot_manifest.json` 默认写入对应 `data/outputs/figures/...`，也可通过 `--figure-dir` 覆盖
- `scripts/analysis/relaxtime/render_combined_meson_density_temperature_scan.py` 可从单 `mu_q` 统一 CSV 渲染温度扫描高 DPI PNG
- `scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py` 可从多 `mu_q` 统一 CSV 渲染 FIG3-like 高 DPI PNG
  - 输出 CSV、README、SVG 与图像 manifest，适合把同一批状态点的多口径介子数密度结果放在同一份可审计产物中

后续 full strict BW 与更完整的 BU 扩展仍应沿同一 workflow 链继续后接，而不是回到脚本层重组流程。
