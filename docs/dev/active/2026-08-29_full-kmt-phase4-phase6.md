# 完整 KMT 路线 Phase 4-6：BQS、BU A/B 与冻结线候选审查

更新日期：2026-08-29

当前状态：in progress。本文承接 `2026-08-29_full-kmt-interaction-phase0-phase1.md`，记录完整 charged KMT 耦合接回现有 BQS quark-only 平衡后的诊断实现。本文不授权修改正式 production baseline，也不把 `x_min_cut`/低节点结果当作实验拟合。

## 1. 路线与边界

本阶段采用以下候选路线：

```text
FixedMuBConservedCharges (quark-only BQS)
    -> same mean-field state and (mu_Q, mu_S)
    -> legacy K123/K4567 BU density
    -> full KMT K12/K45 BU density
    -> K+/pi+ and K-/pi- A/B along baseline_freezeout
```

- Phase 4：复用现有 `FixedMuBConservedCharges`，从其 `x_state` 构造 `FullKMTInteraction`。这里的“接回”是把完整核接到 BQS 平衡背景，不是用介子压力修改 `Omega` 或另解一套非对角夸克自能。
- Phase 5：在同一 `x_state`、同一组 flavor chemical potentials、同一积分网格上，比较旧 effective coupling 与 full charged KMT coupling。`pi^± -> K12`，`K^± -> K45`；当前 `K67` 和中性矩阵只作为输出诊断。
- Phase 6：沿 `default` profile 的 `baseline_freezeout` 稀疏扫描。扫描方向从高 `sqrt(s_NN)` 到低能 continuation，最终 CSV 按能量升序写出。

明确不做：strict-support、节点收敛、完整 `Omega_M`/Sigma_M 反馈、显式 `mu_I` 路线、`K^0` 密度生产、修改 transport 或正式 `MesonDensity` 默认语义。

## 2. 实现

- `src/relaxtime/MesonDensity.jl` 新增 `PhaseShiftInteractionSpec`，默认参数完全保持旧 `2K/(1-4KPi)`；可传入 `FullKMTInteraction` 做带电耦合注入。
- `scripts/analysis/relaxtime/meson_conserved_charge_feedback_runtime.jl` 新增 `solve_quark_only_bu_ab_point`，只求一次 BQS baseline，然后运行 legacy/full 两组四通道 BU 密度。
- `scripts/analysis/relaxtime/scan_full_kmt_bu_freezeout.jl` 是独立 diagnostic 入口，不改写 `Models.run_freezeout_meson_density_scan`。
- A/B 的完整 KMT结果只表示“给定 quark-only 背景上的 full charged coupling 后处理”。Phase 4 不应被表述成“完整 KMT 已进入 PNJL stationarity”。

## 3. 数值配置与门禁

默认沿用当前低成本设置：`p_num=8`、`t_num=4`、`qmax=4 fm^-1`、`q_nodes=4`、`omega=[0.05,3] fm^-1`、`omega_nodes=8`、`density_policy=:x_min_cut`、`bose_x_min=0.05`。BQS gap residual 上限为 `1e-5`；扫描同时保留每个通道的 Bose status、`min_E_minus_mu`、BQS/总守恒荷残差和分项耗时。

只有在以下条件都清楚后才可讨论 production 晋升：

1. full/legacy A/B 的 coupling、归一化和极化输入逐项可追溯；
2. charged BU 的支撑域和节点收敛通过独立门禁；
3. 冻结线所有目标点无静默 fallback，失败点有保留的诊断记录；
4. 与外部文献固定点或实验趋势的比较只作为独立证据，不以调参替代模型闭合。

因此 Phase 6 的默认 `production_candidate_status` 为 `not_authorized`。

## 4. 验证记录

已完成以下验证：

- `julia --project=. tests/unit/relaxtime/test_meson_density.jl`：全部测试集合通过，包含
  `PhaseShiftInteractionSpec` 与 FullKMT 注入测试；
- `julia --project=. tests/integration/relaxtime/test_meson_conserved_charge_feedback_script_contract.jl`：69/69；
- `julia --project=. tests/unit/models/test_fixed_mub_conserved_charges.jl`：36/36；
- `julia --project=. tests/unit/relaxtime/test_meson_interaction_kernel.jl`：63/63；
- `julia --project=. tests/unit/relaxtime/test_meson_rpa.jl`：31/31；
- `julia --project=. tests/unit/relaxtime/test_meson_rpa_adapter.jl`：36/36；
- `git diff --check` 与三个分析脚本的 `Meta.parseall` 均通过；
- full-KMT 稀疏冻结线扫描 7/7 点 `quark_only_full_kmt_bu_ok`，失败行和完整数值配置字段均可序列化。

原始扫描 CSV 保持本地 diagnostic 产物，不提交到版本库。

## 5. 2026-08-29 代表点与冻结线诊断

运行入口：

```text
julia --project=. scripts/analysis/relaxtime/scan_full_kmt_bu_freezeout.jl
```

实际使用 7 点网格
`sqrt(s_NN)=[3,4.5,7.7,11.5,20,62.4,200] GeV`。这是在已有 partial-feedback 热启动中位耗时 `0.2478439 s` 上乘以 A/B 保守因子 `2` 后，通过相同 600 s 网格选择器得到的结果。扫描从 200 GeV 向低能 continuation，7/7 点 quark-only gap 与八个 BU 通道均返回 `ok`。

| `sqrt(s_NN)` [GeV] | `K+/pi+` legacy | `K+/pi+` full | `K-/pi-` legacy | `K-/pi-` full |
|---:|---:|---:|---:|---:|
| 3.0 | 10.2075 | 10.2073 | 0.02457 | 0.02457 |
| 4.5 | 7.7723 | 7.7717 | 0.05204 | 0.05204 |
| 7.7 | 6.1577 | 6.1569 | 0.06327 | 0.06326 |
| 11.5 | 5.0157 | 5.0150 | 0.06231 | 0.06230 |
| 20.0 | 3.5459 | 3.5455 | 0.13056 | 0.13054 |
| 62.4 | 1.8038 | 1.8038 | 0.43748 | 0.43747 |
| 200.0 | 1.2393 | 1.2393 | 0.69488 | 0.69488 |

Full/legacy ratio 相对变化约为 `0.00047%--0.0233%`（绝对值取最大处在 `K-/pi-` 的 11.5 GeV 点）。`K03_P`/`K38_P` 的扫描范围约为 `1.3e-7--8.3e-6 fm^2` 与 `-1.9e-7--1.2e-5 fm^2`；在当前 BQS 背景上，`phi_u-phi_d` 很小，所以交叉耦合对 charged BU ratio 的直接影响暂时低于化学势、质量、相移窗口和 Bose 截断的影响。

quark-only baseline 的 BQS 残差为 `O(1e-15)`，但若把后处理介子荷机械加入总守恒荷，`residual_norm_full` 在低能点可达 `6.68`。这不是求解失败，而是明确显示本路线没有做 meson-charge feedback；它不能被解读为 full-KMT 的 thermodynamic equilibrium。

当前定性判断：`K+/pi+` 随 `sqrt(s_NN)` 单调下降，`K-/pi-` 总体上升但低能段有轻微非单调；没有看到 horn-like 峰，full/legacy 曲线几乎重合。两类 ratio 的数量级跨越 `O(10^-2)--O(10)`，与实验常见 `O(10^-1)` 的比较只能作为口径风险提示，不能作为逐点拟合结论。

## 6. Phase 6 gate verdict

本次实现完成了“候选路线可运行”和“full charged KMT 与旧耦合的同背景 A/B”两项工程门槛，但没有完成 production 认证。charged scalar/matrix 的公式归一化已经由显式 ladder trace 与 Goldstone 条件闭合；尚缺的是 strict ordered-retarded 后端的独立外部数值固定点，以及单电荷 phase-shift 从 `domega/(2pi)` 到 `domega/pi` 的实现迁移和稳定极限回归。其他阻断项仍是：`x_min_cut` 非严格支撑域、低节点未收敛，以及没有 `Omega_M`/`Sigma_M` 热力学反馈。因此 `production_candidate_status=not_authorized`；本地 CSV 仅作 diagnostic 保留。
