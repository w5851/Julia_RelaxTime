# 介子热力学 regression 与 validation 治理方案

更新日期：2026-07-10

当前状态：backlog 治理路线。point-level、短路径和 plot-review regression 均已落地；当前剩余问题是 external validation 是否具备升格条件。本页不属于当前 active 执行批次。

重分类说明：

- 本页包含跨批次的长期 validation 升格门槛，按 `docs/dev/README.md` 规则从 active 移入 backlog；
- 已落地事实继续保留作为下一次任务拉起时的输入；
- 后续执行前必须重新核对文献映射、当前 baseline 和测试入口，不直接照搬 2026-05-09 的任务顺序。

前序任务说明：

- `2026-05-07_介子热力学扩展任务单.md` 已在第一阶段工程落地完成后归档；
- 文献忠实复现问题不在本页内收口，而由独立审计分支继续处理。

---

## 1. 目标

本方案只回答三件事：

1. meson thermo 主线下一步该如何做 point-level regression；
2. canonical `mu_B = 0` 结果资产何时值得升格为 path-level regression；
3. 哪些外部对照可以进入 `tests/validation/`，哪些只能停留在 plot-review 或研究说明层。

本方案不做：

1. 不直接导出新 baseline；
2. 不直接更新 `tests/baselines/`；
3. 不把当前 `Maslov & Blaschke 2023` 方向直接包装成 formal validation gate。

---

## 2. 当前资产

当前 meson thermo 已有：

1. workflow：
   - `Models.solve_gap_and_meson_thermo_point`
   - `Models.solve_gap_and_strict_bw_meson_thermo_point`
   - `Models.solve_gap_and_phase_shift_meson_thermo_point`
2. 合同：
   - `Models.build_meson_thermo_contract_row`
3. canonical case：
   - `scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl`
   - 默认目录 `data/outputs/results/relaxtime/meson_thermo/canonical_muB0_phase_shift_current_pi_sigma_pi/`
4. 最小自动验证：
   - `tests/unit/relaxtime/test_meson_thermodynamics.jl`
   - `tests/unit/models/test_meson_thermo_workflow.jl`
   - `tests/integration/relaxtime/test_phase_shift_meson_thermo_scan_smoke.jl`

当前尚缺：

1. external validation gate。

---

## 3. 分层原则

将 meson thermo 治理拆成四层：

1. `unit`
   - 锁公式局部关系；
   - 例如 `P_meson = P_qp + P_ld`、`ld_cutoff` 影响 `LD` 分区、`thermo_derivation_mode = omega_total_ad`。
2. `integration`
   - 锁 workflow 与结果资产合同；
   - 例如 canonical `scan.csv / README.md / effective_config.json / run_manifest.json`。
3. `regression`
   - 锁项目自产固定点或固定路径不漂移；
   - 不要求 external truth。
4. `validation`
   - 锁外部 reference 或文献 target；
   - 必须先证明路径、口径、模型差异都足够可解释。

当前 meson thermo 已完成前 3 层中的第一版 point-level regression；后续继续补 path-level regression，再决定第 4 层。

---

## 4. point-level regression 方案

### 4.1 目标

point-level regression 的目的不是追求“最终物理真值”，而是先锁：

1. 主 workflow 入口没有漂移；
2. `phase_shift_current` 主口径没有意外变号、爆炸或回退到 legacy 差分；
3. `pi/sigma_pi` 主双通道的量级关系没有被无意改写。

### 4.2 首版 baseline 文件

当前已落地：

- `tests/baselines/relaxtime/baseline_meson_thermo_fixedpoints_v1.csv`

当前已落地 regression 文件：

- `tests/regression/relaxtime/test_meson_thermo_fixedpoint_regression.jl`

### 4.3 建议固定点集合

当前首版固定点只锁最小代表点，不追求大覆盖：

1. canonical `mu_B = 0`, `xi = 0`
2. 温点建议：
   - `T = 170 MeV`
   - `T = 210 MeV`
   - `T = 230 MeV`

原因：

1. `170 MeV` 可覆盖较低温热介子仍有可见贡献的窗口；
2. `210 MeV` 已是当前 phase-shift smoke 常用点；
3. `230 MeV` 更接近高温端，便于暴露 `QP/LD` 分区漂移。

### 4.4 建议锁定字段

point-level baseline 先只锁 phase-shift 主口径最必要字段：

1. 点位 / 入口：
   - `T_MeV`
   - `muB_MeV`
   - `xi`
   - `workflow`
   - `channel_set`
   - `phase_shift_variant`
   - `thermo_derivation_mode`
2. pressure / EOS：
   - `P_meson`
   - `P_meson_qp`
   - `P_meson_ld`
   - `P_quark_meanfield`
   - `P_total`
   - `entropy`
   - `epsilon`
   - `trace_anomaly`
3. 双通道解释字段：
   - `P_primary`
   - `P_secondary`
   - `P_primary_qp`
   - `P_primary_ld`
   - `P_secondary_qp`
   - `P_secondary_ld`
4. 治理字段：
   - `ld_cutoff`
   - `ld_cutoff_mode`
   - `ld_threshold_mode`

当前不建议首版 point-level baseline 直接锁：

1. `P_K_*` 兼容槽位；
2. `phase_structure`；
3. 低价值的诊断性元字段。

### 4.5 建议断言方式

第一版建议分两类断言：

1. 枚举断言：
   - `workflow == phase_shift_current`
   - `channel_set == pi,sigma_pi`
   - `thermo_derivation_mode == omega_total_ad`
   - `ld_cutoff_mode == match_qmax`
   - `ld_threshold_mode == omega_lt_q`
2. 数值断言：
   - 核心 thermodynamic 标量用 `isapprox`
   - 首版建议 `rtol = 5e-3`
   - `atol` 按字段量级设定，优先 `1e-8`

如果后续发现 `entropy / epsilon / trace_anomaly` 比纯 pressure 更敏感，可把：

1. pressure 字段继续保留 `5e-3`；
2. 派生量单独放宽到 `1e-2`；

但这必须伴随误差来源说明，不允许为了过测试临时放宽。

---

## 5. path-level regression 方案

### 5.1 目标

path-level regression 用来锁：

1. canonical 脚本没有改坏结果资产结构；
2. 温度路径上没有明显跳点、缺点、字段漂移；
3. `README` 摘要与 `run_manifest` 语义没有退化。

### 5.2 path-level baseline 与测试文件

当前已落地：

- baseline：
  - `tests/baselines/relaxtime/baseline_meson_thermo_canonical_muB0_path_v1.csv`
- regression：
  - `tests/regression/relaxtime/test_meson_thermo_canonical_muB0_path_regression.jl`

### 5.3 路径范围建议

不建议首版就锁整条重路径的高密度扫描。

建议首版 path-level regression 只锁：

1. `mu_B = 0`
2. `xi = 0`
3. `phase_shift_current`
4. `pi/sigma_pi`
5. `T = 170:20:230 MeV`

原因：

1. 与 point-level 温点复用，便于比较；
2. 点数少，运行时长可控；
3. 已足够覆盖低温、中温、高温三个区段。

### 5.4 路径层锁定对象

path-level regression 不需要复刻全部 CSV 列，只建议锁：

1. 行数与温点序；
2. `equilibrium_converged`；
3. `P_meson`
4. `P_meson_qp`
5. `P_meson_ld`
6. `P_total`
7. `trace_anomaly`
8. `thermo_derivation_mode`

并额外保留两个结构性断言：

1. 每行都满足 `P_meson ≈ P_meson_qp + P_meson_ld`
2. `thermo_derivation_mode` 不得退回 `workflow_fd_legacy`

### 5.5 与现有 integration smoke 的关系

当前的：

- `tests/integration/relaxtime/test_phase_shift_meson_thermo_scan_smoke.jl`

只负责：

1. 脚本能跑；
2. 资产存在；
3. 合同字段没丢。

后续若补 path-level regression，它负责：

1. 锁数值路径；
2. 锁摘要字段；
3. 锁 canonical 结果目录的最小工程含义。

两者不能互相替代。

---

## 6. plot-review 层建议

当前 canonical `mu_B = 0` case 已有首个 plot-review 资产，并已进入独立 regression。

当前目录：

- `data/outputs/results/relaxtime/meson_thermo/plot_review/<case_slug>/`

当前首个 case：

- `canonical_muB0_phase_shift_current_pi_sigma_pi`

当前最小产物：

1. `workflow_scan.csv`
2. `plot_review_summary.csv`
3. `pressure_overlay.png`
4. `trace_anomaly_overlay.png`
5. `qp_ld_split.png`
6. `README.md`

当前已新增：

- baseline：
  - `tests/baselines/relaxtime/baseline_meson_thermo_plot_review_case_v1.csv`
- regression：
  - `tests/regression/relaxtime/test_meson_thermo_plot_review_case_regression.jl`

---

## 7. validation 升格门槛

当前 `Maslov & Blaschke 2023` 方向可以作为 mechanism-aligned reference，但还不应自动升格进 `tests/validation/`。

只有当以下条件同时满足时，才建议新增 meson thermo validation：

1. **通道口径明确**
   - 项目输出与文献通道集合能一一对应；
   - 至少要能解释 `pi/sigma_pi` 与文献 `π/σ` 的对应关系。
2. **路径定义明确**
   - 能明确是 `mu_B = 0` 温扫；
   - 温窗、步长、参考归一化方式可参数化重建。
3. **模型差异可解释**
   - 当前项目与文献在 flavor、相互作用项或截断上的差异已写清；
   - 不会把“不同模型导致的系统差异”误判为回归失败。
4. **观测量稳定**
   - 至少有一组外部 target 可以稳定比较：
     - `P_meson(T)`
     - `P_total(T)`
     - `trace_anomaly(T)`
     - `QP/LD` 分拆趋势中的一个或多个。
5. **结果资产齐备**
   - 已有固定输入路径；
   - 已有 plot-review 资产；
   - 已有人为审阅说明为什么它值得升格为 validation gate。

若任一条件不满足，当前方向应继续停留在：

1. regression；
2. plot-review；
3. 文献机制对照说明。

---

## 8. 下次拉起任务时的推荐执行顺序

1. 重新审计 external reference 的通道、路径、归一化和模型差异；
2. 确认至少一组观测量具备可参数化重建、稳定 target 和明确容差来源；
3. 先生成独立 plot-review / evidence package，再决定是否进入 `tests/validation/`；
4. 若证据不足，保持 mechanism-aligned reference，不创建形式化 validation gate；
5. 只有代码语义或正式默认值确实变化时，才进入 baseline 更新流程。

当前不建议：

1. 因为已有 regression 就自动宣称 external validation 已成立；
2. 为匹配文献图而临时放宽数值容差；
3. 在缺少差异审计时批量重生 baseline。

---

## 9. 本页结论

当前 meson thermo 主线的治理结论是：

1. `phase_shift_current + pi/sigma_pi + mu_B=0` 的 point-level、短路径和 plot-review regression 已建立；
2. 当前剩余治理缺口是 external validation 证据，不是继续重复建设内部 regression；
3. external literature 在完成通道、路径和模型差异审计前，只能作为 mechanism-aligned reference；
4. baseline 更新必须由实际语义变化和数值差异证据触发，不能作为 validation 工作的默认步骤。
