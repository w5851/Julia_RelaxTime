# T150 回退链路修复与 T190/T200 毛刺分析设计

## 1. 背景与目标

本设计面向同一条 `plan_b` 产数链路的两类问题：

1. `T=150, muB=0` 在小 `xi` 前段（重点是 `xi=-0.5`）出现平衡解未收敛并被跳过，导致缺点。
2. `T=190/200` 数据存在明显不光滑波动（毛刺），需给出可审计、可复核的解释。

本轮范围约束（按实际落地更新）：

- 允许改 `scripts`，且优先通过 `scripts` 调用 `src` 中已有 fallback 能力，而不是脚本层重造求解器。
- 允许最小化 `src` 变更：为 workflow pipeline 增补 diagnostics stage 与 manifest 扩展字段、为平均散射率补充 band 统计输出接口。
- 主分析数据口径以 `D:\Desktop\Julia_RelaxTime\.worktrees\repro-main-oldparams\data\outputs\tmp\repro_main_oldparams\results\relaxtime\plan_b\plan_b_merged.csv` 为准。

## 2. 范围与非范围

### 2.1 In Scope

- 调整 `scripts/relaxtime/run_gap_transport_scan.jl` 中平衡求解调用路径：优先走 `Models.solve(...)`（利用 src 既有 multiseed/fallback 机制）。
- 保留当前 `Models.solve_constraint(...)` 兼容兜底路径。
- 增加“跳过点失败机制”的 sidecar 落盘。
- 最小复现实验：`T=150, muB=0, xi=-0.5`（必要时扩到 `-0.48, -0.46`）。
- 基于主 merged 数据完成 `T190/T200` 毛刺链路归因，并输出置信度。

### 2.2 Out of Scope

- 不做全量重跑 `T150/190/200`。
- 不引入新的物理模型参数与新基线准入流程。

## 3. 现状链路与已知缺口

### 3.1 现状主链

- 入口：`scripts/relaxtime/run_manual_relaxation_scan_workflow.jl` 编排 `plan_b`。
- 计算：`scripts/relaxtime/run_gap_transport_scan.jl`。
- 平衡态：当前脚本 `solve_models_equilibrium` 直接调用 `Models.solve_constraint(FixedMu)`。
- 输运：`TransportWorkflow.solve_gap_and_transport`。

### 3.2 关键缺口

- `solve_equilibrium_with_diagnostics` 在无有效候选时抛错，外层只做 `@warn ... skipped`，失败点不入 CSV。
- 现有 CSV 的 quality 字段无法覆盖“被跳过点”。
- 已有 metadata 仅能看到汇总与成功点，不能解释每个失败点为何失败。

## 4. 设计方案

### 4.1 问题1：脚本调用 src 既有 fallback 能力

在 `run_gap_transport_scan.jl` 中修改 `solve_models_equilibrium`：

1. 首选 `Models.solve(model, FixedMu(), T_fm, muq_fm; ...)`。
   - 传入现有 `seed_guess`（沿用 phase-aware 连续性策略）。
   - 传入 `xi/p_num/t_num/residual_norm_max`。
   - 使用 src 既有 `auto_multiseed_fallback`（默认启用）能力。
2. 将 `Models.solve(...)` 返回映射为脚本当前 `eq` 结构：
   - `x_state`, `mu_vec`, `masses`, `iterations`, `residual_norm`, `omega`, `solver_backend`。
3. 若 `Models.solve(...)` 返回不可用或异常，再回落 `Models.solve_constraint(...)`（兼容兜底）。

设计动机：

- `src` 已有 fixed-mu multiseed/fallback 机制，脚本应复用而非重复实现。
- 保持后续 `TransportWorkflow.solve_gap_and_transport` 及 CSV schema 稳定。

### 4.2 问题1：失败机制 sidecar 落盘

新增可选失败落盘文件（append 模式，默认写到 temp/worktree 外路径）：

- 文件：`*_failed_points.csv`
- 字段：
  - `T_MeV`
  - `muB_MeV`
  - `xi`
  - `seed_source`
  - `phase_prev`
  - `phase_curr_hint`
  - `error_type`
  - `error_message`
  - `timestamp`

规则：

- 每当扫描点触发 `SCAN POINT FAILED — skipped` 时写入一行。
- 不改变主 CSV 结构，避免影响现有绘图与合并脚本。

### 4.3 问题2：T190/T200 毛刺分析策略

先做“读数+源码链路”归因，不立即重跑：

1. 从主 merged 提取 `T=190/200` 全 `xi` 段的关键派生量：
   - `eta_over_s`, `sigma_over_T`, `zeta_over_s`, `sigma_over_T_over_eta_over_s`, `zeta_over_s_over_eta_over_s`
   - `tau_*`, `tauinv_*`, `m_*`, `Phi/Phibar`, `s_fm3inv`
2. 对毛刺区做链路映射：
   - 平衡态量是否平滑。
   - 跳变是否主要来自 `tau`/`tauinv` 通道。
   - 比值型观测是否因分母小值放大（数值放大机制）。
3. 若证据不足，再在 temp 路径做定向重跑并输出更多中间量（非默认）。

## 5. 实施步骤

1. 改 `run_gap_transport_scan.jl` 平衡态求解入口（`Models.solve` 优先 + `solve_constraint` 兜底）。
2. 增加失败 sidecar 写入逻辑与参数化输出路径。
3. 最小验证运行 `T150, muB0, xi=-0.5`（必要时 `-0.48,-0.46`）。
4. 基于指定 merged 数据完成 `T190/T200` 毛刺证据链分析。
5. 输出分析文档到 `docs/analysis/relaxtime/t150-convergence-and-t190-t200-spikes-analysis.md`。

## 6. 验证与验收标准

### 6.1 问题1验收

- `xi=-0.5` 点不再被跳过（有有效输出行）。
- 若仍失败，`*_failed_points.csv` 必须完整记录失败机制。
- 修改后不破坏既有 CSV header 与下游合并流程。

### 6.2 问题2验收

- 给出“数值链路证据”与“物理链路证据”两套证据。
- 明确结论分类（数值误差主导 / 物理可解释不光滑 / 混合因素）。
- 为每个结论给出置信度等级与依据。

## 7. 风险与缓解

- 风险：`Models.solve` 与脚本期望字段存在口径差。
  - 缓解：在映射层做显式字段验证与异常兜底到 `solve_constraint`。
- 风险：毛刺归因仅凭现有字段不足。
  - 缓解：启用 temp 路径下的定向重跑与中间量扩展输出。
- 风险：并行进程污染结果目录。
  - 缓解：所有增量验证与诊断输出写入独立 temp 路径。

## 8. 输出物

- 代码修改：`scripts/relaxtime/run_gap_transport_scan.jl`（仅 scripts 层）。
- 失败 sidecar：运行产出 `*_failed_points.csv`（temp 路径）。
- 分析报告：`docs/analysis/relaxtime/t150-convergence-and-t190-t200-spikes-analysis.md`。
- 结论包含：
  - 结论（含是否为数值问题）
  - 证据链
  - 结论置信度（High/Medium/Low）
