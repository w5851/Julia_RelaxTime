# Issue #130：RS transport phase-reference adapter parity 与限定重验

状态：active。runtime-switch PR #260 与 solver-free parity PR #261 已合并；
paired numerical pilot workflow PR #262 已合并到 `main@27e9642d431ba7afd845f2b187f77c0fbbe3be4d`。
本任务承接 RS transport 的 consumer parity；它不删除旧 reference，不覆盖旧 transport
production，也不把 solver-free smoke 当作数值 production 通过。

## 固定输入与边界

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- candidate source run：`32354095831`
- aggregate replay：`32451053476`
- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`
- legacy fallback：`data/reference/pnjl/{boundary.csv,cep.csv,crossover_dense.csv,spinodals.csv}`
- runtime view：strict candidate 的 certified-only 行，缺失/不合格键逐键 legacy fallback
- rollback：`--phase-reference-mode legacy`
- 本机不调用 equilibrium solver；不启动新的 C0/C1/C2、M4 production 或 transport 全量扫描

## 阶段与验收

- [x] 在合并 SHA 上完成 runtime promotion/consumer smoke；candidate/fallback 计数、candidate
  manifest SHA、`solver_called=false` 和 explicit legacy rollback 均有记录。
- [x] 通过实际 `run_gap_transport_scan.jl`/phase-guided reference resolution 入口执行
  `runtime`、`diagnostic`、`legacy` 三模式 source smoke；检查 `mu_q`/`mu_B` 映射、phase kind、
  fallback 语义和 plan anchor。
- [x] 固化 portable evidence：输入哈希、raw dry-run provenance、source summary、parity table、
  claim ledger 和停止边界。
- [x] 作者审核 solver-free parity evidence；作者已授权合并 PR #261。
- [x] 已新增 versioned workflow
  `.github/workflows/relaxtime-issue130-rs-numerical-pilot-v1.yml` 与 paired collector/test。
  numerical scope 固定为 `xi=-0.5,0,0.5`、`muB=150,450,900 MeV`、`alpha_T=1.0,1.2`，
  两套 reference 使用相同 transport 配置并在同一 job 内顺序运行；名义 18 个请求点因
  现有 `direct_coexistence` 锚点的共存侧替换形成 19 个有效计划点，已在 manifest 中显式记录。
- [x] 经作者授权后触发限定 numerical RS pilot `32684074876`；candidate runtime 与显式
  legacy rollback 均完成 19/19 点，solver、非有限、重复键、failed point 和 convergence
  检查通过，实际耗时分别为 315 s / 314 s。
- [x] 对 pilot source artifact 做 solver-free quality audit。原始 in-run verdict
  `pilot_solver_or_curve_failure` 的唯一原因是两套 reference 共同的 5 个
  `tau_u_ubar_ratio_high` 标记；collector replay 将其分类为
  `pilot_pair_complete_with_common_quality_warnings_diagnostic_only`。该分类不放宽
  `scan_quality` 阈值，也不把质量警告升级为正式 RS parity。
- [ ] 作者审核限定 numerical evidence 与共同质量警告；在审核前保持 diagnostic-only、
  fallback、rollback 和旧 production。
- [ ] phase-guided/Paper P1 等消费者完成不再直接依赖旧文件的引用审计后，才创建独立
  `old-reference retirement` PR。

## 当前证据

见：
`docs/analysis/pnjl/phase_reference/issue130_rs_transport_runtime_parity_v1/`。

- runtime view：`certified_candidate_with_legacy_fallback`
- candidate strict：boundary `7162`（certified `3091`）、crossover `1343`、CEP `93`（certified `90`）、
  spinodal `6886`
- fallback：boundary `35`、crossover `315`、CEP `1`、spinodal `31`
- `xi=0, muB=150 MeV, alpha_T=1` 的 phase-guided crossover anchor：
  candidate runtime `198.2066818967471 MeV`，legacy `198.12550509377107 MeV`，差异
  `0.08117680297604579 MeV`。这是 reference 输入差异，不是 RS transport 数值 parity 结论。

## Numerical pilot 结果审计

source run `32684074876` 的 artifact 保留在 Actions；本地派生审计包为
`docs/analysis/pnjl/phase_reference/issue130_rs_transport_runtime_parity_v2_quality_audit/`。

- 5 个 warning key 在 candidate/legacy 中完全相同，原因均为 `tau_u_ubar_ratio_high`。
- 5 个点均 `converged=true`，输运字段 finite 且非负；两套结果的 phase kind/structure 一致。
- 参考输入差异通过 `mode_a_fixed_muB_phase_scaled` 的连续路径传播：先从 reference
  取 `T_phase_base_MeV`，再使用 `T=alpha_T*T_phase_base_MeV` 求 equilibrium 和 transport；
  因此 reference 不仅承担离散标签，也影响连续温度锚点与其后的平衡态/输运量。
- 该包支持“本限定 pilot 未见 candidate-specific 明显 transport regression”，不支持
  全域 numerical parity、旧 reference retirement 或 production promotion。

## CI 留痕

PR #261 合并前最新提交 `eba889feea0d422f20cb33e92ad62d9c8cf377fa` 的 CI 已于 2026-08-23
全部通过（12 项，无失败或跳过）：unit-smoke、integration-core、Julia/Python script smoke、
PNJL benchmark，以及 task-ledger、active-docs、docs consistency、data-output-path、
legacy-switch、dependency audit 和 advisory governance。

这只证明代码、证据 schema 与 solver-free consumer smoke 合同通过；不改变 numerical pilot
仍需独立 dispatch 和数值结果验收的边界。

## 非目标与回退

- 不删除或重写 `boundary.csv`、`cep.csv`、`crossover_dense.csv`、`spinodals.csv`。
- 不把 `unresolved`、`ambiguous`、`interpolated_noncertified` 行升级为 candidate certified。
- 不移除逐键 fallback，不移除显式 legacy rollback。
- numerical pilot 或后续 consumer 审计失败时，维持旧 reference/旧 production，提交具体失败字段和
  provenance；不放宽物理/数值容差。
