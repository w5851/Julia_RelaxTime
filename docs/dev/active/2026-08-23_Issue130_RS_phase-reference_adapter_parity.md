# Issue #130：RS transport phase-reference adapter parity 与限定重验

状态：active（post-repair audit 与 versioned result import 已通过，full candidate/legacy review 已接受）。PR #272 已创建且 11/11 checks 通过，等待作者审核与合并授权。runtime-switch PR #260、
solver-free parity PR #261、paired numerical pilot workflow PR #262 和 provenance repair PR #269
均已合并；RS numerical v2 已完成但仍保持 diagnostic-only。
本任务承接 RS transport 的 consumer parity；它不删除旧 reference，不覆盖旧 transport
production，也不把 solver-free smoke 当作数值 production 通过。

## 固定输入与边界

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- candidate source run：`32354095831`
- aggregate replay：`32451053476`
- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`
- legacy fallback：`data/reference/pnjl/legacy_phase_reference_v1/` versioned snapshot
- runtime view：strict candidate 的 certified-only 行，缺失/不合格键逐键 legacy fallback
- rollback：`--phase-reference-mode legacy`
- 本机不调用 equilibrium solver；本轮 post-repair audit 不重跑 C0/C1/C2、M4 solver 或 transport

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
- [x] 作者已审核限定 numerical evidence 与共同质量警告，并确认新 reference 效果符合预期；
  限定证据仍保持 diagnostic-only，不改写为全域 numerical parity。
- [x] phase-guided/Paper P1 等消费者的 adapter 路径已完成引用审计；PR #263 合并后，作者另行
  授权 `old-reference retirement` 与 RS transport 重验。retirement 只退出 canonical 根路径，
  保留 versioned fallback/rollback snapshot。
- [x] 复用历史 15+15 shard 拓扑完成 RS numerical v2：30 个唯一 run 成功，10 个重复 dispatch
  已取消；calculation SHA、workflow head、direct-coexistence `xi=±0.003` 合同均固定。
- [x] PR #269 合并后，用同一不可变 artifact 完成 post-repair solver-free provenance audit；
  verdict 为 `post_repair_audit_pass_diagnostic_only`，没有意外 artifact/hash、finite/converged、
  failed-point、重复键或 hard-gate 失败。
- [x] 作者审核 full candidate/legacy comparison 与同构审核图；作者确认结果符合预期。
  candidate 图已作为 `author_accepted_formal_layout` 导入正式 figure 路径；数值证据仍保持
  `diagnostic-only`。本次 versioned RS promotion/import PR 已在隔离 worktree 中生成
  `...prod_v2` result trees，保留 legacy fallback 和显式 rollback；PR #272 已提交，CI 11/11 通过，等待作者审核与合并授权。

- [x] solver-free source preflight 与 versioned result import 已完成：aggregate v4、30 个选定
  shard sidecar、candidate figure manifest 和旧 `prod_v1` tree hash 均通过校验；本次 import
  未调用 solver、未切换 runtime、未覆盖旧结果。
- [x] mode-A candidate result 已写入
  `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`：
  910 scan 行、38220 diagnostic 行，manifest SHA256
  `3d0a91170a62f570e119d1189d8e79526c008a1a77a94464772d53c3ca9e5f72`。
- [x] mode-B candidate result 已写入
  `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`：
  909 scan 行、38178 diagnostic 行，manifest SHA256
  `093ff05502f63515ac70b1a3e992971ecbc649c343bc7e591fa4fb1b8fc54921`。
- [x] 两套 result manifest 均记录 `artifact_status=imported_candidate`、
  `numerical_status=diagnostic_only`、`runtime_default_switch=false`、
  `source_solver_called=true`、`aggregate_replay_solver_called=false`、
  `import_solver_called=false`；质量警告和历史 sidecar hash 缺陷均保留。
- [ ] PR #272 的作者审核、合并后的 solver-free consumer smoke 尚未完成；PR checks 已 11/11 通过。
  在此之前不得把 `prod_v2` 设为默认 runtime、不得删除 `prod_v1`。

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
- 该包支持“本限定 pilot 未见 candidate-specific 明显 transport regression”；作者据此授权
  retirement 和后续 production 重验，但证据本身仍不等价于全域 numerical parity。

## 15+15 Numerical v2 与 Post-Repair Audit

source root：
`D:\Desktop\Julia_RelaxTime_issue130_artifacts\rs_sharded_production_v2_20260826`。

- aggregate：`aggregate_replay_20260826_v4`；post-repair audit：`post_repair_audit_20260826_v1`。
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；workflow head：
  `22874505877491754eed27519ad8a7b871c82571`。
- 30/30 unique shards 的 scan、finite/converged、failed points、重复键、hard gate 和 timing
  全部通过；mode-A 910 行、mode-B 909 行。
- 60 个 manifest mismatch 全部是历史 `effective_config.json`/`failed_points.csv` sidecar 写入顺序
  缺陷；PR #269 修复未来 producer，旧 artifact 不改写。mode-B `alpha_T=NaN` 是固定-T 路由的
  预期字段，不是非有限失败。
- 质量警告保留为诊断：mode-A `tau_u_ubar_ratio_high` 230 条，mode-B 240 条；候选/legacy
  phase-field mismatch 分别为 15 和 47 个共同键。它们进入作者审核，不自动判定为 solver regression。
- 当前 verdict：`post_repair_audit_pass_diagnostic_only`。本包只证明执行与 provenance 可审计，
  不晋升 RS production、不移除 fallback、不修改旧结果。

作者已确认 candidate/legacy 差异及同构图像。新增正式 figure 导入证据位于
`docs/analysis/relaxtime/issue130_rs_sharded_production_v2_post_repair_audit_v1/manifest.json`
的 `figure_formalization`，目标为两棵 `...prod_v2` 图树。该导入只复制 PNG 并写入逐图 hash，
`solver_called=false`、`production_write=false`；它不等价于数值 production promotion。

下一步是提交并审核独立的 versioned RS promotion/import PR，再在合并 SHA 上执行 solver-free
consumer smoke；只有 smoke 和作者授权均通过，才另立 runtime-default switch / old-reference
retirement 后续任务。

## CI 留痕

PR #261 合并前最新提交 `eba889feea0d422f20cb33e92ad62d9c8cf377fa` 的 CI 已于 2026-08-23
全部通过（12 项，无失败或跳过）：unit-smoke、integration-core、Julia/Python script smoke、
PNJL benchmark，以及 task-ledger、active-docs、docs consistency、data-output-path、
legacy-switch、dependency audit 和 advisory governance。

这只证明代码、证据 schema 与 solver-free consumer smoke 合同通过；不改变 numerical pilot
仍需独立 dispatch 和数值结果验收的边界。

## 非目标与回退

- 不重写 legacy 数值；canonical 根路径退役时以 byte-preserving snapshot 保留。
- 不把 `unresolved`、`ambiguous`、`interpolated_noncertified` 行升级为 candidate certified。
- 不移除逐键 fallback，不移除显式 legacy rollback。
- numerical pilot 或后续 consumer 审计失败时，维持旧 reference/旧 production，提交具体失败字段和
  provenance；不放宽物理/数值容差。
