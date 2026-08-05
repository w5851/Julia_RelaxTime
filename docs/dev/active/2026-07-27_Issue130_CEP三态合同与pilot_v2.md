# Issue #130：CEP 三态合同与窄窗口 pilot v2

创建日期：2026-07-27

状态：三态合同 PR #147 已 squash merge；pilot v2 数值 Actions、aggregate replay 和
evidence 导入已完成。hybrid shadow replay 与 Stage-C offline feasibility 已完成；随后
Maxwell tolerance contract、三交点候选 feasibility 与单点 endpoint refinement 已完成。
作者已接受有界零密度端点极限作为一阶证书；当前先沉淀 endpoint-limit evidence，再进入
公共 Maxwell/hybrid production integration，不启动全温区 phase-reference 或 transport。
PR #146 已 squash merge（merge SHA
`5071f6c80bb10c812358b855e8da2bde8a758f9d`；最终 Actions run `30199406478`，
9/9 matrix + aggregate success），v1 evidence 结论保持 `diagnostic_only`，不晋升
reference。三态合同 merge SHA 为
`0603e4df656ac610e12d51f9c810ec93be6e8f14`；pilot 分支为
`codex/issue-130-cep-narrow-pilot-v2`。

三态合同 PR：<https://github.com/w5851/Julia_RelaxTime/pull/147>，最终 head
`c418e1af`，merge 后远程分支已删除。

## 顺序与边界

1. [x] 合并并清理 PR #146，保留 v1 evidence 与 diagnostic-only 结论。
2. [x] 完成三态 CEP 合同实现与 focused 验证，并 squash merge PR #147。
3. [x] 从 `0603e4df` 创建 `codex/issue-130-cep-narrow-pilot-v2`，完成 runner、collector、
   plot、workflow 和 focused contract 验证；完整数值由 GitHub Actions 运行 cascade
   discovery、dense baseline 和独立 oracle。
4. [x] pilot v2 数值 run 与 aggregate replay 完成，证据导入
   `docs/analysis/pnjl/cep_narrow_pilot_v2/`；当前停在作者对 CEP 区间、原始 rho–mu
   曲线、准确度和性能证据的物理审核。
5. [ ] production integration 通过后，才重放全温区 phase-reference；reference 审核后再启动 transport。

## 三态合同验收

- [x] `confirmed_first_order`、`confirmed_monotone`、`ambiguous_near_critical` 切片语义。
- [x] `unknown_budget` 仅用于停止/诊断，预算超限不强制改写物理标签。
- [x] CEP 两端 evidence、实际 ambiguity interval 与 `temperature_resolution_target_MeV`。
- [x] `CEPResult` 新字段、旧构造兼容和 ambiguous/not-found 单点空值合同。
- [x] Dense CEP CSV 新旧 schema 兼容、每个 xi 保留一行、merge/validator 支持。
- [x] phase summary/report/manifest/grid convergence 消费 `result_status`。
- [x] research interpolate/direct 与 production 采用同一三态语义；不接入 cascade、不改变 equilibrium/mixed-meson/non-fixedmu 语义。
- [x] focused regression/validation、CLI 和治理检查全部通过。
- [x] 提交、推送并创建 ready PR #147，更新 Issue #130，停在作者审核。

## 三态合同本地证据

- Julia unit：CEP detector `42/42`、production helpers `25/25`、phase artifacts
  `28/28`、phase grid convergence `20/20`、dense builder `23/23`。
- Julia integration/regression/validation：phase core `14/14`、artifact promotion
  `14/14`、old/new consistency `10/10`、mode compare `26/26`、phase regression
  `49/49`、literature targets `5/5`。
- Python：dense merge/validator `41 passed`；修改脚本 `py_compile` 通过。
- Governance：docs/active-docs/script-entry/model-entry/solver-leakage/PNJL migration/
  data-path/agent-instructions、unit-skip、relaxtime governance、profile matrix、
  dependency、skill governance 和 Claude skill sync 全部通过。
- 兼容边界：默认 production 仍要求两层 rho geometry 才确认一阶/单调；显式
  `rho_geometry_convergence=false` 仅保留旧单层一阶边界候选，不提升单层
  `no_s_shape` 为 `confirmed_monotone`。

## Pilot v2 固定口径

- 新 schema/tag：`cep_narrow_pilot_v2`，不覆盖 v1、canonical reference 或 C0/C1/C2。
- `xi=-0.5,0,0.5`；先运行三个 cascade discovery，再生成不可变验证窗口，随后六个独立
  dense/oracle 验证 job 和 aggregate 审计。
- dense `Δrho=0.00625`；oracle `Δrho=0.003125` 且独立重跑局部温度搜索；端点搜索目标保留
  `0.125 MeV`，不把它当作实际 ambiguity 宽度。
- `xi=0.5` 最多向高温侧扩展至 canonical `+32 MeV`；没有双端证据则为
  `oracle_inconclusive`，不扩成全相图。

## Pilot v2 实现记录

- [x] 新增两层 rho 证据 runner：cascade `0.05→0.025`、dense `0.0125→0.00625`、
  oracle `0.00625→0.003125`；保留 exact `(T,xi,rho)` memoization 和 request-scoped
  solver telemetry。
- [x] 新增三态双端搜索、`validation_windows.json` 冻结脚本、v2 collector/plotter，
  以及 v1/v2 workflow 分支和 aggregate replay 合同。
- [x] 新增 Python v2 contract tests；`py_compile` 和 v1/v2 Python focused tests 通过。
- [x] Actions numerical run `30339094251` 完成（3 discovery、freeze、6 validation、
  aggregate）；calculation SHA `9d325b62ef366f1de70cd3b2a3555e71bd12664a`。
- [x] aggregate replay `30342409727` 使用 collector 修复 SHA
  `d2f8b1d0c2fb7a6cec7d21a1076a1d24639c9ca5` 完成；`status=pass`、
  `oracle_status=stable`、`cascade_status=within_oracle`，Actions critical path
  `1247 s`、runner-minutes `91`。
- [x] 导入 `docs/analysis/pnjl/cep_narrow_pilot_v2/` 的聚合表、图像、manifest、
  plot manifest、curve index 和审计文档；原始 `curve_points.csv` 保留在 Actions/
  本地 artifact，不复制进仓库。
- [ ] 物理审核通过后才进入 production integration、phase-reference replay 和 transport。

本机无法可靠加载完整 `src/models/Models.jl`（Julia 进程在预编译阶段被环境终止），
因此本轮仅记录 parser/smoke/Python 治理结果；Julia runtime contract 需在 CI Actions
中验证，不在本机执行完整 PNJL 数值。

## 后续 production integration（当前分支）

当前从 `main@a8dc06ec077bb4f38c34dd162467be140631613e` 开发
`codex/issue-130-cep-cascade-production-integration`。本阶段只加入显式
`rho_support_cascade` 与 request-scoped point cache/telemetry，以及独立 9-job
production shadow workflow；`uniform_nested` 默认路径、旧 evidence、旧 reference、
C0/C1/C2 artifacts 和 transport 均保持不变。shadow 的物理 verdict 仍需 Actions
证据与作者审核；未得到 `full_cascade_candidate` 前不转 ready、不晋升 reference，
低温失败默认停在 `hybrid_required` 等待判断。

### Production integration shadow（当前 Draft 阶段）

- [x] 从 `main@a8dc06ec077bb4f38c34dd162467be140631613e` 提取 opt-in
  `rho_support_cascade`、request-scoped exact `(T,xi,rho)` cache 和 solver telemetry；
  默认 `uniform_nested` 保持原路径。
- [x] 新增 9-job workflow
  `.github/workflows/pnjl-cep-production-shadow.yml`，固定 calculation SHA、
  `p_num=24/t_num=8`、`rs_reduced_adaptive` 和三种 rho 分辨率口径；不写 reference。
- [x] 新增 collector/plotter、curve index、provenance/hash/重复键/finite 合同与
  Python contract tests；完整 PNJL 数值尚未在本机运行。
- [x] PR #149 已合并并从 main 触发 shadow；source run 与 replay 均保留 immutable
  calculation SHA、provenance 和 raw curve index。
- [x] shadow 主 verdict 为 `oracle_inconclusive`；不自动落地 hybrid、不重放 C0/C1/C2、
  不启动 transport。

## Hybrid shadow v2 诊断收口（2026-07-31）

- [x] PR #153 已合并为 `main@fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc`，并触发
  hybrid shadow run `30601857594`；9 个数值 job、aggregate 和 provenance 均成功。
- [x] 已下载并审计 source artifact：58,875 个曲线点 finite/converged，重复键为 0，
  solver failure/fallback/retry 为 0；当前仍为 diagnostic-only。
- [x] 当前 run 暴露三类问题：oracle 低温必需 anchors 不稳定、Stage C 分类与 oracle
  存在差异、hybrid solver work 高于 dense；另有 `mu_transition_MeV` 和 Actions
  provisional/final snapshot 合同缺陷。
- [x] PR1 `codex/issue-130-hybrid-shadow-evidence-replay`：修正字段、provenance、
  strict oracle gate 和 replay；检查通过后自动合并，并从 source run `30601857594`
  生成 final aggregate evidence。
- [x] PR2：利用现有完整曲线做 Stage C 全局曲线合并、多 Maxwell 候选和 adaptive-cap
  离线 feasibility；未得到 `feasible_candidate` 前不修改 production。
- [ ] PR3：仅在 feasibility 通过后实现自适应 Stage C，自动合并后从新的 main SHA
  运行 targeted/full shadow；通过并经作者审核前不启动 C0/C1/C2。

## Hybrid shadow replay 与 Stage-C offline feasibility（2026-07-31）

- [x] PR1 evidence/replay 修正已 squash merge：PR #154 merge SHA
  `04a84ed4b98f990248cadb8b4a474d6b83a2813c`；历史 Stage-C support 标签兼容修正
  PR #155 merge SHA `982fed04783c64c431b1a4745e88ab43e2040dc3`。
- [x] source run `30601857594` 已用 aggregate replay `30609095886` 收口；postprocess
  SHA 为 `982fed04783c64c431b1a4745e88ab43e2040dc3`，evidence state 为 `final`，主
  verdict 为 `oracle_inconclusive`。保留 5 个分类差异、6 个性能超限项和 11 组旧
  Stage-C 标签 compatibility warnings；未重算 PNJL。
- [x] PR2 离线 Stage-C feasibility 已用 replay CSV 完成 deterministic replay：完整
  Stage-B 0.00625 global reference 与局部 0.003125 点合并，枚举稳定 `+→−→+`
  候选，并测试 caps `48/96/160/224`；solver_called=false。
- [x] PR #157（merge SHA `ded56b776557528c8bdfa5bf9f41ff30f73a7dda`）收紧跨分辨率
  candidate overlap 与 position/density/area geometry gate；重放结论保持不变。
- [x] 当前 PR2 预期 verdict 为 `maxwell_candidate_inconclusive`：所有 cap 的成本
  估计均低于 dense 的 15384 unique solves，但有两个弱 S/CEP 邻域 first-order
  oracle anchors（`xi=-0.5,T=147.0947265625` 与 `xi=0.5,T=106.9599609375`）
  无法形成唯一稳定 topology；没有选出 `selected_policy`。
- [ ] 因 feasibility 未通过，不创建 production adaptive Stage-C PR，不启动
  targeted/full shadow、C0/C1/C2、reference replay 或 transport；等待作者判断弱 S
  与 Maxwell 候选的后续物理/数值方案。

## Stage-C certificate feasibility v2（2026-07-31）

- [x] 从 `main@b02336802fdb6e15a0592cd152049e82941afcd6` 创建
  `codex/issue-130-stagec-certificate-feasibility-v2`。
- [x] 新增 solver-free Julia replay
  `scripts/analysis/pnjl_cep_hybrid_stagec_certificate_feasibility_v2.jl`：直接复用
  `PhaseCore`、Maxwell construction 和 phase geometry comparison；Stage-B semantic
  grid 不再免费计入成本，oracle row 不参与选点。
- [x] 新增 density-support/drawdown 特征路由和 `12/24/48/96/160` cap frontier；
  v2 不修改 equilibrium solver、`detect_s_shape`、默认 uniform/cascade 或历史 evidence。
- [x] 两个作者确认的一阶点重算为 `confirmed_first_order`，两个高温 control 重算为
  `confirmed_monotone`；说明此前弱 S 假阴性来自 v1 离线候选规则，而非 production
  detector 在这四个点上的结果。
- [x] v2 成本 frontier 为 `6588–6901` unique solves，dense reference 为 `15384`；
  Maxwell/geometry 和候选唯一性 gate 均通过。
- [x] v2 verdict 为 `oracle_inconclusive`：5 个低温 oracle-ambiguous 点被离线
  Stage-C 证书确认，已写入 `deep_oracle_required.csv`，不自动晋升为物理结论。
- [x] 新 evidence 写入
  `docs/analysis/pnjl/cep_hybrid_stagec_certificate_feasibility_v2/`，保留 v1
  evidence 原样；代表性曲线图、作者裁决 provenance、manifest 和 claim ledger 已生成。
- [x] focused Julia helper `13/13`、Julia parser、Python `py_compile` 和
  `git diff --check` 通过；完整 PNJL solver 未在本机调用。
- [ ] 因当前 verdict 为 `oracle_inconclusive`，暂不创建或合并 adaptive production PR，
  不运行 targeted/full shadow，也不启动 C0/C1/C2、reference replay 或 transport；下一步
  仅在作者允许后，对 `deep_oracle_required.csv` 中 5 个点运行 `0.003125→0.0015625`
  deep oracle。

## Deep oracle（作者已批准，2026-08-01）

- [x] 固定范围为 `(-0.5,20)`、`(0,5)`、`(0,20)`、`(0.5,5)`、`(0.5,20)` 五个低温
  anchor；不扩展到全相图。
- [x] 新增独立 workflow `.github/workflows/pnjl-cep-deep-oracle.yml`，每个 anchor
  一个独立 job，使用 `0.003125→0.0015625` 的 independent-oracle rho 分辨率。
- [x] workflow 使用独立的 workflow-head SHA 与 calculation SHA provenance；不写 reference。
- [x] 修复 `nothing` 无法由 CSV.jl 输出的问题（PR #160，merge
  `592374ff23c101575587d1164e6aeca9a7231fc1`），保留 JSON `null` 语义并通过 focused CI。
- [x] 从新 main SHA 触发并完成 GitHub Actions deep-oracle run `30676440627`；workflow-head
  SHA 为 `592374ff23c101575587d1164e6aeca9a7231fc1`，calculation SHA 仍为
  `fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc`。
- [x] 下载并审计五个 job 与 aggregate artifact 至
  `D:\Desktop\Julia_RelaxTime_issue130_artifacts\cep_deep_oracle_20260801`：12,462 条
  curve rows、曲线 key 唯一、所有最终 equilibrium 解 finite/converged、solver failures=0；
  每点 `unique_solves=2561`、`point_requests=3842`、`cache_hits=1281`，哈希和双 SHA
  provenance 全部一致，`reference_write=false`。
- [x] 结果保持 diagnostic-only：5 个点的 `cep_accuracy.csv` 均为 `ambiguous`/`found=false`，
  aggregate physical verdict=`author_review_required`；`xi=0,T=5` 的局部 slice 虽为
  `confirmed_first_order` 且 geometry 通过，但单点不构成 CEP resolved 证书。
- [ ] deep-oracle 结果不自动创建 adaptive production PR，不启动 targeted/full shadow、
  C0/C1/C2、reference replay 或 transport；等待作者对 ambiguous/geometry 证据和后续
  Stage-C 方案作物理判断。

## Maxwell tolerance-contract feasibility（2026-08-02）

- [x] 从 `main@b91158795c0b7d65e319c9c09a67915f3eef28a1` 创建
  `codex/issue-130-maxwell-tolerance-feasibility`；本阶段只做 solver-free replay。
- [x] 新增 `scripts/analysis/pnjl_maxwell_tolerance_feasibility.jl`，读取 deep-oracle
  aggregate `30676440627`，按实际 `(xi,T,rho)` 曲线重建 coarse/fine 层，显式传递
  `tol_area` 到 Maxwell bisection，并独立记录 outer geometry gate。
- [x] 扫描内部 area tolerance `1e-4/5e-5/1e-5/5e-6`。五个低温 anchor 在严格
  `5e-6` 下均保持两层 `+→−→+` topology、Maxwell residual、position/density geometry
  和迭代预算通过；`1e-5→5e-6` 未发生候选切换。
- [x] 输出 `docs/analysis/pnjl_maxwell_tolerance_contract_feasibility_v1/`，含逐层
  frontier、逐点 summary、输入曲线索引、manifest、claim ledger 和 audit；
  `solver_called=false`，不复制原始曲线、不修改 production/reference/history。
- [x] focused Julia `14/14`、Python `2 passed`、parser 和 `git diff --check` 通过。
- [ ] 当前 feasibility verdict 为 `feasible_candidate`；下一步从新的 main merge SHA
  实现 production tolerance contract（内部 solver tol 派生自 active acceptance tol，
  不新增独立用户参数），再只对这五点通过 Actions 重验。此前不启动全相图、C0/C1/C2、
  reference replay 或 transport。

## Production Maxwell tolerance contract（后继本地分支，2026-08-02）

- [x] 在 feasibility 通过后创建 `codex/issue-130-maxwell-tolerance-contract`，没有改写
  feasibility 分支或历史 evidence。
- [x] 新增内部派生合同：`tol_solver = 0.1 × minimum(active acceptance tolerances)`；
  `area_tol_good` 始终 active，rho/temperature geometry gate 按启用状态加入。
- [x] production 的 cascade、memoized uniform、hybrid 和默认 rho refinement 均显式传递
  派生 `tol_area`；classification 与密度提取复用同一个 `MaxwellResult`，solver tolerance
  未满足时不发布 valid 结果；默认 `maxwell_construction` API 和 uniform 路径保持兼容。
- [x] `maxwell_solver_tol`、派生因子及 active gates 写入 config snapshot、config hash
  和 diagnostics；API 文档补充 inner solver / outer geometry gate 的职责边界。
- [x] focused contract 11/11、CEP detector 42/42、PhaseCore 12/12、production helpers
  33/33、phase geometry 20/20、Python contract 2/2；完整 PNJL 未在本机运行。
- [x] production contract PR #162 已 squash merge，merge SHA 为
  `467be1fce847a9c991ec362c3335be07fccbe604`；默认 API、uniform 路径、equilibrium
  solver 与历史 evidence/reference 均未改写。内部 solver tolerance 由 active acceptance
  gates 派生，不新增用户可调独立容差。

## Maxwell tolerance-contract 五点 Actions revalidation（2026-08-02）

- [x] 从 `main@467be1fce847a9c991ec362c3335be07fccbe604` 触发固定五点 run
  `30730990835`：<https://github.com/w5851/Julia_RelaxTime/actions/runs/30730990835>。
  tag 为 `cep_maxwell_tolerance_contract_v1_20260802`；workflow head SHA 与
  calculation SHA 均为 `467be1fce847a9c991ec362c3335be07fccbe604`。范围严格为
  `(-0.5,20)`、`(0,5)`、`(0,20)`、`(0.5,5)`、`(0.5,20)`，rho 为
  `0.003125→0.0015625`；未启动全相图或其他 production。
- [x] 5 个数值 job 与 aggregate 全部 success，无需 failed-job rerun。artifact 已下载至
  `D:\Desktop\Julia_RelaxTime_issue130_artifacts\cep_maxwell_tolerance_revalidation_20260802`。
  aggregate 与 job manifest 的文件 SHA 全部匹配；`curve_points.csv` 共 12,462 行，
  `(xi,method,T,rho_level,rho)` 重复键为 0，曲线无 NaN/Inf。
- [x] 五个切片最终 equilibrium 解全部 finite/converged；solver failure、fallback、
  governed rescue、retry、exception 和 nonconverged attempt 均为 0。每点成本守恒为
  `unique_solves=2561`、`point_requests=3842`、`cache_hits=1281`；总计
  `unique_solves=12805`、`point_requests=19210`、`cache_hits=6405`。
- [x] config/diagnostics 实际记录 `maxwell_solver_tol=5e-6`，其来源为
  `0.1 × min(1e-4,5e-5)`；外层 position/density/Maxwell-area gate 仍为
  `0.025 MeV/0.0025/5e-5`。五个切片均为 `confirmed_first_order`，geometry 均通过；
  `maxwell_area_gate` 最大为 `4.8444e-6`，position error 最大为 `0.009704 MeV`，
  density error 最大为 `0.001414`，均在 gate 内。
- [x] 与此前 deep-oracle artifact 对齐后，12,462 个共同曲线键的 `muq`、pressure、
  residual、iterations、converged/finite 状态和 sampling role 均未发生数值漂移；本次
  变化来自 Maxwell tolerance contract 的后处理/证书判定与 provenance，不是 equilibrium
  solver 重算漂移。
- [x] `cep_accuracy.csv` 中五个单点仍为 `result_status=ambiguous`、`found=false`，
  原因是每个 job 只有一个温度，缺少高温 confirmed-monotone 端点；aggregate physical
  verdict 为 `author_review_required`。因此本 run 只证明五个固定切片的 Maxwell/geometry
  证书改进，不构成 CEP resolved，也不晋升 reference。
- [x] revalidation 结论为 candidate/diagnostic-only：保持不启动 C0/C1/C2、C3/O1、
  formal production、phase-reference replay/promotion 或 transport。下一步由作者决定是否
  进行更大范围或温度两端验证。

## Stage-C tolerance-contract solver-free replay（2026-08-02）

- [x] 新增版本化 replay 脚本
  `scripts/analysis/pnjl_cep_hybrid_stagec_tolerance_replay.jl`，并将 v2 replay 的
  Maxwell 语义拆为：CEP acceptance `area_tol_good=1e-4/area_tol_bad=5e-4`、内部
  `maxwell_solver_tol=5e-6`、外层 rho geometry area gate `5e-5`。已有 v2 evidence 不改写。
- [x] 输入严格锁定 base final aggregate `30601857594`（calculation
  `fde2b929…`）与五点 tolerance revalidation aggregate `30730990835`（calculation
  `467be1fc…`）；两套 manifest/hash、曲线键、finite/converged 和 provenance 均验证通过，
  `solver_called=false`，`reference_write=false`。
- [x] 五点重验局部证书全部通过：`(-0.5,20)`、`(0,5)`、`(0,20)`、`(0.5,5)`、
  `(0.5,20)` 均为两层 `confirmed_first_order`，geometry 通过，最大 area gate
  `4.6916e-6`，低于外层 `5e-5`。
- [x] replay 不能把两个 calculation SHA 的所有 curve row 静默视为同一数据集：共同
  positive-rho 键 5,978 个，`muq` 最大差 `1e-6 MeV`；rho=0 的 5 个端点存在最大
  `4.54764 MeV` 的非唯一根差异，已作为 provenance warning 单独记录并由 tolerance
  artifact 覆盖。旧/新共同点还有 127 个 solver iteration-count 差异，但 finite/converged、
  pressure 和 residual 均一致。
- [x] 24-anchor Stage-C frontier 在 caps `12/24/48/96/160` 下的模拟 unique solves
  为 `6588/6797/6972/6972/6972`，均低于 dense `15384`；但仍有 3 个 oracle-ambiguous
  anchor 被 hybrid 证书确认，以及 5 个旧 oracle/hybrid 分类 mismatch，因此主 verdict
  保持 `oracle_inconclusive`，没有 `selected_policy`。
- [x] 新 evidence 写入
  `docs/analysis/pnjl/cep_hybrid_stagec_tolerance_replay_v1/`，包含 replay、candidate、
  cost frontier、五点 tolerance certificate、curve identity audit、双 run provenance、
  manifest 和 claim ledger；不复制原始曲线。
- [ ] 因 replay 未得到 `feasible_candidate`，不创建 adaptive Stage-C production PR，
  不运行 targeted/full shadow，不启动 C0/C1/C2、phase-reference replay、reference
  promotion 或 transport；当前需要作者判断是否接受新的 calculation SHA 下重新做 24-anchor
  shadow，或继续只审计五点 tolerance 证据。

## Full 24-anchor shadow revalidation（2026-08-02）

- [x] 作者批准接受新的 calculation SHA 后重做完整 24-anchor shadow；从 `main` 触发
  workflow run `30736597984`：<https://github.com/w5851/Julia_RelaxTime/actions/runs/30736597984>。
- [x] workflow head SHA 为 `1cd84f02363ce0346b984dd524575579d8ad5d55`；数值 checkout
  calculation SHA 固定为 `467be1fce847a9c991ec362c3335be07fccbe604`，tag 为
  `cep_hybrid_stagec_tolerance_full_20260802`。矩阵为 3 个 xi ×
  `production_hybrid`/`memoized_dense`/`independent_oracle`，每个 xi 包含 8 个 anchors。
- [x] 运行期间只监控 9 个 numerical jobs 与 aggregate；失败只重跑 failed jobs。不得启动
  adaptive production、C0/C1/C2、phase-reference replay 或 transport。

- [x] 9 个 numerical jobs 全部 success；aggregate 首次因 Stage-C legacy label contract、
  oracle 低温必需锚点和 hybrid 成本 gate 失败。按约定仅执行
  `gh run rerun 30736597984 --failed`，没有重算数值 jobs；aggregate 重跑仍为 deterministic
  failure。
- [x] 本地下载并审计全部 9 格 artifacts：calculation SHA、finite/converged、curve keys 和
  telemetry 文件完整。`--legacy-replay` solver-free collector 将同一数据收口为
  `oracle_inconclusive`，并保留 6 个 Stage-C 标签 compatibility warnings；性能仍为
  hybrid `15954` requests / `15953` unique solves / `303.11 s`，dense `15384` / `15384` /
  `270.26 s`，classification mismatch 为 `(xi=0,T=60)`，oracle 低温必需 anchors 为
  `(-0.5,5)`, `(-0.5,20)`, `(0,5)`。

- [x] 为满足严格双 SHA provenance，另从临时 ref
  `codex/issue-130-shadow-calculation-467`（指向 calculation SHA）触发 run
  `30737739707`：<https://github.com/w5851/Julia_RelaxTime/actions/runs/30737739707>；
  workflow head 与 calculation SHA 均为 `467be1fce847a9c991ec362c3335be07fccbe604`。
  9 个 numerical jobs 全部 success，aggregate 及 failed-only rerun 均 deterministic failure。
- [x] 严格 run 的 solver-free `--legacy-replay` 收口仍为 `oracle_inconclusive`：同一
  3 个 oracle 低温锚点、`(xi=0,T=60)` 分类 mismatch 和 hybrid 成本超 dense gate；
  hybrid `15954` requests / `15953` unique solves，dense `15384` / `15384`，无 fallback/retry
  恶化。严格 run artifacts 已下载至
  `D:\Desktop\Julia_RelaxTime_issue130_artifacts\cep_hybrid_stagec_tolerance_full_calcsha_20260802`。
- [ ] 因 oracle、classification 和 performance gates 未通过，保持 diagnostic-only；不创建
  adaptive production PR，不启动 C0/C1/C2、phase-reference replay 或 transport，等待作者
  决定是否先做低温/`(xi=0,T=60)` 曲线物理审计及针对性算法修正。

## Stage-C 离散极值 guard production revalidation（2026-08-03）

- [x] 在现有 `main@dce7936888e57ddf2b3d6231ca4ad159fb552728` 上完成 solver-free
  extrema-guard feasibility PR #165；selected cap 为 `12`，模拟 unique solves 为
  `12148`，低于 memoized dense `15384`。`(-0.5,5)` 缺少左侧严格外点，按合同保持
  `ambiguous_near_critical`，不硬插入 `(0,0)`。
- [x] production 分支
  `codex/issue-130-stagec-extrema-guard-production` 实现 opt-in
  `rho_support_hybrid` 的 `extrema_outer_samples_v1`：Stage-C 始终保留完整 Stage-B
  全域曲线，只在两个 μ 极值外侧首个严格 Stage-B 采样点形成的 guard 内按 Stage-B
  特征排序补 `0.003125` 点；默认 uniform/cascade、equilibrium solver、旧 evidence
  和 reference 不变。
- [x] 新增 `RhoHybridVerificationConfig`、guard μ/状态/来源诊断、scope-aware shadow
  runner/collector/workflow 合同；targeted 结果使用 `targeted_hybrid_candidate`，完整
  24-anchor 结果使用 `full_hybrid_candidate`。
- [x] 本地 deterministic/focused 验证：production helper `45/45`，hybrid collector
  Python `7/7`，相关 Python shadow 合计 `18/18`，parser、docs consistency、active-docs、
  script entrypoints、Models entry contract、solver leakage、relaxtime governance、
  data-output-path 和 unit-skip 全部通过；未在本机运行 PNJL 完整数值。
- [x] production PR #166 已 squash merge，merge SHA 为
  `3217bed3635574f00c04cbee75e843b4c49451db`。随后 collector 温度序列化修正 PR #167
  merge SHA 为 `8c49742f6f17061041405384c02a334ec54f7be6`，replay provenance 修正 PR #168
  merge SHA 为 `782e08fdf1a5d7e041811d6ec947866ff70fbd47`；这些 PR 均未改变 calculation
  SHA 或数值输入。
- [x] 从 calculation SHA `3217bed3635574f00c04cbee75e843b4c49451db` 运行 targeted shadow
  `30805637032`，9 个 numerical jobs 全部 success；aggregate replay
  `30808473818` 仅下载并重建 evidence，postprocess SHA 为
  `782e08fdf1a5d7e041811d6ec947866ff70fbd47`。双 SHA、source run、曲线哈希、唯一键和
  finite/converged 审计通过；targeted 主 verdict 为 `oracle_inconclusive`。
- [x] targeted oracle 仅在 `(-0.5,5)`、`(-0.5,20)`、`(0,5)` 未形成必需一阶证书；hybrid
  fixed-rho `9407`、dense `11538`、oracle `23058`，fallback/retry 均为 `0`。
- [x] 为严格冻结这 3 个点，deep-oracle scope PR #169 已 merge，workflow head SHA 为
  `3777810c2a58683277882fe6178683ba28558779`；required-three deep run
  `30809754119` 使用 `0.003125→0.0015625` 独立 oracle，3 个 jobs 全部 success。
  `(-0.5,20)` 与 `(0,5)` 两层均为 `confirmed_first_order` 且 geometry 通过；
  `(-0.5,5)` 两层均为 `maxwell_solver_tolerance_not_met`，保持
  `ambiguous_near_critical`。deep aggregate 为 `author_review_required`，不构成自动
  晋升。
- [ ] 因 `(-0.5,5)` deep oracle 仍未闭合，当前 verdict 保持 `oracle_inconclusive`；不运行
  完整 24-anchor shadow，不创建 evidence promotion PR，不启动 C0/C1/C2、phase-reference
  replay、reference promotion 或 transport。下一步需作者决定是否对该低温点进行物理/数值
  审计（尤其 Maxwell tolerance 与候选区间），或接受其 ambiguous 语义后再调整 gate 计划。

## Maxwell 三交点候选与 endpoint-limit 证书（2026-08-05）

- [x] PR #171 已 squash merge，merge SHA 为
  `df1fcdece7bc0888dd57c465c0015828743596c5`。该 PR 加入 solver-free 三交点 Maxwell
  candidate feasibility、单点 endpoint refinement runner/collector/workflow 和 focused tests；
  没有修改 equilibrium solver、reference 或 transport。
- [x] 从 workflow head `df1fcdec…`、固定 calculation SHA
  `3217bed3635574f00c04cbee75e843b4c49451db` 运行单点 Actions run
  `30980094983`，范围严格为 `(xi,T)=(-0.5,5 MeV)`。1281 个基础点加 12 个近零 targeted
  点全部 finite/converged，solver failure、fallback、retry 和 exception 均为 0。
- [x] 13 个 refinement levels 始终保持唯一三交点 Maxwell candidate；最终
  `mu_M=331.5739309844038 MeV`、`rho_q=2.8043699741472983`、面积残差
  `-2.9237516382062305e-7`，左共存密度被约束为
  `0 <= rho_h < 7.62939453125e-7`。source verdict 保持
  `candidate_only_endpoint_inconclusive`，不改写历史 artifact。
- [x] 作者明确接受上述有界零密度端点极限作为一阶相变证据，不要求严格正密度下界。
  三态分类因此仍为 `confirmed_first_order`，内部证书类型定义为
  `endpoint_limited_first_order`；兼容 scalar `rho_hadron` 使用 `0.0`，并同时保存真实上界，
  不新增第四种物理状态，也不硬编码特定 `(xi,T)` 例外。
- [x] 当前分支 `codex/issue-130-maxwell-endpoint-limit-contract` 以 solver-free 方式生成
  `docs/analysis/pnjl_maxwell_endpoint_limit_contract_v1/`，固定校验 source run、双 SHA、两个
  manifest SHA、三交点唯一性、逐级 bracket 二分、Maxwell/geometry、finite/key uniqueness
  和成本守恒；derived verdict 为 `endpoint_limited_first_order_candidate`。
- [ ] endpoint-limit evidence PR 合并后，另开公共 Maxwell/hybrid production PR，将该证书
  作为通用诊断合同落地，再按 targeted → 必要 deep oracle → full 24-anchor shadow 顺序重验。
  在 `full_hybrid_candidate` 与作者曲线/成本审核前，仍不启动 C0/C1/C2、phase-reference、
  reference promotion 或 transport。

## Maxwell Endpoint-Limit Production Integration（2026-08-05）

- [x] PR #173 已 squash merge；当前 production 分支从 `main@a80cf5573adec9a6061cc433dd1d1083d7b82b24`
  创建，endpoint-limit evidence 与历史 pilot 保持不可变。
- [x] 公共 `PhaseCore.maxwell_construction` 已改为唯一三交点候选合同：枚举全部去重交点，
  拓扑间隙重置面积变号，支持多 candidate 诊断，并在 `max_iter` 耗尽时保持未收敛。
- [x] hybrid 增加 opt-in `bounded_zero_density_v1` endpoint route：Stage-B 唯一 endpoint-dependent
  candidate 且右外点存在时加入 `rho=0.003125` anchor，最多 12 次 `[0, anchor]` 二分；成功时
  发布 `endpoint_limited_first_order`、`rho_hadron=0` 和真实上界/插值诊断，否则保留 ambiguous。
- [x] config snapshot/hash、production sweep diagnostics、shadow v3 schema 与 endpoint gate 已同步；
  默认 uniform/cascade、equilibrium solver、reference 和 transport 均未改变。
- [ ] focused CI 全绿后 squash merge production PR；随后用该 merge SHA 运行 endpoint-limit targeted
  shadow，必要 deep oracle 仅限已批准集合，targeted 通过后再运行 full 24-anchor shadow。

## Endpoint-Local Geometry 合同 v2（2026-08-05）

- [x] 保持 v3 数值 artifact、`ceec2295…` calculation SHA 和历史 evidence 不变；collector 增加
  `run_mode=numerical|aggregate_replay`，数值 aggregate 明确为 provisional，只有带认证的
  source-run replay 才生成 final Actions 成本快照。
- [x] deep overlay 只在 standard oracle 为 ambiguous 的 approved required-three 点生效；输出
  `standard_oracle_status`、`deep_oracle_status`、`final_oracle_status` 和 `oracle_source`，不把
  deep 标签用于 support/route 选择。
- [x] solver-free feasibility 输出到
  `docs/analysis/pnjl_cep_endpoint_local_contract_feasibility_v2/`。三个低温点均有唯一三交点
  与实际右外支 bracket，保留完整 Stage-B 曲线，左侧只做 active-bracket midpoint replay；
  既有 position `0.025 MeV`、density `0.0025`、area `5e-5` 门禁均未放宽。derived verdict 为
  `feasible_candidate`，覆盖 targeted 18 + approved deep required-three；完整 24-anchor shadow
  仍是后续硬要求。
- [x] 从 feasibility merge SHA `650816e9…` 创建
  `codex/issue-130-endpoint-local-production-v2`，将
  `three_crossing_endpoint_local_v2` 作为显式 hybrid policy 落地；保留旧
  `bounded_zero_density_v1`、完整 Stage-B 曲线与默认 uniform/cascade 语义，并将 endpoint
  route/bracket、证书类型、补点 trace 与成本字段写入 shadow schema v4。
- [ ] PR2 focused CI 全绿后 squash merge；随后以新的 production merge SHA 串行运行 focused →
  必要 deep overlay → targeted 18-anchor → full 24-anchor shadow。full evidence 经作者审核前，
  不启动 C0/C1/C2、phase-reference 或 transport。
