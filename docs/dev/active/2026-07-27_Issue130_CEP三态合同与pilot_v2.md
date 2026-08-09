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

## Endpoint-local v4 focused replay 合同修正（2026-08-05）

- [x] PR #175–#180 已合并；当前 `main` 为
  `a93d6de19cd0a76c44eb488515c4596adcc74e4d`，focused calculation SHA 为
  `aded4a51d25fec6b6c121ce78ed3f0cf7f327b5a`，numeric run 为 `31027329971`。
- [x] focused numeric 的三个问题点均产生预期 endpoint 证书：`(-0.5,5)` 为
  `endpoint_limited_first_order`，`(-0.5,20)` 与 `(0,5)` 为
  `endpoint_local_geometry_first_order`；hybrid/dense/oracle 的 solver 与成本数据均已生成。
- [x] replay `31031162303` 已进入物理/合同 gate；失败原因为 endpoint-local
  `support_low/high` 错把最终收缩 bracket 当作声明 support，导致 anchor/midpoint 被判定为越界，
  不是 solver、Maxwell 或容差失败。
- [x] 修正 support 字段为初始左 bracket 与全部 endpoint anchor/midpoint 的 envelope；最终收缩值
  继续只写入 `endpoint_lower_bound/endpoint_upper_bound`（PR #181）。修正后才继续 targeted/full。

## Endpoint-local v4 full shadow 与 evidence 收口（2026-08-06）

- [x] PR #181 修正 endpoint-local `support_low/high` 为初始 bracket 与实际 anchor/midpoint
  的 envelope；PR #182 修正 deep artifact pattern；PR #183 区分 workflow head SHA 与
  calculation SHA。当前 workflow head 为 `6ceaa1ea…`，固定 calculation SHA 为
  `8fb86255c7894004cb1c52bbb03f9f53a4828411`，旧 endpoint-limit evidence 保持不可变。
- [x] focused numeric `31064585632` 与 required-three deep oracle `31065343978` 完成；
  focused final replay `31066464218` 的 verdict 为 `focused_hybrid_candidate`。
- [x] targeted numeric `31066589411` / replay `31067875491` 通过，verdict 为
  `targeted_hybrid_candidate`；随后 full numeric `31067984700` 与 final replay
  `31068641378` 通过，verdict 为 `full_hybrid_candidate`。本次 full replay 是同一
  calculation SHA 的 final aggregate，不调用本机 solver。
- [x] full matrix 共 72 条 slice rows：hybrid 与 oracle 均为 18 个
  `confirmed_first_order`、6 个 `confirmed_monotone`；oracle/classification/endpoint/
  coverage/performance errors 均为空。54659 条 raw rho–mu 曲线 finite/converged 且 key
  唯一；fallback/retry 均为 0。
- [x] 三个 endpoint 证书均通过：`(-0.5,5)` 为
  `endpoint_limited_first_order`，support `[0,0.00625]`、活动上界
  `[0,7.62939453125e-7]`；`(-0.5,20)` 为
  `endpoint_local_geometry_first_order`，上界 `[9.765625e-5,1.953125e-4]`；
  `(0,5)` 为同类证书，上界 `[4.8828125e-5,9.765625e-5]`。
- [x] 成本 gate 通过：hybrid unique solves `12845`（聚合按 xi 为 `4291/4287/4267`），
  dense `15384`，oracle `30744`；residual/Jacobian、fixed-rho requests 和 runner cost
  均不高于 dense，Actions wall-time 处于允许噪声范围。
- [x] evidence 已生成至
  `docs/analysis/pnjl_cep_endpoint_local_production_shadow_v4/`，包含聚合表、source/
  evidence manifest、curve index、claim ledger、plot manifest 和 9 张代表性 PNG；完整
  `curve_points.csv` 仅保留在 Actions/local artifact，外部 SHA 为
  `176f530b8f09b345ddd6d6ece40b55f1442f7c589c873885c6e924eeb7190dc8`。派生表允许的
  `NaN` 仅表示 schema 中“不适用”，未声明 NaN/Inf 均被 importer 拒绝。
- [x] 根据作者审核反馈，重新生成全部代表图：右侧面板改为独立 `rho–mu` 纵轴，增加
  Maxwell `mu_M`、共存密度和 spinodal 标记，并用线型/marker 区分三种方法；输入曲线、
  numerical source run、calculation SHA 和数值 gate 均未改变。plot policy 为
  `independent_rho_mu_zoom_with_phase_markers_v2`。
- [x] 对没有 Maxwell/spinodal/support 的 `confirmed_monotone` 切片，右侧改用原始
  production 曲线最长低斜率区间作为 display-only smooth window，并独立收紧 μ 轴；
  不参与三态分类或 gate。对弱 S 形一阶图同步减小显示 padding，以避免局部起伏被纵轴压扁。
- [ ] 当前 evidence PR 仅作为诊断候选提交，停在作者审核代表性曲线、三个 endpoint
  证书、三态一致性和成本；`full_hybrid_candidate` 不自动晋升 production/reference。
  作者审核完成前不启动 C0/C1/C2、phase-reference replay/promotion、C3/O1 或 transport。

## C0 production 输出物化失败与历史初步诊断（2026-08-07）

- [x] 新 C0 run `31103188845` 绑定 calculation SHA
  `fa390d89c836598e37b9fb75e7f8e7368283df1e`，41 个数值 shard 均在 `solve_points` 阶段
  失败；artifact 中 coarse/cascade/Stage-B 曲线已经生成，失败属于 production sweep
  summary 的 `Float64(::Nothing)` 物化错误，未形成 C0 数值 verdict。
- [x] 诊断边界收窄为 phase-result materialization contract：三个可空 geometry summary
  字段被当作必填数值转换。该修复不修改 equilibrium solver、Maxwell、geometry 容差、
  hybrid policy、旧 reference 或历史 evidence。
- [x] PR #186 已以 `main@d6cfa33edeace5aab2eb69f401ff543611a573c9` 合并，完成
  `nothing_to_inf_v1` summary 归一化、缺失字段 provenance 和 focused regression；嵌套
  convergence records 保留 `null`。
- [x] 随后 C0 run `31139870481` 仍以同一 `Float64(::Nothing)` 失败；此前将失败归因于
  sweep summary 的初步诊断不完整，不能作为当前失败位置的最终结论。
- [ ] 当前 repair 分支修复 previous-T rho-support prior 的可空过滤和 deterministic
  nearest-temperature 选择；focused CI 通过并合并后，以新的 calculation SHA 重跑 C0；
  C1/C2、reference 和 transport 继续暂停。

根因是 previous-T prior 过滤条件把缺失字段时的 `continue` 逻辑写反：缺失
`cascade_support_center/gap` 时继续对 `nothing` 执行 `Float64`，而两个字段有效时反而
丢弃了有效 prior。当前 `codex/issue-130-c0-rho-prior-contract-repair` 只修复该内部
路由并补 deterministic regression；新 calculation SHA 前不得复用上述失败 run。

## C0 accepted; cross-stage convergence route (2026-08-08)

- C0 `31149826740` is accepted as a diagnostic-only candidate at calculation
  `ffa816df0a145f73d7490db1ed9ff10c92e017a4`; it is not promoted or recomputed.
- The accepted audit distinguishes successful finite/converged final curves from internal
  geometry/interpolation records that remain unresolved. C1/C2 therefore compare stable
  physical outputs across stages rather than requiring the historical 1,650 records to be
  zero. The allowed residual class is geometry/interpolation only; solver failure, payload
  incompleteness, NaN/Inf or non-finite final curves still fail the gate.
- The telemetry repair branch adds `pnjl_phase_shard_diagnostics_v1` and
  `pnjl_phase_telemetry_aggregate_v1`. A missing historical `phase_summary.json` is recorded
  as `telemetry_unavailable` with null counters. It does not change any phase classification,
  Maxwell/geometry tolerance or rho-support policy.
- C1 (`issue130_endpoint_hybrid_c1_20260808_integral_tight`) changes only thermal quadrature
  to `rtol=1e-10, atol=1e-12`. C2 (`issue130_endpoint_hybrid_c2_20260808_grid_tight`) restores
  C0 quadrature and tightens rho/T/xi geometry/refinement contracts. C2 must also close the
  CEP bracket to 0.1 MeV and keep public anchors/first-order geometry stable before a new
  convergence evidence package is opened.

## C2 convergence comparator replay v1 (2026-08-09)

- [x] 固定 C0/C1/C2 证据：runs `31149826740`、`31235607046`、`31258823755`，calculation
  SHA `ffa816df0a145f73d7490db1ed9ff10c92e017a4`；C2 仍为 `diagnostic-only`，不重算。
- [x] comparator 改为显式区分新 xi、共享 xi 的新增/缺失 T 网格和已实际求值的
  first-order→ambiguous 分类回归；boundary 比较加入 `area_residual`，并输出 public
  anchor 三态表、geometry gate 表和 CEP endpoint-resolution/ambiguity 双宽度表。
- [x] crossover 改在双方共同 monotone/crossover 温区的物理 `mu` 并集上双向线性插值，
  不再按行号配对；物理区间由各自 CEP 的 `T_first_monotone_MeV` 约束，避免把一阶区间
  的平台误报为 crossover 风险。
- [x] solver-free evidence package 已生成到
  `docs/analysis/pnjl/c2_convergence_audit_v1/`：输入 SHA、C0/C1/C2 replay、9 个
  public classification regressions、16 个 CEP bracket failures、`xi=0.2875` crossover
  局部风险、代表图、manifest 和 claim ledger 均已保留；原始全量曲线不入仓库。
- [x] PR1 已通过 focused CI 并合并为 `main@a3dc55a9672bedba5716b18be371374515109ec1`；
  C2 comparator/evidence replay 保持 `classification_regression`，固定 calculation SHA
  `ffa816df0a145f73d7490db1ed9ff10c92e017a4`，不重算 C0/C1/C2。
- [ ] PR2 `codex/issue-130-stagec-density-certificate-feasibility` 正在实现限定范围
  Stage-C density-certificate feasibility workflow；固定 cap=12 失败时只能作为诊断记录，
  不能自动升 cap 或放宽容差。只有 Actions replay 得到 `feasible_candidate` 才能创建
  production PR。

## Stage-C density-certificate feasibility v1（2026-08-09）

- [x] solver-free postprocessor `scripts/analysis/pnjl_stagec_density_certificate_feasibility.py`
  已加入三种确定性路由：`stage_b_features_v1`、`balanced_density_features_v2` 和
  `geometry_feedback_v2`；路由只读取完整 Stage-B 曲线及当前请求内诊断，oracle 标签只在
  事后 gate 使用。
- [x] Actions numerical runner `scripts/analysis/pnjl_stagec_density_certificate_job.jl`
  已固定 9 个密度锚点和 6 个 first-order/monotone controls，并显式记录 calculation SHA、
  postprocess SHA、曲线哈希、cache/solver telemetry 和 `reference_write=false`。
- [x] workflow `.github/workflows/pnjl-stagec-density-certificate-feasibility.yml` 已加入
  10 个 xi × 3 个 method 的 30-job matrix，以及 `aggregate_replay`；数值 job 从 immutable
  calculation worktree 运行，collector 从 workflow head 运行。
- [x] focused Python/replay/governance checks 已通过：新合同 `7 passed`，相关 shadow/
  endpoint/Maxwell replay `31 passed`，Python/Julia syntax、script/model/solver/data-path/
  PNJL migration、docs 和 active-docs governance 均通过。Actions 仍只消费固定 calculation
  SHA。
- [x] PR #190 已创建并合并为 `main@2edfd8c8760f1c03d65c82c520a356b8d6c5c3ba`；
  feasibility 未通过前不实现 production ranking policy，不启动新的 C0/C1/C2、reference 或 transport。

## Stage-C density feasibility numerical run 31295557583（2026-08-09）

- [x] 以 production merge SHA `2edfd8c8760f1c03d65c82c520a356b8d6c5c3ba` 作为 workflow
  head、以固定 calculation SHA `ffa816df0a145f73d7490db1ed9ff10c92e017a4` 启动 30-job
  numerical matrix；前置 checkout、双 SHA provenance 和依赖安装均通过。
- [x] run 在首批 independent-oracle job 的 runner 输出阶段暴露统一合同错误：摘要使用了
  未定义的 `STAGE_C_STEP`，不是 equilibrium solver、Maxwell 或曲线数值失败。已取消该
  必然继续失败的 run，原始失败日志保留在 Actions。
- [x] 窄范围 runner-contract repair 已完成；只用新的 postprocess SHA 重跑 failed jobs，
  未重算 calculation SHA。
- [x] PR #191 已合并为 `main@bdc580b78916b4428c980caf76cc74f101573854`；run `31295557583`
  的失败 job keys 已保留为 failed-only rerun 输入，workflow 现在支持显式
  `rerun_failed_only=true` 与 `failed_job_keys`，并在该模式跳过 aggregate，避免用不完整
  artifact 伪造全矩阵 verdict。
- [x] PR #192 已合并为 `main@4db46e1ef1694b28171ff79d6c00700a507b35ce`；run `31296088809`
  的 10 个 failed-only artifact 全部成功，双 SHA、CSV hash、finite/converged、重复键、成本
  守恒和 cap=12 审计通过。随后完整 run `31296511813` 的 30/30 numerical jobs 成功。

## Stage-C density feasibility replay 收口（2026-08-09）

- [x] 完整 run `31296511813` 的第一次 aggregate 只因 plotter 使用非法 MathText 标签失败；数值
  artifact 未失败。PR #193 修复标签、增加 plot smoke test，并允许“source 数值矩阵完成但 aggregate
  总结论为 failure”的 artifact-only replay；focused CI 7 项全绿后 squash merge 为
  `main@11fbf9a71575ea207ed3233608c2bc2ade1ee346`。
- [x] main head 的 aggregate replay `31300285869` 成功，固定 calculation SHA
  `ffa816df0a145f73d7490db1ed9ff10c92e017a4`，source postprocess SHA 为 `4db46e1e…`，
  replay/producer SHA 为 `11fbf9a7…`，`source_job_count=30`、`solver_called=false`，图像和
  150 个输入文件 hash 均已记录。
- [x] 最终 verdict 为 `integration_failed`：三条路由在 cap `12/16/24` 都有 12 个
  classification mismatch；geometry、finite/converged、candidate uniqueness 和 cost gate
  通过，但 7 个 CEP bracket（均超过 `0.1 MeV`）与 `xi=0.2875` 的 3 个 crossover risk 未通过。
- [ ] 因未得到 `feasible_candidate`，保留 replay evidence 为 diagnostic-only；不创建 Stage-C
  production PR，不启动新的 C0/C1/C2、reference promotion、C3/O1 或 transport。后续需先由作者
  判断 classification/CEP/crossover 阻塞的处理路线。
## Stage-C feasibility contract v2（2026-08-09）

- [x] 从本地 `main@88336a46` 创建分支
  `codex/issue-130-stagec-feasibility-contract-v2`；保留 v1 runner/workflow/evidence 不变，
  不向本地 `main` 推送。
- [x] 新增 solver-free production-parity evaluator
  `scripts/analysis/pnjl_stagec_density_certificate_feasibility_v2.jl`：固定 source run
  `31296511813`、calculation SHA `ffa816df…`，直接调用 `Models` 的 S-shape/Maxwell/geometry
  合同，路由不读取 oracle 标签，且明确 `solver_called=false`。
- [x] evaluator 按 numerical job 的 `rho_level` 分离 Stage A（production hybrid
  `0.05→0.025`）和同一 independent-oracle session 的 Stage B（全域
  `0.0125→0.00625`, `level=0`）与不重合的 `0.003125` Stage-C pool（`level=1`）；
  memoized dense 仅作为成本基线。完整 Stage-B 曲线与实际补点逐点合并重算。这修正了
  v1 replay 将 oracle fine 点提前消费、把 Stage-C additions 误记为零的风险。
- [x] 三条 route 的候选序列、逐项 geometry component、完整 Stage-A/B/实际 Stage-C key 并集
  成本和 cap frontier 均写入 versioned 输出；manifest 最后写入，覆盖 README/AUDIT/plot
  metadata 和所有输入/派生文件 hash。
- [x] 本地外部 artifact replay 已完成：`feasible_candidate`、cap `12`、route
  `stage_b_features_v1`，hybrid unique solves `7985`，dense baseline `9615`，15 个锚点
  classification/geometry/finite/cost gates 通过，`solver_called=false`。该结果只用于开发
  自检，不是 Actions evidence 或 production 授权。
- [x] 新增固定 matplotlib plotter、aggregate-only workflow
  `.github/workflows/pnjl-stagec-density-certificate-feasibility-v2.yml` 和 focused
  Julia/Python contract tests；workflow 从固定 calculation worktree 运行 evaluator，并校验
  source run 的 30 个 numerical jobs、双 SHA 与 artifact hashes。
- [ ] PR focused CI 通过后，运行 v2 aggregate replay；无论 verdict 如何导入独立不可变 v2
  evidence。只有 Actions replay 再次得到 `feasible_candidate` 才创建 production PR；此前
  不运行新的 C0/C1/C2、CEP/crossover numerical workflow、reference、C3/O1 或 transport。

## Stage-C feasibility contract v2 Actions 收口与 production integration（2026-08-09）

- [x] PR #194 已 squash merge，workflow head 为
  `72c49c26a061cac7b71e66543ebdfd2c91ec3fca`；aggregate replay run `31308813742`
  成功完成 30 个 source artifact 校验、Julia evaluator、绘图和 evidence contract。
- [x] 正式 verdict 为 `feasible_candidate`：selected policy 为
  `stage_b_features_v1`、cap `12`，hybrid unique solves `7985`，dense baseline `9615`；
  15 个 anchors 无 classification/geometry/finite/candidate/cost failure，且
  `solver_called=false`。source run `31296511813` 的原始 `failure` 结论只作为 provenance，
  不改变 replay 的 solver-free 证据。
- [x] aggregate evidence 已下载并保存在仓库外
  `D:\Desktop\Julia_RelaxTime_issue130_artifacts\stagec_feasibility_v2_actions_31308813742`，
  保持 v1 evidence 不变。
- [ ] 从 merge SHA 创建 `codex/issue-130-stagec-density-production-v2`，原样落地
  selected policy；production focused CI、ready PR、targeted/CEP/crossover/full shadow
  仍未运行，且在 `full_hybrid_candidate` 与作者审核前不启动新的 C0/C1/C2、reference、
  C3/O1 或 transport。
