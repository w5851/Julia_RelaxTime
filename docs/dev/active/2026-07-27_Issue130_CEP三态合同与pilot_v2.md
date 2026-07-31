# Issue #130：CEP 三态合同与窄窗口 pilot v2

创建日期：2026-07-27

状态：三态合同 PR #147 已 squash merge；pilot v2 数值 Actions、aggregate replay 和
evidence 导入已完成。hybrid shadow replay 与 Stage-C offline feasibility 已完成，
当前因 `maxwell_candidate_inconclusive` 停在作者物理审核。PR #146 已 squash merge（merge SHA
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
- [x] 当前 PR2 预期 verdict 为 `maxwell_candidate_inconclusive`：所有 cap 的成本
  估计均低于 dense 的 15384 unique solves，但有两个弱 S/CEP 邻域 first-order
  oracle anchors（`xi=-0.5,T=147.0947265625` 与 `xi=0.5,T=106.9599609375`）
  无法形成唯一稳定 topology；没有选出 `selected_policy`。
- [ ] 因 feasibility 未通过，不创建 production adaptive Stage-C PR，不启动
  targeted/full shadow、C0/C1/C2、reference replay 或 transport；等待作者判断弱 S
  与 Maxwell 候选的后续物理/数值方案。
