# Issue #130：CEP 三态合同与窄窗口 pilot v2

创建日期：2026-07-27

状态：三态合同已实现并完成本地 focused 验证，已提交为独立 ready PR #147，等待作者审核；基于
`main@5071f6c80bb10c812358b855e8da2bde8a758f9d` 的分支
`codex/issue-130-cep-three-state-contract`。PR #146 已 squash merge（merge SHA
`5071f6c80bb10c812358b855e8da2bde8a758f9d`；最终 Actions run
`30199406478`，9/9 matrix + aggregate success），v1 evidence 结论保持
`diagnostic_only`，不晋升 reference。

三态合同 PR：<https://github.com/w5851/Julia_RelaxTime/pull/147>，head
`8f3b1110`。

## 顺序与边界

1. [x] 合并并清理 PR #146，保留 v1 evidence 与 diagnostic-only 结论。
2. [x] 完成三态 CEP 合同实现与 focused 验证；ready PR #147 已创建，等待作者审核再合并。
3. [ ] 三态 PR 合并后，从 merge SHA 创建 `codex/issue-130-cep-narrow-pilot-v2`，只经
   GitHub Actions 运行 cascade discovery、dense baseline 和独立 oracle。
4. [ ] pilot v2 停在作者对 CEP 区间、原始 rho–mu 曲线、准确度和性能证据的物理审核。
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
