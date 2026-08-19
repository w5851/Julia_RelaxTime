# Issue #130：Maxwell CEP 近端局部补点与派生补全任务单

状态：active；solver-free preflight 已完成，作者已授权 11-target numerical pilot；workflow/collector 正在独立分支 focused 验证。
父任务为 `issue130-phase`。这是与 crossover μ endpoint refinement 分开的
required follow-up，只处理 Maxwell 侧在 CEP 附近的 support/geometry 缺口。

## 目标

在固定 C2 v5 证据上，先用 solver-free 方式建立 Maxwell 近 CEP 的 target list
和成本 envelope；只有预检通过并获得单独授权后，才对少量代表性 `xi` 运行
numerical pilot。目标是判断缺失来自 rho geometry/refinement 还是端点采样，
不是把 unresolved 记录直接改成通过。

## 固定输入与预检结果

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- v5 证据：`docs/analysis/pnjl/c2_phase_surfaces_diagnostic_v5_no_triangulation/`。
- 目标表：`maxwell_surface_point_status.csv`、`grid_unresolved_diagnostics.csv`、
  `crossover_maxwell_endpoint_separation.csv`。
- preflight 输出：`docs/analysis/pnjl/issue130_endpoint_refinement_preflight_v1/maxwell_local/`。
  窗口为 `T <= T_CEP_bracket_low` 且距离不超过 `8 MeV`；共 760 个候选，484 个
  已有 boundary 行，276 个为 `input_incomplete`，代表性 xi 共 11 个 pilot 目标。

## 预检与 numerical pilot 合同

- 预检按固定 `(xi,T)` 识别 CEP 下方 unresolved Maxwell 行；缺少 boundary 行的
  目标只保留为 `input_incomplete`，不得自动进入 numerical。
- 每个数值目标必须重新取得完整 `rho-mu` 曲线，执行公共 Maxwell、三态和 geometry；
  candidate 必须唯一，三个交点 finite/有序，且不放宽 position/density/area 门禁。
- 成本只报告 `rho=0:0.00625:4` 的 641 点单位；wall time 未测量。当前 pilot 为
  `authorized_pending_workflow_merge`，预检通过不等于 Maxwell certificate 或
  production feasibility。pilot 额外记录最多 12 个当前 crossing bracket midpoint，
  并按每级 strict candidate 结果计费。

## Numerical pilot implementation

- workflow：`.github/workflows/pnjl-issue130-maxwell-cep-local-pilot.yml`；
  runner：`scripts/analysis/pnjl_issue130_maxwell_cep_local_pilot.jl`；
  collector：`scripts/analysis/collect_pnjl_issue130_maxwell_cep_local_pilot.py`。
- 固定 11 个 `pilot_candidate` target；每个目标先完整计算
  `rho=0:0.00625:4`，随后最多补 12 个当前三交点 bracket midpoint；不读取 oracle
  标签，不写 phase-reference。
- focused CI 通过并合并后，使用 calculation SHA
  `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48` 触发 numerical pilot；若有 solver、
  finite、重复 key、provenance 或 candidate/geometry failure，停止并保留
  `pilot_inconclusive` evidence。

## Derived completion 边界

`certified_layer` 保留原始 C2 证书；`completed_layer` 只能在同一 Maxwell 物理区、
相邻已证实行之间做显式 `interpolated_noncertified` 派生，不能伪造 candidate、
`maxwell_area` 或 strict geometry certificate，也不能自动晋升 phase-reference。

不修改 equilibrium solver、Maxwell、三态规则、endpoint policy 或容差；不重跑 C0/C1/C2。
