# Issue #130：versioned phase-reference import 任务单

状态：review；promotion gate v1 已在 PR #250 合并，新的 import candidate 已在本分支生成，
等待独立 PR 审核。该任务只导入一个新的 canonical sibling，不切换默认 runtime reference。

## 固定输入

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- promotion gate merge SHA：`10231d83f8ed56b684f6acbd74a71ce85fdb47cf`。
- source numerical run：`32354095831`；solver-free replay：`32451053476`。
- source package：`docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/`。
- gate package：`docs/analysis/pnjl/phase_reference/issue130_phase_reference_promotion_gate_v1/`。

## Import 范围

导入脚本 `scripts/analysis/pnjl/import_issue130_phase_reference.py` 将已审核的三层证据
逐字节复制到：

- `data/reference/pnjl/issue130_phase_reference_v1/strict/`；
- `data/reference/pnjl/issue130_phase_reference_v1/derived/`；
- `data/reference/pnjl/issue130_phase_reference_v1/render/`；
- `data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v1/`。

新 sibling 的 `manifest.json` 保存 source/imported hash、旧 reference 快照、统一 xi 网格、
三层语义和 `runtime_consumption=false`。旧的 `boundary.csv`、`cep.csv`、`spinodals.csv`、
`crossover_dense.csv` 与 `phase_reference_dense_manifest.json` 不被覆盖。

## 验收与边界

- strict 层保留 native unresolved geometry；derived 层的非原生行仍为
  `interpolated_noncertified`；render 层不三角化缺失单元。
- import 阶段不调用 solver，不重新计算 Maxwell/CEP，不修改 tolerance 或 production policy。
- import 阶段不切换 runtime consumer；`data/reference/pnjl/issue130_phase_reference_v1`
  先作为 versioned candidate，后续另开 runtime-consumer 评估任务。
- focused Python tests、CSV schema、旧 reference SHA、manifest/plot hash 和仓库治理检查
  必须通过后才可合并。

## 后续

1. 审核并合并本 import PR。
2. 在新 sibling 上做只读 runtime consumer compatibility audit；确认是否需要显式配置选择。
3. runtime audit 通过并经作者授权后，单独评估 RS transport；不得因 import 成功自动启动 transport。
