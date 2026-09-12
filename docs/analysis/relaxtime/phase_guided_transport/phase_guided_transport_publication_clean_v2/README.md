# Issue #130 RS `publication_clean_v2` 语义修正版

## 目的与边界

本包是 `publication_clean_v1` 的版本化、solver-free 显示语义修正。v1 保持不变；v2 不修改 `data/outputs/results/**`、production registry 或任何 raw CSV，也不调用 equilibrium/transport solver。

本轮只处理两件事：

1. 一阶相变端点之间不再用实线桥接；端点使用当前 raw `phase_curr` 标注为手征恢复相（quark）或手征破缺相（hadron）。
2. 旧的 CEP/phase-switch 中点星标不再渲染；中点和旧 bracket 仅保留在 provenance 表中，图例不再出现 `CEP`。

`mode_b, T=120 MeV, mu_B=900 MeV` 的历史 phase-kind bracket 为 `[-0.14,-0.13]`，但当前 raw 扫描的 `phase_curr` 实际由 quark 切换为 hadron 的相邻端点为 `[-0.13,-0.12]`。v2 用后者断线，并在 `boundary_gap_map.csv` 保留前者作为来源 bracket；这是语义对齐，不是新的 solver 复核。

## 输入 provenance

| mode | scan rows | diagnostic rows | scan SHA256 | diagnostic SHA256 |
| --- | ---: | ---: | --- | --- |
| mode_a | 910 | 38220 | `144f59dd6d74b1d13b24b793328e620e796fb645b7af9418174ba5a55e360587` | `6eaf02101e330762569d028fe1398c94184e243ea6d7e39a8db69718a6e4e246` |
| mode_b | 909 | 38178 | `3b73b937f1f09cf7112b3a42219088f98b852e5ef0354174b66f71823f7bb29c` | `ada5d7c204a02e705cfbb5744bedbcc68158d54cce55faa0c8b92d0ba50d0d0f` |

- source case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2`；calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；workflow head：`22874505877491754eed27519ad8a7b871c82571`。
- source solver 已调用；本次派生 `solver_called=false`。
- v1 快照仍保留；v2 另建目录，避免破坏既有 manifest/hash。
- 本包生成图：72 张 PNG（6 个 panel × 12 个 observable）；曲线索引 216 条。
- 本轮平滑候选：4 条；仍是 display-only，raw 值和现有 provenance 不变。

## 断线与端点合同

| boundary | gap | endpoint semantics |
| --- | --- | --- |
| mode-A μB=900, αT=1 | `[-0.003,+0.003]` | quark → 手征恢复相；hadron → 手征破缺相 |
| mode-B T=120, μB=900 | `[-0.13,-0.12]` | quark → 手征恢复相；hadron → 手征破缺相 |

`publication_marker_map.csv` 中的历史 midpoint 行全部 `render_marker=false`；它们不再是 CEP 图形标记。`boundary_gap_map.csv` 是绘图实际使用的断线/端点证据表。
`phase_switch_inventory.csv` 列出 raw `phase_curr` 的全部相邻切换，并区分 crossover/continuation bookkeeping 与实际纳入 v2 断线的两处一阶端点。

## 写作层 caveat

mode-B 相变线左侧的高 `eta/s` 区域可靠度较低。这一条是论文分析/写作时的解释边界：v2 不把它当作自动误差条，不自动删点，也不修改数值。若论文需要定量误差声明，应另立诊断任务。

## 未做事项

- 本轮不做独立 bulk 全局分支复核；历史 `bulk_derivative_branch_audit.csv` 继续作为历史证据。
- 不把断线或端点标记写回 raw/reference；不生成新的 CEP 数值。
- `manuscript_eligible=false`，待作者审核 v2 图后再决定是否作为公开候选。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v2.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v2.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v2.py
```
