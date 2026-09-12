# RS `publication_clean_v2` figure-layer snapshot

本目录是 `publication_clean_v2` 语义修正版的 figure-layer provenance 入口；实际 PNG 位于
`data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v2/`。

## 语义修正

- 一阶相变端点之间不画连续桥接线；断线区间和端点相别见分析包的
  `tables/boundary_gap_map.csv`。
- 历史 CEP/phase-switch 中点只保留为 provenance，全部 `render_marker=false`，图例不再标为 CEP。
- 端点 marker 使用当前 raw `phase_curr`：quark 为手征恢复相，hadron 为手征破缺相。
- mode-B `T=120 MeV, μB=900 MeV` 使用当前 raw 的 phase_curr 切换端点 `[-0.13,-0.12]`；旧
  `[-0.14,-0.13]` 只作为历史 phase-kind bracket 保留。
- mode-B 相变线左侧高 `η/s` 可靠度较低是写作层 caveat，不改变 figure 数值或筛点。
- 存在已审计一阶断线且正值动态范围达到 100 倍的图使用 log-y；其余图保持线性坐标。该变换只改善显示，不改变数据或相变语义。
- 图例使用 `α_T`/`μ_B`，MeV 显示值四舍五入到个位；原始 CSV、目录键和 provenance 精度不变。

## Provenance 与边界

- 源分析包：`docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2/`；
- source case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2`；
- figure-layer manifest：`figure_layer_manifest.json`；公共图层 manifest：
  `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v2/plot_manifest.json`；
- 本次只做 byte-preserving PNG mirror；没有调用 solver，也没有修改 raw result、registry 或 v1 快照；
- `manuscript_eligible=false`，v2 仍需作者审核后再决定是否作为论文公开候选。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v2.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v2.py
```
