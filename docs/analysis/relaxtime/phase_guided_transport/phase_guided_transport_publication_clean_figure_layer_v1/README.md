# RS `publication_clean_v1` figure-layer snapshot

本目录是 Issue #130 RS `publication_clean_v1` 派生显示层的 figure-layer 说明与 provenance
入口；对应的 mode_a/mode_b 图像位于 `data/outputs/figures/` 下的实际目录中。

## 目录与内容

- `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/mode_a_fixed_muB_phase_scaled/`：3 个 `muB` panel × 12 个 observable，共 36 张 PNG；
- `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/mode_b_fixed_T_sparse_muB/`：3 个 `T` panel × 12 个 observable，共 36 张 PNG；
- `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/plot_manifest.json`：由分析层 manifest 经 byte-preserving sync 生成，路径映射到 figure-layer 目录；
- `figure_layer_manifest.json`：源包、计算/workflow SHA、逐文件 SHA-256 和字节保持记录。

## Provenance 与边界

- 派生源包：`docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1/`；
- source case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2`；
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；
- workflow head：`22874505877491754eed27519ad8a7b871c82571`；
- 本次仅复制已审核的派生 PNG，内容按字节保持；没有调用 solver，也没有修改 raw result、
  `production_registry.json` 或分析包中的表格/图像；
- 新增的 `tau_*`、`eta`、`sigma`、`zeta` 图由当前 `prod_v2` raw 值直接绘制；tau-specific 清理配方尚未建立，
  因而不能将这些图表述为 tau 数值清理或新的收敛证明；字段政策见分析包的 `tables/observable_policy.csv`；
- `publication_clean_v1` 仍是 display-only 派生层，`manuscript_eligible=false`，不等同于
  raw production/reference 晋升。

分析包中的原始 72 张图仍保留为可追溯证据；后续如需改变图像语义，应创建新的版本目录，
 不要覆盖本快照。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v1.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer.py
```

第二条命令只同步 PNG 和 figure-layer manifests，并验证源图与目标图的 SHA-256 一致；不会写入 raw result。
