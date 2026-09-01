# RS `publication_clean_v1` figure-layer snapshot

本目录是 Issue #130 RS `publication_clean_v1` 派生显示层的 figure-layer 说明与 provenance
入口；对应的 mode_a/mode_b 图像位于 `data/outputs/figures/` 下的实际目录中。

## 目录与内容

- `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/mode_a_fixed_muB_phase_scaled/`：3 个 `muB` panel × 3 个 observable，共 9 张 PNG；
- `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/mode_b_fixed_T_sparse_muB/`：3 个 `T` panel × 3 个 observable，共 9 张 PNG；
- `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1/plot_manifest.json`：沿用 PR275 的图层 manifest，仅将图像路径改为 figure-layer 目录；
- `figure_layer_manifest.json`：源包、计算/workflow SHA、逐文件 SHA-256 和字节保持记录。

## Provenance 与边界

- 派生源包：`docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1/`；
- source case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2`；
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；
- workflow head：`22874505877491754eed27519ad8a7b871c82571`；
- 本次仅复制已审核的派生 PNG，内容按字节保持；没有调用 solver，也没有修改 raw result、
  `production_registry.json` 或分析包中的表格/图像；
- `publication_clean_v1` 仍是 display-only 派生层，`manuscript_eligible=false`，不等同于
  raw production/reference 晋升。

分析包中的原始 18 张图仍保留为可追溯证据；后续如需改变图像语义，应创建新的版本目录，
不要覆盖本快照。
