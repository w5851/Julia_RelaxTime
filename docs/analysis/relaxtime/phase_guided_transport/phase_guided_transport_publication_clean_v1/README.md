# Issue #130 RS `publication_clean_v1` 派生层

## 目的与边界

本包从作者接受的两套 `prod_v2` 原始 RS transport 结果生成论文显示候选层。它只做 solver-free 后处理：不修改 `data/outputs/results/**`、`production_registry.json` 或正式 figure，不调用 equilibrium/transport solver，也不把派生值写回 raw CSV。

`publication_clean` 的含义是：沿用已审核的极点/小分母显示配方，并在当前 `prod_v2` 网格上重新取左右真实邻点计算显示插值；一阶相变点继续保留 raw 值。本轮另记录作者提出的局部平滑候选，但它们仍待本轮作者确认。该层默认 `manuscript_eligible=false`，须经作者审核后才能用于论文。

## 输入 provenance

| mode | scan rows | diagnostic rows | scan SHA256 | diagnostic SHA256 |
| --- | ---: | ---: | --- | --- |
| mode_a | 910 | 38220 | `144f59dd6d74b1d13b24b793328e620e796fb645b7af9418174ba5a55e360587` | `6eaf02101e330762569d028fe1398c94184e243ea6d7e39a8db69718a6e4e246` |
| mode_b | 909 | 38178 | `3b73b937f1f09cf7112b3a42219088f98b852e5ef0354174b66f71823f7bb29c` | `ada5d7c204a02e705cfbb5744bedbcc68158d54cce55faa0c8b92d0ba50d0d0f` |

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow head：`22874505877491754eed27519ad8a7b871c82571`
- source result case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2`
- source solver：已调用；本次派生：`solver_called=false`
- marker contract：`confirmed_interval_midpoint_v1`（严格 CEP/相标签切换区间的中点 display proxy）
- 旧显示配方：`docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_v2_pole_sensitive_rendering/tables/paper_display_replacements.csv`（SHA 记录于 manifest）
- 本包生成图：18 张 PNG，见 `figures/plot_manifest.json`
- 历史 marker：9 条审计记录，其中 0 条历史 raw marker 实际渲染；T=120、μB=900 的 3 条历史星标被抑制。
- CEP-slice 审计：1 个固定切片，严格固定切片交点星标为 0 条。
- 论文显示 marker：6 条（3 个 observable × 2 个已确认区间），均为区间中点的 display-only proxy。

## 派生规则

1. 继承配方中的 19 个非一阶显示替换点从当前 raw 的 `left_xi < target_xi < right_xi` 三点线性插值；`raw_value` 与左右端点同时保留。
2. 旧配方中 `mode_a, μB=900, αT=1.0, ξ=0` 当前不存在：direct-coexistence 合同保留 `ξ=−0.003/+0.003` 两侧 raw 端点，并在图上只显示其区间中点 `ξ=0` 的 display-only marker；不把中点写成输运 raw 值。
3. 一阶点不做平滑；所有派生行带 `canonical_data_modified=false`。`quality_flag`/`quality_reason` 原样带入。
4. 已知 `mode_a, μB=900, αT=1.0, ξ=−0.01` bulk 分支问题仍列为排除项；本包不会用插值掩盖它。
5. 本轮作者审阅新增 4 个 T=200 MeV 局部平滑候选：`(mode B, μB=900, ξ=−0.10)` 与 `(mode B, μB=0, ξ=0.36)` 的 `η/s`、`ζ/s`；全部按当前 raw 左右邻点线性重算，单列于 `review_adjustment_map.csv`，不改变原始结果。
6. 历史配方中 T=120 MeV、μB=900 MeV、ξ=−0.09 的星标只作为审计输入，不进入 publication-clean 图；图上改为显示作者确认的 `phase_reference_kind: crossover → first_order` 最小切换区间 `ξ∈[−0.14,−0.13]` 的中点 `ξ=−0.135`。严格 CEP 表仍作为独立坐标审计，区间中点是论文显示 proxy，不宣称零宽度精确 CEP 坐标。
7. `publication_marker_map.csv` 保存两个已确认区间的左右端点 phase label、质量标志和 observable 值；`ξ=0` 与 `ξ=−0.135` 的星标都不进入 `publication_clean_points.csv`。
8. mode-A μB=450 MeV、αT=1.0、ξ=−0.20 的非一阶斜率变化与既有 `simple_1m4KΠ` 小分母机制窗口一致；这属于机制归因证据，不是新的相变标签。

## 结果文件

- `tables/input_inventory.csv`：输入路径、行数、哈希及 solver provenance。
- `tables/replacement_map.csv`：19 个继承配方加 4 个本轮作者审阅候选的当前邻点插值映射。
- `tables/review_adjustment_map.csv`：4 个本轮作者请求的 T=200 局部平滑候选及局部残差。
- `tables/first_order_marker_map.csv`：旧标记与当前 raw/direct-coexistence 端点合同的逐项对齐，并记录是否渲染。
- `tables/marker_semantics_audit.csv`：历史 marker 语义和 publication-clean 渲染决定。
- `tables/cep_marker_audit.csv`：固定输运切片与 strict CEP 表的坐标配对；严格交点与论文区间中点是分开的语义。
- `tables/cep_marker_map.csv`：仅在存在 strict 固定切片交点时生成的 CEP 星标及其显示值；当前为空。
- `tables/publication_marker_map.csv`：两个作者确认的最小相标签切换区间中点 marker，包含左右端点 provenance 和 display 值。
- `tables/publication_clean_points.csv`：三种论文展示 observable 的长表，含 raw/clean/status。
- `tables/curve_index.csv`：18 条 panel/series/observable 曲线的覆盖和替换计数。
- `tables/claim_ledger.csv`：证据强度、范围限制和未声明事项。
- `figures/plot_manifest.json` 与 `figures/**.png`：同构的 18 张显示候选图，未声明为 manuscript-ready。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v1.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v1.py
```

执行脚本只写入本目录；任何 raw、registry、旧 `prod_v1` 和 legacy fallback 的变更都属于越界操作。
