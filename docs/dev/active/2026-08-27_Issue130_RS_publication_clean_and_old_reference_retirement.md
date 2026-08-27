# Issue #130：RS `publication_clean` 派生层与 old-reference retirement

状态：active；publication-clean v1 已生成，等待作者审核；old-reference retirement 尚未开始。

## 范围

本任务承接 Issue #130 RS transport 的两步后处理治理：

1. 从已登记 `approved`、但仍为 `manuscript_eligible=false` 的两棵 `prod_v2` raw result tree 生成可审阅的 `publication_clean_v1` 派生层；
2. 在派生层审核完成后，单独审计并决定旧 RS `prod_v1` 的 retirement。

本任务不重新调用 equilibrium/transport solver，不修改 raw CSV、canonical figures、production registry 或 phase-reference legacy fallback。publication-clean 的插值只改变派生 `clean_value`，不构成新的数值结果或收敛证书。

## 固定输入

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow head：`22874505877491754eed27519ad8a7b871c82571`
- mode-A：`data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`
- mode-B：`data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`
- 已审核配方：`docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_v2_pole_sensitive_rendering/tables/paper_display_replacements.csv` 与 `paper_first_order_markers.csv`

## publication-clean v1 验收

- [x] 读取两棵 `prod_v2` raw tree，校验 scan/diagnostic 行数、重复键、非有限值、failed points、双 SHA 和 source/derived solver provenance。
- [x] 19 个非一阶显示替换键按当前 raw 左右邻点重算；保留旧配方值用于迁移审计，不能直接复制旧 `prod_v1` 数值。
- [x] 旧配方 `mode_a, μB=900, αT=1.0, ξ=0` 在当前 direct-coexistence grid 中缺失时，保留 `ξ=−0.003/+0.003` 两侧 raw 端点，并以其确认区间中点 `ξ=0` 作为 display-only 星标；不生成唯一 `ξ=0` raw 输运值。
- [x] 生成 `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1/`：长表、替换表、marker 对齐表、曲线索引、18 张图、manifest、plot manifest 和 claim ledger。
- [x] 派生层所有行写入 `canonical_data_modified=false`；manifest 明确 `solver_called=false`、`manuscript_eligible=false`。
- [ ] 作者审核 18 张同构图、19 个继承显示替换值、两个有限区间中点 marker、两侧端点 provenance 和已知 `xi=-0.01` bulk 分支排除项。
- [ ] 作者审核本轮新增的 4 个 T=200 局部平滑候选：mode-B `(μB=900, ξ=-0.10)` 与 `(μB=0, ξ=0.36)` 的 `η/s`、`ζ/s`；候选只写入 clean layer。
- [x] 修正 T=120、μB=900 的历史星标：`xi=-0.09` 只保留在审计记录；依据当前 `phase_reference_kind: crossover → first_order` 的作者确认区间 `[-0.14,-0.13]`，publication-clean 改渲染中点 `xi=-0.135`。严格 CEP 表仍独立记录固定切片坐标审计。
- [x] 回溯 mode-A `μB=450, αT=1.0, ξ=-0.20`：当前点为 crossover/no-transition，既有 mechanism window 支持 `simple_1m4KΠ` 小分母归因；不修改 raw。
- [x] 统一论文 marker 语义：严格 CEP/相标签切换本质上保留有限宽度 bracket；`tables/publication_marker_map.csv` 记录两侧 raw 端点、phase label、质量标志及中点 display 值，图中只绘制中点，不把中点写入 raw/clean 长表。

## old-reference retirement 后续门禁

作者完成 publication-clean 审核后，另立 retirement audit/PR，至少完成：

- [ ] 枚举所有代码、配置、脚本、文档、workflow 和 figure manifest 对旧 `prod_v1` 路径的引用；区分 active consumer、历史证据和 rollback/fallback。
- [ ] 在不调用 solver 的前提下验证默认 `prod_v2` 分析路径、显式 `--phase-reference-mode legacy` 回退和 versioned legacy snapshot 均可解析。
- [ ] 固化旧 mode-A/mode-B `prod_v1` tree hash、引用清单、fallback smoke、删除/保留 allowlist 与回退步骤。
- [ ] 未经作者单独授权不物理删除旧 `prod_v1`；若仍有 active consumer 或 manuscript/复现依赖，保持 versioned snapshot。

## 证据与复现

- 生成器：`scripts/analysis/relaxtime/build_phase_guided_publication_clean_v1.py`
- Python tests：`tests/unit/python/test_phase_guided_publication_clean_v1.py`
- 派生包：`docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1/`
- 复现命令：

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v1.py
python -m pytest -q tests/unit/python/test_phase_guided_publication_clean_v1.py
```

## 停止边界

publication-clean v1 若校验失败、输入 hash 漂移、替换键缺失或图层无法重建，则停止并保留 raw。old-reference retirement 不得借此自动完成；任何删除、registry retirement 或 manuscript eligibility 晋升都需要新的作者审核与独立 PR。
