# `issue130_phase_reference_v1` 去重说明（2026-11-27）

本目录保留为**数据来源记录**，其表数据已去重到 `issue130_phase_reference_v2`。

## 为什么可以删表

v1 与 v2 是**同一次计算**的两种封装：

| 项 | v1 | v2 |
| --- | --- | --- |
| `calculation_sha` | `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48` | 同 |
| `source_run_id` | 32354095831 | 同 |
| 生成时间 | 2026-08-21 | 2026-08-29 |
| 状态 | `imported_candidate` | `accepted_for_downstream` |

四张面（cep 161 行、crossover 3 135 行、maxwell 12 537 行、spinodal 11 989 行）
逐格数值最大差 **0.000e+00**。

## 删除了什么

| 对象 | 体积 | 判据 |
| --- | --- | --- |
| `render/tables/`（5 文件） | 3.90 MB | 与 `issue130_phase_reference_v2/render/tables/` **字节相同**（sha256 逐一核对） |
| `strict/tables/`（4 文件） | 3.60 MB | 与 `issue130_phase_reference_v2/strict/tables/` **字节相同** |
| `data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v1/` | 14.75 MB | 2026-08-22 的图审阅件，已退役 |
| `scripts/analysis/pnjl/render_issue130_phase_surface_v{2,3,9}.py` | 70 KB | 只服务上述图件 |
| `tests/unit/python/test_pnjl_phase_surface_v{2,3,9}.py` | — | 只覆盖上述渲染链 |

## 保留了什么（不可删）

| 对象 | 原因 |
| --- | --- |
| `manifest.json` / `README.md` / `IMPORT_AUDIT.md` | v2 的 `source_manifest` 按 sha256 `7d37aedd8ec7a1c4…` 指向本目录 `manifest.json` |
| `derived/manifest.json`、`derived/README.md` | v2 的 `internal_build_inputs` 指向 `v1/derived`；该 manifest 现已记录被删表的 sha256 与行数 |
| `render/manifest.json`、`strict/manifest.json` 及其 README | 记录被去重表的历史哈希与用途；台账 `task_tracks.toml` 直接引用 `render/manifest.json` |

## 如果需要复核

`derived/tables/` 的数值可由 v2 重建：

```python
# v2 accepted 表与 v1 derived 表逐格相同，列名不同
'cep_boundary_derived_reference_v1.csv'  <-> 'cep_boundary_accepted_phase_map_v1.csv'
'crossover_surface_derived_reference_v1.csv' <-> 'crossover_surface_accepted_phase_map_v1.csv'
'maxwell_surface_derived_reference_v1.csv' <-> 'maxwell_surface_accepted_phase_map_v1.csv'
'spinodal_surface_derived_reference_v1.csv' <-> 'spinodal_surface_accepted_phase_map_v1.csv'
```

## 第二轮：`derived/tables/`（2026-11-27）

第一轮保留 `derived/tables/` 的理由是"该封装的唯一副本"。第二轮删除，因为那个唯一性已失去意义：

- 它唯一的消费者是 `render_issue130_phase_surface_v{2,3,9}.py` 渲染链，该链已在第一轮退役；
- 5 张表与 `issue130_phase_reference_v2/accepted/tables/` **逐格数值相同**（四张面最大差 `0.000e+00`），
  仅列名与格式不同；
- 删除前已把每张表的 `sha256`、字节数与行数写入 `derived/manifest.json` 的 `retired_tables`。

| 表 | 行数 | 现行等价物（v2 accepted） |
| --- | --- | --- |
| `cep_boundary_derived_reference_v1.csv` | 161 | `cep_boundary_accepted_phase_map_v1.csv` |
| `crossover_surface_derived_reference_v1.csv` | 3 135 | `crossover_surface_accepted_phase_map_v1.csv` |
| `maxwell_surface_derived_reference_v1.csv` | 12 537 | `maxwell_surface_accepted_phase_map_v1.csv` |
| `spinodal_surface_derived_reference_v1.csv` | 11 989 | `spinodal_surface_accepted_phase_map_v1.csv` |
| `surface_coverage_mask.csv` | 483 | `surface_coverage_mask.csv` |

删除后本层只剩 **30.2 KB 簿记**（原始 33.8 MB）。复原某一版表格的方式：从 v2 accepted 取同名面，
按 `_derived_reference_v1` 的列名重命名；数值不再逐字节可复原，这是本次去重的代价，
已由 `retired_tables.sha256` 记录身份。

## 证据树保持原样（判断说明）

`docs/analysis/pnjl/**` 下的审计树仍引用 v1 的原始路径。这些是**历史审计记录**——记录了当时
审查过什么、哈希是多少。把它们改写成墓碑等于修改审计结论，比保留悬空路径更糟。因此：

- 审计证据树不动；
- 台账 `file:` 证据里指向已删文件的条目改为 `note:` 墓碑（台账会校验文件存在性）；
- 本文件是这批退役的唯一权威记录。
