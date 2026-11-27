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
| `derived/tables/`（5 文件，6.09 MB） | v2 的 `internal_build_inputs` 指向本目录 `derived`，且这是该封装的**唯一副本**（v2 的 `accepted/` 是另一种列名与格式，数值相同但非字节副本） |
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
