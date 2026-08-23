# PNJL Phase Surface Series

Issue #130 的 phase-surface 版本系列索引。此目录是 docs-only 的统一 namespace；它不替代原始证据包、phase-reference package 或正式 runtime reference。

## 目录语义

```text
phase_surface_series/
├── SERIES.md
├── series_manifest.json
├── analysis/
│   ├── v1/ ... v7/
│   └── v4_display16/
└── figure_layer/
    ├── phase_surface_render_v1/
    ├── phase_surface_render_v2/
    └── phase_surface_render_v3/
```

`analysis/v1` 至 `analysis/v7` 是 C2 surface 诊断线的连续版本。`analysis/v4_display16` 是 v4 的独立显示变体，不能与 `analysis/v4` 合并或互相覆盖。`figure_layer` 是显示层版本，不是新的数值层：

| 新 namespace | 历史语义 | 内容边界 |
| --- | --- | --- |
| `figure_layer/phase_surface_render_v1/` | render v1 | 从 derived layer 生成的原始 figure-layer 包 |
| `figure_layer/phase_surface_render_v2/` | render v8 | 精修前公开图语义 v8 的 byte-preserving snapshot |
| `figure_layer/phase_surface_render_v3/` | render v9 | 低温 display continuation 和 CEP closure 的显示层候选 |

## Source package contract

为保持现有稳定入口，以下 source package 原位保留，不在本次 namespace PR 中移动：

`docs/analysis/pnjl/c2_surface_views/`

`docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/`

其中 phase-reference package 仍由 `import_issue130_phase_reference.py`、promotion gate 和对应测试直接消费。新目录中的内容是逐文件、逐字节相等的 snapshot；包内 `manifest.json`、checksum、execution log 和生成时 provenance 不重写。源/目标文件数、字节数和 inventory SHA-256 对照见 `series_manifest.json`。

## Reading and promotion boundaries

1. 先读取具体版本目录中的 `README.md`、`manifest.json`、`decision.json` 和 claim ledger。
2. `analysis` 层保留 native unresolved、interpolated/non-certified 等诊断语义；图面闭合不等于数值证书闭合。
3. `figure_layer` 只描述可复现的显示输出；它不修改 strict/derived CSV，也不触发 solver、Maxwell、reference promotion 或 RS transport。
4. 新 namespace 的存在不改变历史路径的 provenance；若未来切换 canonical source，必须另开 PR，同时更新 importer、promotion gate、测试和所有消费者。

## Migration record

- migration type: byte-preserving namespace snapshot
- base: `origin/main@3ee9e03c86924137dfb4df922b5f792bd22fe901`
- numerical solver called: `false`
- scientific payload regenerated: `false`
- source packages removed: `false`
- per-package source/destination inventory equality: `true`
- detailed mapping and hashes: `series_manifest.json`
