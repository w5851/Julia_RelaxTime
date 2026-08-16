# Figure Asset Registry v1

这是历史 PNG/PDF/SVG 退役流程的 PR A dry-run 报告。它只盘点 Git 已跟踪的 `data/outputs/figures` 资产；当前工作树中的未跟踪 C1/C2/pilot 文件被排除且未修改。

## 文件

- `asset_registry.json`：逐文件 hash、格式元数据、manifest/provenance 关联、仓库内文本引用和 review group。
- `cleanup_candidates.csv`：按 case/variant 聚合的人工审核项。
- `cleanup_preflight_v1.json`：作者审核后的执行前快照，仍保留 `delete_performed=false`、`move_performed=false` 和 `rename_performed=false`。
- `retirement_execution_v1.json`：PR B 实际 allowlist 执行结果，记录删除前 hash、迁移/改名前后路径与 hash，以及引用更新验证。
- `relocation_execution_v1.json`：PR B 后续 namespace 整理记录，记录 P1 Mott production case 的源/目标路径、manifest 路径更新和字节级 hash 保持。

## Namespace follow-up

P1 Mott/isentropic production case 的正式图像已从
`data/outputs/figures/relaxtime/paper_p1_mott_phase_isentropic/production_20260531/`
迁移到
`data/outputs/figures/relaxtime/mott_phase/paper_p1_mott_phase_isentropic/production_20260531/`。
这是图像 namespace 整理，不是数值迁移：`data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/production_20260531/`
和 `data/reference/pnjl/paper_p1_mott_phase_isentropic_20260531/` 保持原路径和内容不变。

`asset_registry.json`、`cleanup_candidates.csv`、`cleanup_preflight_v1.json`、`author_review_decisions.md`
和 `retirement_execution_v1.json` 是清理流程的历史快照，保留迁移前的源路径；当前 canonical 路径及迁移 hash 以
`relocation_execution_v1.json` 和目标目录下的 `plot_manifest.json` 为准。

## 语义

- `contract_case`：已有 `plot_manifest_v1` 的新合同 case，默认保留。
- `legacy_manifest`：有历史 manifest 但不符合新 schema，默认保留为 `legacy`。
- `legacy_unregistered`：没有可识别的新 manifest，必须先确认论文/报告引用。
- `owner_review_only`：只表示需要作者判断，不表示允许删除或移动。

仓库内引用由脚本扫描 Git tracked 文本文件得到；论文、supplement、报告和幻灯片等仓库外引用仍为未知。registry 明确记录 `delete_performed=false`、`move_performed=false`，本报告不提供删除接口。

PR B 的执行 manifest 由 `scripts/plotting/build_figure_retirement_manifest.py` 依据 preflight 快照验证生成；它不执行文件操作。

## 重建

在仓库根目录运行：

```powershell
python scripts/plotting/inventory_figure_assets.py `
  --root data/outputs/figures `
  --registry docs/analysis/figure_asset_registry_v1/asset_registry.json `
  --candidates docs/analysis/figure_asset_registry_v1/cleanup_candidates.csv `
  --overwrite
```

如需仅在本机检查未跟踪 pilot，显式加入 `--include-untracked`，但不得把该结果直接作为 PR B 的删除清单。
