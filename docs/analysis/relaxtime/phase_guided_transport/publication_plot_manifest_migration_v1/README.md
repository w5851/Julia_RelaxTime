# v10/v11 Manifest 存储迁移

本包记录作者授权的 metadata-only 迁移。三个图层目录现在各只保留一份
`plot_manifest.json`：v10 PNG、v11 PNG、v11 PDF。共同 provenance 存一次，
222 条逐图记录完整保留，未重绘图像或调用求解器。

`original_manifest_graph.zip` 保留迁移前 225 份 manifest/index 的原始字节，
由 `config/plotting/historical_snapshots.toml` 登记 SHA-256。旧验收与 analysis
package 的 manifest 引用继续核对这些冻结字节；PNG、PDF、CSV 和数值 provenance
仍核对现存文件，不能回退到历史数据。

`manifest.json` 记录三个新总 manifest 的 hash、旧 index 引用和未改变产物清单。
检查同时证明全部逐图字段与旧 manifest 等价，包括尚未解决的 PDF 字形尺寸门禁。
原 v11 阶段性接受记录、`manuscript_eligible=false` 和 current=v5 均保持不变。

```powershell
python scripts/analysis/relaxtime/migrate_phase_guided_plot_manifest_bundles.py --check
python scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py --check
```

恢复历史证据时读取 ZIP 内原仓库相对路径；不要覆盖现行图包。恢复旧 renderer 的
合同另见代码 snapshot registry。存储迁移不是新的视觉、数值或论文资格接受。
