# v10/v11 Manifest 存储迁移

本包记录作者授权的 metadata-only 迁移。三个图层目录现在各只保留一份
`plot_manifest.json`：v10 PNG、v11 PNG、v11 PDF。共同 provenance 存一次，
222 条逐图记录完整保留，未重绘图像或调用求解器。

`original_manifest_graph.zip` 保留迁移前 225 份 manifest/index 的原始字节，
由 `config/plotting/historical_snapshots.toml` 登记 SHA-256。旧验收与 analysis
package 的 manifest 引用继续核对这些冻结字节；PNG、PDF、CSV 和数值 provenance
仍核对现存文件，不能回退到历史数据。

`original_package_metadata.zip` 另保留六份原 package、验收、current 和 placement
元数据的原 Git 字节（约 54 kB），与原代码快照共同闭合源依赖。即使 squash 合并、
删除分支并清理旧 Git 对象，v11 snapshot 检查也不依赖可达的旧 commit。
这不是历史数值数据 fallback；当前数值和图像仍须通过现存字节检查。

`manifest.json` 记录三个新总 manifest 的 hash、旧 index 引用和未改变产物清单。
检查同时证明全部逐图字段与旧 manifest 等价，包括尚未解决的 PDF 字形尺寸门禁。
原 v11 阶段性接受记录、`manuscript_eligible=false` 和 current=v5 均保持不变。

```powershell
python scripts/analysis/relaxtime/migrate_phase_guided_plot_manifest_bundles.py --check
python scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py --check
```

恢复历史证据时读取 ZIP 内原仓库相对路径；不要覆盖现行图包。恢复旧 renderer 的
合同另见代码 snapshot registry。存储迁移不是新的视觉、数值或论文资格接受。
