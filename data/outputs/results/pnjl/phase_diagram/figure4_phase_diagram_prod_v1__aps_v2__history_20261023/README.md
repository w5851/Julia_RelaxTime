# PNJL 相图 APS 前期图包归档

2026-10-23 留存三个 2026-10-08 APS 阶段：首轮 PNG、r2 PNG、r2 PDF，共 48 文件（含 PDF）。
它们是 PR #326 已采用的 PRD v3 图包之前的历史阶段；原 PNG/PDF、manifest、接受说明和实际源码快照全部保留。
三个包均为 audit、manuscript_eligible=false、current_publication_layer=false；历史作者接受记录不代表当前投稿资格。

`original_aps_packages.zip` 保存所有原始文件，以仓库相对路径为成员名。
`archive_manifest.json` 列出 ZIP 和每个成员的 SHA-256、字节数、原目录与恢复方式；完整 48 文件已在隔离目录恢复并逐字节核验。
已入库 PRD v3 的两个 ZIP 仅覆盖这些包中的 9 个文件字节；另 39 个文件不能作为重复物删除，本归档完整保留。
核验通过后只移除三个原始 APS 目录的明确文件清单及空目录；未重绘，未修改已采用 PRD v3 图包。

恢复时先核对 archive_manifest.json 中的 ZIP hash，再解压到新的隔离目录并逐项核对成员；
不要覆盖当前 figures case。原 manifest 的路径与源码 hash 保持原字节，解释以当时阶段为准。
