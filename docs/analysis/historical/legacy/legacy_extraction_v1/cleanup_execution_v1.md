# Legacy 树删除执行记录（2026-11-27）

`D:\Desktop\legacy` 下的八个外部验证树已由作者手工删除。本文件记录删除前后的核查，
以及本批**未能完成**的一项验证。

## 删前核查

按 `tables/legacy_extraction_index_v1.csv` 逐目录比对实际文件：

| 目录 | 实际文件 | 记录文件 | 字节相符 |
| --- | --- | --- | --- |
| `fortran-archive-validation` | 488 | 488 | ✅ |
| `fortran-gas-liquid-validation` | 488 | 488 | ✅ |
| `fortran-julia-transport-validation` | 26 | 26 | ✅ |
| `fortran-magnetic-field-validation` | 416 | 416 | ✅ |
| `fortran-meson-validation` | 21 | 21 | ✅ |
| `fortran-mott-validation` | 17 | 17 | ✅ |
| `fortran-relaxtime-validation` | 1663 | 1663 | ✅ |
| `fortran-rmft-validation` | 65 | 65 | ✅ |

八项文件数与字节数全部相等，与 `manifest.json` 记录的
`inventory = {directories: 8, files: 3184, bytes: 717610994}` 一致。

## 未能完成：树哈希核对

`manifest.json` 给出的规则是「按相对路径排序的 UTF-8 行
`relative_path|file_sha256|file_size`，以 LF 连接，无尾随 LF」。按该规则以及 32 种
常见构造变体（相对根/相对子目录、正反斜杠、带目录前缀、带尾换行、是否二次排序）
计算，均无法复现 `tables/legacy_extraction_index_v1.csv` 中任何一条 `tree_sha256`。

因此本批的核查**只到文件数与字节数**这一层，不声称做过哈希级核对。若今后再有同类
提取包，需要先修好这条规则的可复现性，否则 cleanup gate 的「verify the recorded
inventory」无法真正满足。

## 删除后

- 八个子目录与其父目录 `D:\Desktop\legacy` 均已不存在；
- 本提取包保留为唯一溯源记录，`admission.cleanup_status` 已更新为 `deleted_2026-11-27`；
- 残余风险：旧二进制与未入册扫描不再可复查。`tables/` 中已准入的结论不再依赖这些树，
  现行生产输入也没有引用它们。
