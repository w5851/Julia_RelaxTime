# 旧公共显示副本退役（2026-11-04）

本批在工程交付提交`4b7c285226d6d0546fb3f03007856c1976074c4f`完成后执行（提交回执见`engineering_commit_receipt.json`，可恢复Git bundle见`engineering_git_bundle_receipt.json`），范围只包含
`D:/w/meson-freezeout-defaults/data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1`至`publication_clean_v4`的288张PNG。
每张图都在同版本分析包`figures/mode_a|mode_b/`有逐字节相同的保留副本；公共目录原manifest保持字节不变。
执行结果见`execution.json`，完整逐文件范围及保留源映射见`preflight.json`。

恢复包为`D:/DevData/Julia_RelaxTime/artifacts/transport_recalculation_20261104_v1/legacy_public_mirrors_v1_v4.zip`，
包含288张PNG和4份原manifest。`archive_receipt.json`绑定ZIP与每个成员，`restore_verification.json`记录实际恢复验证。
如需运行旧公共图层验证或历史formalization，先用已提交的`transport_delivery_archive.py restore`解压到不存在的新目录，
按清单核对后恢复旧公共路径。正常阅读历史图片直接使用分析包中的同字节副本。
这是本机副本和Git历史恢复，不声称异机备份；旧manifest中的原路径描述退役前状态，不是当前完整公共图包。

保留范围：v5公共PNG被v13的输入manifest显式绑定；v10–v13图包、全部分析包、旧raw网格、收敛证据、
数值基线、Mott诊断、桌面另一个checkout及论文附件保持原位。后续若迁出这些有依赖的资产，须先完成消费者切换和恢复验证。
当前数值为1867点正式版本，正式图层仍为v14，v15 PDF暂缓；本批不改变上述资格。

清理后43项历史绘图与产物依赖检查全部通过，台账路径检查通过；current、v14正式化与v14/v15图包哈希不变，见`validation.json`、`governance.json`和`integrity.json`。
