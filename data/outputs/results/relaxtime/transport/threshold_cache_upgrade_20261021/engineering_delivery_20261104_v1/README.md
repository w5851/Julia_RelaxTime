# 输运重计算与结果分析：工程交付（2026-11-04）

本目录整理已经验收的计算和分析，不生成新物理结果。当前数值、图层和研究未决问题统一由
[生产计划](../../../../../../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)维护。
计算源码固定为 `2e5c411c716d5cbc524c3fc6cfb335fe05b39444`；工程收尾前分支 HEAD 为
`3a07a0e2a6082011bc72ea3f9e02cdf55e1bb4b5`。

## 交付范围

- 源码与测试：整批队列、来源复核、η 加密组装、Mott 窗口和 σ 渐近分析、原始图绘制及历史图包兼容。
  `git_paths.txt`列出本次提交的新增/修改/删除路径，`file_manifest.json`记录其中现存文件的Git内容哈希和SHA-256；其余工作树改动保留在原位。
- 数值与图像：原1819 点和η加密1867 点、正式晋升/回放依据、Mott诊断、v14/v15及相关PNG审阅包。
  v14是当前正式图层；v15 PDF暂缓；诊断图不因此获得论文资格。
- 恢复归档：当前输运证据及历史试算、重复源码快照、失败记录均按原路径和字节保存。
  archive manifest逐文件记录SHA-256，恢复命令拒绝覆盖已有目录。
- 不在本次提交内：介子冻结线说明、其他活动任务状态、历史图资产退役改动及benchmark环境的本地Manifest。
  旧资产的退役和清理安排在本阶段收尾之后，现有文件不因归档而删除。

## 验证与恢复

归档位置：`D:/DevData/Julia_RelaxTime/artifacts/transport_recalculation_20261104_v1/transport_evidence.zip`。
这是worktree外的本机持久副本；不声称异机备份。归档摘要见 `archive_receipt.json`：43,645个文件、原始2,716,313,796字节，ZIP为1,278,373,244字节，
SHA-256为`0cffb26b3fc06b724da31be80ddde636a1f9719a62f56bc2aa273b07624aabc8`。
`restore_verification.json`记录43,645个文件已经实际解压并逐文件验证；恢复目录在同一归档目录的`restored_verification/`。

在仓库根目录验证或恢复到一个不存在的新目录：

```powershell
python scripts/analysis/relaxtime/transport_delivery_archive.py verify `
  --archive D:/DevData/Julia_RelaxTime/artifacts/transport_recalculation_20261104_v1/transport_evidence.zip `
  --receipt data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/engineering_delivery_20261104_v1/archive_receipt.json

python scripts/analysis/relaxtime/transport_delivery_archive.py restore `
  --archive D:/DevData/Julia_RelaxTime/artifacts/transport_recalculation_20261104_v1/transport_evidence.zip `
  --receipt data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/engineering_delivery_20261104_v1/archive_receipt.json `
  --destination D:/DevData/Julia_RelaxTime/artifacts/transport_recalculation_restore_new
```

恢复目录以仓库相对路径保存文件；用于完整重放时，与计算提交的独立检出及工程源码包配合。
`engineering_source_receipt.json`记录同目录`engineering_source_v2.zip`的哈希、文件数及适用基底；
它保存本次最终Git内容的源码、测试、配置、说明及恢复工具，不包含其他任务未提交改动。
包内`__transport_delivery__/manifest.json`记录逐文件身份；使用上述工具对源码包执行同样的`verify/restore`命令即可。
`engineering_source_restore.json`记录源码包的实际恢复检查。
`git_paths.txt`和`file_manifest.json`不互相自校验；清单及回执由本次Git提交绑定，回执不会反向写入它所描述的ZIP。
现有图包的代码/合同按冻结快照验证，数值、图像及manifest仍须满足原始哈希。
不以当前可变README或SOP字节替换历史快照，不以重新标记修复历史证据。

`validation/`保存本次治理、文件身份、Python绘图及quote-remote定向分析测试记录。

- quote-remote：59项合成差分/分解及一次性I/O夹具检查通过，168份源码运行前后及交付时一致。
- 本机Python绘图/I/O：首次92项中91项通过，旧raw图包因系统字体更新失败；
  取回同字节历史字体后发现冻结测试源码未被识别，修复后19项定向复测全部通过。原失败日志按原字节保留在`validation/initial_failure_logs.zip`。
- 兼容改动：验证器允许显式`--font-sources`映射到相同字节字体，继续校验原哈希、大小与字体族；
  测试脚本可使用精确冻结源码快照，基线CSV和图像仍按原数据身份校验。
- 原始3,490份Git产物与外部归档逐字节一致；任务台账、活动文档、分析脚本治理和文档链接完成复核。

旧字体的本机保留路径见`validation/historical_font_receipt.json`，只作为本机恢复依赖，不提交或分发字体。
旧raw图包的检查使用`PLOT_FONT_SOURCES`指向`validation/historical_font_sources.json`；
重定位仅适用于期望该字体哈希的历史包，不能用旧字体覆盖采用新字体的图包。
生产核和基线未改变，沿用已通过的分层、代表点、314项回归和整批数值验收；本次没有重跑物理扫描。

## 后续顺序

当前状态统一、恢复核验、源码/测试与证据的选择性提交构成本次本地工程交付；
本记录不表示已推送、建立PR或合并。2026-11-04的SSH fetch因publickey权限失败；
只读GitHub API确认远端同名分支仍为`16885cfeb3e01351064e5eb2f2e06acd57aa83e5`，没有对应PR。
本地提交和来源恢复包就绪后，再处理旧资产。
旧资产已有盘点仅是候选范围；执行时仍须复核引用、实际路径和恢复链，不能借本交付清单直接删除。
完整输运量的严格光滑性、π最窄尖峰、混合谱支排除、A解析降维与v15 PDF保留各自范围。
