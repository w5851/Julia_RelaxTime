# 输运重计算 PR 交付整理（2026-11-04）

本页记录原工程交付之后的分支隔离和存储修正。科学结论与后续范围仍由
[生产计划](../../../../../../../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)维护；
[原工程交付](../README.md)及其清单描述当时的原始提交，不回写其来源身份。

## 分支和工作范围

PR 使用 `codex/transport-recalculation-wrapup`，对应工作树为
`D:/w/transport-recalculation-wrapup`。它从已发布的 `16885cfeb3e01351064e5eb2f2e06acd57aa83e5`
继续，保留计算提交 `2e5c411c716d5cbc524c3fc6cfb335fe05b39444`，并转入两份输运交付提交：

| 原提交 | PR 分支提交 | 内容 |
| --- | --- | --- |
| `4b7c285226d6d0546fb3f03007856c1976074c4f` | `8c7dfded` | 工程交付、生产和诊断证据 |
| `08db84bbe26256e5c93760dcb394bb4569b1ddb6` | `1ea2778e9ed9474592b8a3bfae5a17c63321b905` | v1–v4 公共重复 PNG 退役 |

独立的 Issue130 工作流退役提交 `3a07a0e2a6082011bc72ea3f9e02cdf55e1bb4b5` 未带入本 PR。
`docs/guides/scripts/README.md` 的合并保留对应分支内容；原交付的 `file_manifest.json` 继续用于验证原工程源码包，
不是新分支全部文件的清单。新分支的生产 `src` 树与计算提交相同。

原工作树 `D:/w/meson-freezeout-defaults` 保留在 `08db84bb`，含 8 个未提交 tracked 路径和
39,787 个 untracked 文件：39,767 个输运证据、18 个介子数密度文件、2 个其他文件。
其中 39,765 个输运文件已在外部证据包中验证，另外两份失败日志已在
`../validation/initial_failure_logs.zip` 逐字节核对。归档不改变它们在原树中的未跟踪状态。
介子冻结线／两族参数包络的主要实现已由 [PR #327](https://github.com/w5851/Julia_RelaxTime/pull/327)
合并；其后续本地说明、其他任务状态与该独立提交均不借本 PR 提交或清除。

2026-11-04 再检查时，GitHub SSH 的 `ls-remote`、`fetch` 和目标分支的 `push --dry-run` 均成功。
这更新了原交付记录中当时的 `publickey` 失败状态；实际写入及最终合并状态以对应 PR 和远端 ref 为准。

## 历史图片定位与恢复

PR 范围的路径门禁发现 70 张历史审阅／诊断 PNG 位于 results，另有 4 份退役 README 位于 figures。
本次保持门禁规则，将 PNG 按原路径保存至 [review_path_migration.zip](review_path_migration.zip)，
并将退役说明迁到 `data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v1` 至 `publication_clean_v4`。
这只是证据存储整理，没有重绘、数值删改、新增视觉接受或图层晋升。

ZIP 共 76 个成员：70 张 PNG、4 份原退役说明及 2 份需要修复图片链接的原分析文档。
每个成员的原始路径、大小和 SHA-256 均保存在 ZIP 内的 `__transport_delivery__/manifest.json`；
[archive_receipt.json](archive_receipt.json)绑定整个 ZIP，[restore_verification.json](restore_verification.json)
记录实际解压并逐成员比对成功。旧审阅 JSON、DS 返回记录和工程清单保留原字节。
[migration.json](migration.json)列明图片成员、所查直接引用及 35 份保持原字节的公共 plot manifest。
相对路径或绝对本机路径出现在旧审阅记录中时，含义仍是生成当时的位置；恢复包使用仓库相对路径。

两份分析文档只将已封存的内嵌图片链接改为本页恢复入口；数值、解释及原文件字节分别保持现行正文和包内快照。
常用的三张历史诊断图成员是：

- `data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/analysis_v2/ablation_tau.png`
- `data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/analysis_v2/weighted_threshold_shape.png`
- `data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/followup_analysis_v1/finite_x_resonance.png`

在仓库根目录验证或恢复到尚不存在的目录：

```powershell
$package = 'data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/engineering_delivery_20261104_v1/pr_wrapup'
python scripts/analysis/relaxtime/transport_delivery_archive.py verify --archive "$package/review_path_migration.zip" --receipt "$package/archive_receipt.json"
python scripts/analysis/relaxtime/transport_delivery_archive.py restore --archive "$package/review_path_migration.zip" --receipt "$package/archive_receipt.json" --destination 'tmp/transport_review_restore_new'
```

恢复工具拒绝覆盖已有目录。需要重放旧审阅时，在独立恢复目录中组合原记录和这些原路径成员；
不要把恢复出的旧 README 覆盖当前说明。正式数值 CSV、v14/v15 正式或待审图包、current 指针及依赖它们的来源记录未迁出。

## 验证与下一阶段

原始数值资格和工程测试见 [validation_summary.json](../validation_summary.json)。PR 包含生产源码相对 main 的变化，
因此仍需 GitHub 的 unit、integration、regression 和 validation 检查；本次存储整理没有重跑物理数值。
路径整理的本地验证记录见 [validation.json](validation.json)；最终 CI 与合并状态以 GitHub PR 为准。

后续先完成 v15 的 PDF 交付检查与图层晋升，再更新 current、论文插图及相应说明并编译核对，
最后复核旧资产的消费者和恢复链。旧数值、收敛证据或代码依赖仍被 v15 引用时继续保留，不能按版本号整体删除。
