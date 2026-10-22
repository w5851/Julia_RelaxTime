# 化学冻结线研究链的 PR 证据

本批收纳默认比较、近期文献筛选、多项式族与 Andronic 族的有限参数探索、
standard/fine 核验和 PNG 审阅。研究结论见
[主报告](../../../../../../docs/analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/README.md)。
证据仍为 `diagnostic_only`，图像为 `png_review`、作者审阅 `pending`。

## 保存与恢复

[archive_manifest.json](archive_manifest.json) 记录原工作树、完整基线提交、
526 个直接保留文件及五个 ZIP 的 SHA-256、字节数和成员数。
四个数值 ZIP 保留 5806 个原始逐点记录、逐点 checksum 和源码快照；第五个 ZIP 保存
三份完整 INSPIRE 检索响应。成员路径相对仓库根目录。
所有成员与原工作树逐字节比较后才打包；CSV、参数表、失败日志、summary、核验记录、
图注、PNG 和图像 manifest 直接入库。原始记录中的绝对路径是历史运行身份，保持原字节。

两组早期 standard 图件的评分、残差和绘图源码已移到对应 `results/` 目录，
`figures/` 仅保留 PNG 与 `plot_manifest.json`。移动的六个文件保持原 hash；
两个当前图像 manifest 仅更正三处路径字符串，原始 manifest 与原归档清单另存原字节。
[布局迁移记录](layout_relocation_20261022.json) 逐项关联旧路径、现路径和 hash，
[原归档清单](archive_manifest.original.json) 保留迁移前的 522 个文件记录。
已有验证记录描述迁移前布局，其数值与图像检查结论仍适用。

在仓库根目录执行：

```powershell
python scripts/analysis/relaxtime/restore_freezeout_review_evidence.py verify
python scripts/analysis/relaxtime/restore_freezeout_review_evidence.py restore
```

`verify` 检查直接保留文件与 ZIP 的字节/hash、成员数、CRC 和路径。
`restore` 先检查所有已有文件冲突，再解压缺失成员；重复执行保持已有相同字节。
存在不同字节时明确报错，不覆盖。恢复的 `points/` 与 `source_snapshot/` 已按批次忽略，
其持久副本为上述 ZIP。复算时用当前仓库下对应的同名批次目录替换历史绝对前缀。

多项式文献 PDF 和检索记录位于
[仓库内来源副本](../../../../../inputs/meson_freezeout_polynomial_literature_20261017/README.md)；
Andronic 来源位于
[原有来源清单](../../../../../inputs/meson_freezeout_andronic_literature_20261018/recent_sources_manifest.json)。
历史 design 保留原 catalog hash；当前 catalog 的 source_root 改为仓库相对路径，
参数值、区间、种子和选择规则一致。

## 本次收尾验证

本次定向软件检查、证据恢复、数值复算和 PNG 检查见
[PR 验证记录](../../../../../../docs/analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/pr_validation_20261022.json)。
恢复后可执行本批保留的复核脚本，重新核对运行 manifest、源码快照、密度商、
固定九点评分、整条曲线成员、standard/fine 排序与解析 T–muB 带：

```powershell
python scripts/analysis/relaxtime/review_freezeout_evidence_20261022.py
```

实际结果写入 [evidence_review.json](evidence_review.json)；原计算、summary 和 PNG 保持原字节，
当前图像 manifest 的路径更正与原始记录通过上述迁移清单关联。
本阶段不安排新增扫描；后续问题由作者组会讨论后独立确定。
