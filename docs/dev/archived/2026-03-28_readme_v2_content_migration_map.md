---
title: README v2 内容迁移映射（实施工件）
archived: true
original: docs/dev/active/2026-03-28_readme_v2_content_migration_map.md
archived_date: 2026-03-30
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# README v2 内容迁移映射（实施工件）

> 用途：在 README 重构期间记录“保留/迁移/下沉”决策，避免把历史长内容重新堆回首页。

## 1) 旧内容 -> 新去向映射

| 旧 README 内容块 | 新处理方式 | 目标位置 |
| --- | --- | --- |
| 大段单位约定与物理量说明 | 下沉（保留入口） | `docs/reference/`（README 仅链接） |
| 长篇“当前功能概览”逐模块细节 | 收敛为能力边界 + 稳定入口矩阵 | `README.md`（精简版）+ `docs/guides/scripts/README.md` |
| 计算链路长解释 | 下沉 | `docs/reference/` / `docs/api/` |
| 历史更新日志（2025/2026 条目） | 下沉 | Release 页面 + `docs/guides/STATUS.md` |
| 架构迁移长说明 | 下沉（README 仅导航） | `docs/architecture/` |
| 目录迁移对照大表 | 下沉 | `docs/dev/` / `docs/guides/` |

## 2) README v2 八段预算（<=250 行）

| 段落 | 预算上限 |
| --- | --- |
| 1. 项目定位 + Release | 15 行 |
| 2. 三步 Quickstart（脚本优先） | 55 行 |
| 3. 稳定脚本入口矩阵 | 35 行 |
| 4. 能力边界与非目标 | 30 行 |
| 5. 推荐验证命令 | 30 行 |
| 6. 文档导航 | 25 行 |
| 7. 仓库结构极简图 | 35 行 |
| 8. 贡献与治理入口 | 20 行 |
| **总计** | **245 行** |

## 3) AC 覆盖核对矩阵

| AC | 核对项 | 证据 |
| --- | --- | --- |
| AC-1 | 八段结构 + <=250 行 | README 行数统计 + 章节标题 |
| AC-2 | 1 主命令块 + 1 可选 server 命令 | Quickstart 小节 |
| AC-3 | 稳定入口白名单 + 非用户目录声明 | Section 3/4 |
| AC-4 | `releases/latest` 在前 20 行 | README 顶部链接 |
| AC-5 | 验证区包含 artifact + maintainer checks | Section 5 |
| AC-6 | 架构仅导航（短块）并外链深文档 | Section 7 + docs 链接 |
| AC-7 | quickstart 的 git-hygiene（清理或忽略） | Quickstart/验证区说明 |
| AC-8 | 主命令写明 5-10 分钟预期 | Quickstart 文案 |
