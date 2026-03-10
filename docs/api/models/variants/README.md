# Models 变体主题

本目录承接 `Models` 公开表面中不适合直接平铺到 `docs/api/models/` 一级目录的“模型家族/模型变体”主题。

当前已规划主题：

- [magnetic/README.md](magnetic/README.md)：外磁场 PNJL / Landau 能级离散化相关能力。

设计原则：

- 变体主题优先从 `Models` 聚合导出视角组织，而不是继续沿用旧 `PNJL.*` 兼容路径做主导航。
- 变体主题仍遵循三层视图：面向用户入口、职责核心、导出 API 全集。
- 旧 `docs/api/pnjl/` 页面在完成吸收后应降级为迁移说明或跳转页。