# Rotation 变体主题 API

本目录是 `Models` 变体层中 rotation 主题的主入口，面向旋转背景下的 `RotationModel` 与单点 workflow 公开能力。

推荐阅读顺序：

1. [Overview.md](Overview.md)：先看 `RotationModel`、`solve_rotation_point` 与最短调用路径
2. [CoreConcepts.md](CoreConcepts.md)：理解模型适配层与旋转热力学核心层职责边界
3. [generated/Exports.md](generated/Exports.md)：自动生成的公开导出全集

本主题覆盖的 `Models` 公开表面包括：

- `RotationModel`
- `solve_rotation_point`

说明：

- `RotationThermo` 核心实现位于 `src/models/variants/rotation/core/RotationThermo.jl`，但该页不作为用户首入口。
- 用户入口优先通过 `Models` 聚合导出访问，而不是直接 include core 文件。
