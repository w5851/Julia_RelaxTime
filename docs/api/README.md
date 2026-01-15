# docs/api

这里是项目的 API 文档入口。

## 组织方式

- `docs/api/*.md`：跨模块或 relaxtime 主链路相关的 API 文档（保持历史链接兼容）。
- `docs/api/pnjl/*.md`：PNJL 求解/扫描/工作流相关 API 文档。
- `docs/api/relaxtime/`：预留子目录（见该目录下 README）。

## 命名约定

- 文档文件名尽量与 `src/` 中的 Julia 源文件同名（如 `TotalCrossSection.jl` → `TotalCrossSection.md`），便于从名字直接定位实现。
- 如果某 API 属于 PNJL 子模块，优先放入 `docs/api/pnjl/`，避免根目录膨胀。

如果后续需要做目录重构，建议采用“新旧路径并存 + 旧路径保留跳转说明”的方式，避免破坏外部引用。
