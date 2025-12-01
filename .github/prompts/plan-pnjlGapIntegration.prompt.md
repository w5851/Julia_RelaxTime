## 计划：集成 PNJL 能隙方程 UI

将外部 PNJL 能隙方程后端引入 `src/`，与现有的松弛时间模块并置，通过当前的 HTTP 服务器进行暴露，并扩展 `web/` 前端，使用户可以在传统界面和新的 PNJL 视图之间切换，同时共享底层 API。

### 步骤
1. 审计 `src/relaxtime/*.jl`、`src/simulation/HTTPServer.jl` 和 `server/server.jl` 中现有的松弛时间流程，绘制请求处理流程并定位用于新增 PNJL 模块的挂钩点。
2. 定义外部 PNJL 后端（输入/输出、依赖项、数据布局）在 `src/` 内的落地方式（例如 `src/pnjl/`），并确保复用共享的常量与工具函数。
3. 扩展服务器路由（`server/server_full.jl`、`server/test_server.jl`），以公开新的 PNJL 接口，将响应连接到能隙求解器，并在 `docs/reference/domain-knowledge/` 中记录接口契约。
4. 更新前端资源（`web/index.html`、`web/js/*`、`web/css/style.css`），添加 UI 模式选择器和专用的 PNJL 面板，这些面板调用新接口且不破坏现有工作流。
5. 添加文档与测试（`docs/guides/`、`tests/unit/`），涵盖 PNJL 求解器的使用、API 载荷，以及针对松弛时间计算的回归检查。

### 需要进一步考虑的问题
1. 是否需要外部 PNJL 后端的仓库路径以及任何非标准依赖？选项 A：作为子模块暴露 / 选项 B：内联移植 / 选项 C：以包（package）形式消费。
2. 新的 UI 应该沿用现有样式，还是采用独立的布局？
