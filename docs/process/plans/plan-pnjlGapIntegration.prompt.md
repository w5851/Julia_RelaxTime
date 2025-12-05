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

### 前端需求描述（Excel 风格表格与模块选择器）

需求概述：在前端实现一个 Excel-like 的多表格视图：页面底部显示标签页（sheet tabs），每个标签对应一组结果/表格；用户可通过底部标签切换表格、增加/删除/重命名表格。与此同时，提供一个模块选择器（Module Selector），用户可以按需勾选或排列计算模块（例如：PNJL gap、松弛时间、差分截面等），系统根据所选模块触发相应计算并将结果保存到可切换的表格中。

关键交互点：
- 底部标签（Sheets）：显示表名、可新建/重命名/删除、拖拽排序，显示是否有未保存或计算中状态。
- 模块选择器：以侧边栏或顶部下拉形式出现，可展开模块的参数面板，用户能设置参数并触发“运行”，运行结果会新增或更新当前 Sheet。
- 表格功能：支持虚拟化以处理大表、列排序、过滤、单元格复制/粘贴、导出 CSV/Excel、可选只读或可编辑模式。
- 状态与持久化：当前工作空间（Sheets + 选中模块 + 参数）保存在浏览器 `localStorage` 或后端用户会话，以便刷新后恢复。

API 约定（示例）：
- 列出可用模块：GET `/api/modules` -> `[{id, name, description, params_schema}]`
- 触发计算：POST `/api/modules/{id}/run`
	- 请求体示例：
		```json
		{
			"params": { "T": 0.15, "mu": 0.1 },
			"target_sheet": "PNJL_result_1",
			"save_as_new_sheet": true
		}
		```
	- 返回示例：`{ "status":"queued", "job_id":"abc123", "sheet_id":"PNJL_result_1" }`
- 查询结果：GET `/api/results/{sheet_id}` -> `{ "columns": [...], "rows": [...], "meta": {...} }`
	- 可选：使用 WebSocket `/ws/results` 推送进度与流式数据。

前端实现建议（技术/库）：
- 表格：考虑 `ag-Grid`, `Handsontable`, `Tabulator` 等支持虚拟化与导出的库。
- 底部 tabs：一组可编辑标签实现新建/重命名/删除（双击重命名）。
- 参数表单：可使用 JSON schema 表单库以自动生成参数界面。
- 持久化：`localStorage` 用于快速原型；如需多人或跨设备，后端保存配置。

验收条件（Acceptance Criteria）：
- 用户能在页面底部切换至少 3 个 Sheet 并能新建/删除/重命名。
- 从模块选择器触发计算后，能看到“计算中”状态，并在完成后把结果写入选定 Sheet 或创建新 Sheet。
- 表格支持导出 CSV/Excel，并能在刷新后恢复上次会话状态。

简短需求语句（用于 issue）：
"实现底部标签式 Sheet 切换（Excel 风格）与模块选择器：用户可新建/重命名/删除表格，选择并配置计算模块，运行后将结果保存到选定表格。支持导出、状态持久化与计算进度显示（轮询或 websocket）。"

### 外部仓库接入建议（问题一）

结论建议：为了让 VS Code agent 能直接读取并分析外部项目源码，推荐在开发与集成阶段将外部 PNJL 仓库以 `git submodule` 或本地 path 依赖（multi-root workspace / `dev` 模式）纳入工作区；其次是将代码内联到 `src/pnjl/`。仅通过 registry/远程包（不把源码放到工作区）不利于 agent 静态分析与快速调试。

原因简述：VS Code agent 能访问的是工作区的文件系统。子模块或 path 依赖会把源码放在仓库或本地路径下，agent 可以打开、检索符号、运行本地测试或做重构；而仅通过远程包管理器安装（若未把源码检出到工作区）时，agent 无法直接审阅实现细节。

优缺点对比：
- 子模块（推荐在外部仓库有独立生命周期时）：
	- 优点：保留独立提交历史、便于同步上游、源码可被 agent 直接读取与 diff、便于单独更新。
	- 缺点：需要维护 submodule 更新流程（`git submodule update --init --recursive`）。
	- 示例命令（PowerShell）：
		```powershell
		git submodule add <git-url> src/pnjl
		git submodule update --init --recursive
		```
- 本地 path 依赖 / dev 模式（推荐并行开发）：
	- 优点：外部项目作为工作区的一部分，agent 可跨项目引用并进行分析；在 Julia 中可使用 `pkg> dev /path/to/pnjl`。
	- 缺点：需要管理路径与 workspace 配置。
	- Julia 示例：
		```julia
		pkg> dev /absolute/path/to/pnjl
		```
	- 或在 `Project.toml` 中使用 `path` 依赖：
		```toml
		[deps]
		PNJL = { path = "../pnjl" }
		```
- 内联移植（把外部代码拷贝进 `src/pnjl/`）：
	- 优点：集成最直接、agent 访问无障碍。
	- 缺点：造成代码重复，需要自行维护与 upstream 的同步策略。
- 仅通过包管理器安装（Option C）：
	- 优点：依赖管理清晰、易于版本化。
	- 缺点：若源码不在工作区，agent 无法直接读取实现细节（除非将包以 `dev` 或源代码形式放入 `.julia/dev` 等可见路径）。

私有仓库注意事项：若外部仓库为私有，子模块或 path 方式仍可行，但需提前准备凭证（SSH key 或 Personal Access Token）；计划中应注明仓库可见性以便配置 CI/开发环境。

可直接复制到 issue/计划中的一句话（简短版本）：
"为便于 VS Code agent 阅读与调试，开发阶段请将 PNJL 源码以 `git submodule` 或本地 path 依赖方式纳入工作区；长期可同时发布为独立包以利版本管理。"
