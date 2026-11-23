# 指南总入口

面向前两类读者的操作手册：**人类使用者**与**运行者/维护者**。这里只回答「怎么跑通」，不写科研诊断过程（[docs/analysis](../analysis/README.md)）与接口契约（[docs/api](../api/README.md)）。

## 首读

- [QUICKSTART.md](QUICKSTART.md)：最小跑通路径。
- [USER_GUIDE.md](USER_GUIDE.md)：功能与参数说明。
- [STATUS.md](STATUS.md)：当前链路状态与已知限制。

## 稳定入口与 SOP

- [scripts/README.md](scripts/README.md)：稳定脚本白名单（唯一权威）；全量脚本目录见 [scripts/run_script_catalog.md](scripts/run_script_catalog.md)。
- [sop/README.md](sop/README.md)：正式生产与绘图 SOP 总入口；具体流程见 `sop/workflows/`（相图、输运、介子密度、介子热力学、绘图、GasLiquid/RMF）。
- [examples/](examples/)：可直接运行或改参数的最小示例（`phase_diagram.md`、`web_demo.md` 等）。

## 主题指南

- 配置与环境：[configuration.md](configuration.md)
- 错误处理约定：[error_handling.md](error_handling.md)
- 前端调试：[FRONTEND_DEBUG.md](FRONTEND_DEBUG.md)
- 基准与趋势：[benchmark_trends.md](benchmark_trends.md)
- 基线版本管理：[BASELINE_VERSION_MANAGEMENT.md](BASELINE_VERSION_MANAGEMENT.md)
- 参数结构迁移：[PARAMETER_STRUCT_MIGRATION.md](PARAMETER_STRUCT_MIGRATION.md)
- 缺陷与修复记录：[FIXES.md](FIXES.md)
- 发布治理（带版本号，注意时效）：[RELEASE_GOVERNANCE_v0.1.1.md](RELEASE_GOVERNANCE_v0.1.1.md)

## 边界

- 某次跑批的状态快照与诊断细节属于 `docs/analysis/`；本域只在 SOP 中引用结论，不重复维护状态。
- 命令的精确参数以 [命令目录](../dev/agent_command_reference.md) 为准。
