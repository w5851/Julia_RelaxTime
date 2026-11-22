# 架构索引

本域只回答「代码怎么组织、依赖朝哪走」。接口契约去 [docs/api](../api/README.md)，选型理由去 [docs/decisions](../decisions/README.md)。

## 首读

- [agent_overview.md](agent_overview.md)：include-driven 装配方式、模型类型层级、求解器、相图流水线、RelaxTime 链路、工作流引擎与配置系统。

## 模块契约

- [models_solver_contract.md](models_solver_contract.md)：求解器与约束求解的模块契约。
- [models_param_container_contract.md](models_param_container_contract.md)：参数容器的构造与校验约定。
- [models_ad_implicit_contract.md](models_ad_implicit_contract.md)：隐式微分与 AD 后端约定。

- [Julia模块配置文件设计指南.md](Julia模块配置文件设计指南.md)：2026-01 的设计参考；现行加载与配置契约以 [agent_overview.md](agent_overview.md) 为准。

## 依赖

- [dependency_rules.md](dependency_rules.md)：依赖方向的硬规则。
- [dependency_review.md](dependency_review.md)：依赖评审记录。
- [dependencies.md](dependencies.md)、[dependencies.manual.md](dependencies.manual.md)：依赖清单（生成与手工补充）。
- `dependencies.mmd`、`dependencies.svg`：依赖图的源与渲染产物。
- [qcd_refactoring_dependencies.md](qcd_refactoring_dependencies.md)：重构期的依赖迁移记录。

## 边界

- 需要「为什么这样选」→ `docs/decisions/`（ADR）。
- 需要「接口怎么调用、字段是什么」→ `docs/api/`。
- 对外稳定契约进 `docs/api/`；跨子系统内部契约留在本域。
- 一次性排障过程不要写在这里，留在对话或 PR（`docs/dev/archived/` 是只读历史区，不接新内容）。
