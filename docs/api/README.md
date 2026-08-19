# docs/api

这里是 API 文档总入口。对外导航只保留“从稳定入口到主题文档”的最短路径；脚本白名单与当前状态分别以 `README.md` 和 `docs/guides/` 为准。

需要按完整计算流程查看“已实现能力、公式信源、代码入口、测试层级与验证边界”时，先读 [已实现计算能力与方法追踪清单](../reference/implemented_capabilities.md)。

## 1. 先按你要做的事找入口

### 跑相图 / 看相结构工作流

- 运行入口：`scripts/pnjl/calculate_phase_structure.jl`
- 主题入口：`docs/api/models/phase/README.md`
- 首读页面：`docs/api/models/phase/Overview.md`

### 跑 T-μ / T-ρ 扫描

- 运行入口：`scripts/models/run_unified_scan.jl`
- 主题入口：`docs/api/models/scans/README.md`
- 首读页面：`docs/api/models/scans/Overview.md`

### 看 `Models` 统一编排与 workflow

- 主题入口：`docs/api/models/workflows/README.md`
- 首读页面：`docs/api/models/workflows/Overview.md`

### 看核心求解器与约束

- 主题入口：`docs/api/models/solver/README.md`
- 首读页面：`docs/api/models/solver/Overview.md`

### 看输运 / 弛豫时间链路

- 主题总入口：`docs/api/relaxtime/README.md`
- transport 主题：`docs/api/relaxtime/transport/README.md`

### 看模型变体或衍生量

- 模型变体总入口：`docs/api/models/variants/README.md`
- 衍生量总入口：`docs/api/models/derived/README.md`

## 2. 当前主导航口径

- `docs/api/models/phase/`、`docs/api/models/scans/`、`docs/api/models/workflows/`、`docs/api/models/solver/` 是 `Models` 主链的核心主题。
- `docs/api/relaxtime/transport/` 是 transport 算法与 provider 的主主题；不再额外平铺为新的 `docs/api/models/*` 一级主题。
- `docs/api/pnjl/PNJL.md` 仅保留兼容层聚合说明；现行求解、扫描、相图入口以 `docs/api/models/*` 为准。

## 3. 每个主题内部怎么读

主题目录默认遵循三层视图：

1. `README.md`：主题入口与阅读顺序
2. `Overview.md`：面向调用方的用户入口
3. `Algorithms.md` 或 `CoreConcepts.md`：职责核心层
4. `generated/Exports.md`：自动生成的导出 API 全集

## 4. 维护约定

- 同一事实尽量只保留一个 authoritative 页面；根 README 负责“跑什么、从哪进”，API 文档负责“接口做什么、怎么组织”。
- 新 `Models` 能力优先吸收到对应主题，不再把旧 `docs/api/pnjl/` 页面当并列主入口。
- 模型变体类与衍生量类主题优先进入分类总层，例如 `models/variants/`、`models/derived/`，避免 `docs/api/models/` 一级目录继续横向膨胀。
