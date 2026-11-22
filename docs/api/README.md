# docs/api

这里是 API 文档总入口。对外导航只保留“从稳定入口到主题文档”的最短路径；当前状态以根 [README.md](../../README.md) 为准，脚本白名单见[稳定脚本指南](../guides/scripts/README.md)。

需要按完整计算流程查看“已实现能力、公式信源、代码入口、测试层级与验证边界”时，先读[已实现计算能力与方法追踪清单](../reference/implemented_capabilities.md)。

## 1. 先按你要做的事找入口

### 跑相图 / 看相结构工作流

- 运行入口：[`scripts/pnjl/calculate_phase_structure.jl`](../../scripts/pnjl/calculate_phase_structure.jl)
- 主题入口：[phase/README.md](models/phase/README.md)
- 首读页面：[phase/Overview.md](models/phase/Overview.md)

### 跑 T-μ / T-ρ / 磁场扫描

- 运行入口：[`scripts/models/run_unified_scan.jl`](../../scripts/models/run_unified_scan.jl)（子命令与参数见[命令目录](../dev/agent_command_reference.md)）
- 主题入口：[scans/README.md](models/scans/README.md)
- 首读页面：[scans/Overview.md](models/scans/Overview.md)

### 看 `Models` 统一编排与 workflow

- 主题入口：[workflows/README.md](models/workflows/README.md)
- 首读页面：[workflows/Overview.md](models/workflows/Overview.md)

### 看核心求解器与约束

- 主题入口：[solver/README.md](models/solver/README.md)
- 首读页面：[solver/Overview.md](models/solver/Overview.md)

### 看输运 / 弛豫时间链路

- 主题总入口：[relaxtime/README.md](relaxtime/README.md)
- 链路核心：[transport/README.md](relaxtime/transport/README.md)、[transport/CoreConcepts.md](relaxtime/transport/CoreConcepts.md)
- 分主题：[散射与截面](relaxtime/scattering/ScatteringAmplitude.md)、[极化](relaxtime/polarization/Polarization.md)、[传播子与有效耦合](relaxtime/propagator/EffectiveCouplings.md)、[介子热力学](relaxtime/meson_thermo/MesonThermodynamics.md)、[介子数密度](relaxtime/meson_density/MesonDensity.md)

### 看数值积分

- 主题入口：[integrals/OneLoopIntegrals.md](integrals/OneLoopIntegrals.md)（另含 `GaussLegendre`、`OneLoopIntegralsAniso`、`CauchyPV`、`CausalSpectralBubble`）

### 看模型变体或衍生量

- 模型变体总入口：[variants/README.md](models/variants/README.md)
- 衍生量总入口：[derived/README.md](models/derived/README.md)
- 兼容层说明：[pnjl/PNJL.md](pnjl/PNJL.md)

### 对接 Web / HTTP 服务

- 启动与使用：[稳定脚本指南](../guides/scripts/README.md)、[Web 示例](../guides/examples/web_demo.md)。
- 服务装配：[启动入口](../../scripts/server/server_full.jl) 调用 [ServerLauncher](../../src/simulation/ServerLauncher.jl)；[FullServerApp](../../src/simulation/FullServerApp.jl) 组织处理器与路由。
- 请求到计算的映射：[routing.jl](../../src/simulation/fullserver/routing.jl) 分发到同目录处理器；请求解析、单位转换和返回字段需对照相应处理器，异步任务另看 `pnjl_scan_jobs.jl` 或 `script_task_jobs.jl`。
- 计算层数据结构见[数据契约](data_contracts.md)及对应计算主题；Julia 函数可用不表示已有对应 HTTP 路由。稳定请求/响应说明放在对应接口页面，未实现的适配方案放在已有开发设计中。

### 使用公共绘图 helper

- [曲线身份、局部样式上下文与字体策略](plotting/PlotStyle.md)

### 查参数容器与类型

- [PARAMETER_TYPES_API.md](PARAMETER_TYPES_API.md)

## 2. 当前主导航口径

- `models/phase/`、`models/scans/`、`models/workflows/`、`models/solver/` 是 `Models` 主链的核心主题。
- `relaxtime/transport/` 是 transport 算法与 provider 的主主题；不再额外平铺为新的 `models/*` 一级主题。
- `pnjl/PNJL.md` 仅保留兼容层聚合说明；现行求解、扫描、相图入口以 `models/*` 为准。

## 3. 每个主题内部怎么读

按主题规模组织用户入口、职责核心与自动导出索引。小主题可在一个页面内呈现前两种视图；较大主题可以拆为：

1. `README.md`：主题入口与阅读顺序
2. `Overview.md`：面向调用方的用户入口
3. `Algorithms.md` 或 `CoreConcepts.md`：职责核心层
4. `generated/Exports.md`：自动生成的导出 API 全集

## 4. 维护约定

- 同一事实尽量只保留一个 authoritative 页面；根 README 负责“跑什么、从哪进”，API 文档负责“接口做什么、怎么组织”。
- 新 `Models` 能力优先吸收到对应主题，不再把旧 `pnjl/` 页面当并列主入口。
- 模型变体类与衍生量类主题优先进入分类总层，例如 `models/variants/`、`models/derived/`，避免 `models/` 一级目录继续横向膨胀。
- 本页所有入口写成可点击的相对链接；新增主题时同时在本页登记，避免出现“文档存在但索引不到”。
