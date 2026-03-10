# docs/api

这里是项目的 API 文档入口。

## 组织方式

- `docs/api/models/phase/`：相图主题在新架构下的主目录。
	- 主题入口：`docs/api/models/phase/README.md`
	- 用户入口：`docs/api/models/phase/Overview.md`
	- 算法核心：`docs/api/models/phase/Algorithms.md`
	- 导出 API 全集：`docs/api/models/phase/generated/Exports.md`
- `docs/api/models/workflows/`：`Models` 统一入口下的 workflow 主题主目录。
	- 主题入口：`docs/api/models/workflows/README.md`
	- 用户入口：`docs/api/models/workflows/Overview.md`
	- 职责核心：`docs/api/models/workflows/CoreConcepts.md`
	- 导出 API 全集：`docs/api/models/workflows/generated/Exports.md`
- `docs/api/models/scans/`：`Models` 统一入口下的扫描主题主目录。
	- 主题入口：`docs/api/models/scans/README.md`
	- 用户入口：`docs/api/models/scans/Overview.md`
	- 职责核心：`docs/api/models/scans/Algorithms.md`
	- 导出 API 全集：`docs/api/models/scans/generated/Exports.md`
- `docs/api/models/solver/`：`Models` 统一入口下的核心求解与约束主题主目录。
	- 主题入口：`docs/api/models/solver/README.md`
	- 用户入口：`docs/api/models/solver/Overview.md`
	- 职责核心：`docs/api/models/solver/CoreConcepts.md`
	- 导出 API 全集：`docs/api/models/solver/generated/Exports.md`
- `docs/api/models/variants/`：`Models` 模型家族/模型变体主题总层。
	- 当前主题：`docs/api/models/variants/magnetic/`
	- 主题入口：`docs/api/models/variants/magnetic/README.md`
	- 用户入口：`docs/api/models/variants/magnetic/Overview.md`
	- 职责核心：`docs/api/models/variants/magnetic/CoreConcepts.md`
	- 导出 API 全集：`docs/api/models/variants/magnetic/generated/Exports.md`
- `docs/api/models/derived/`：`Models` 衍生量/响应量主题总层。
	- 当前主题：`docs/api/models/derived/susceptibility/`
	- 主题入口：`docs/api/models/derived/susceptibility/README.md`
	- 用户入口：`docs/api/models/derived/susceptibility/Overview.md`
	- 职责核心：`docs/api/models/derived/susceptibility/CoreConcepts.md`
	- 导出 API 全集：`docs/api/models/derived/susceptibility/generated/Exports.md`
	- 当前主题：`docs/api/models/derived/derivatives/`
	- 主题入口：`docs/api/models/derived/derivatives/README.md`
	- 用户入口：`docs/api/models/derived/derivatives/Overview.md`
	- 职责核心：`docs/api/models/derived/derivatives/CoreConcepts.md`
	- 导出 API 全集：`docs/api/models/derived/derivatives/generated/Exports.md`
- `docs/api/models/` 下后续主题不再默认全部平铺并列。
	- 模型家族或模型变体类主题，应优先放入带额外分类层的目录，而不是直接新增一级目录；magnetic 属于这一类。
	- susceptibility、cumulant、导数量等能力，优先按 `Models` 的“衍生量”主题簇治理，而不是与 phase/workflows/scans/solver 直接同层并列；当前已落地 `derived/susceptibility/` 与 `derived/derivatives/` 两条主子线。
	- transport 相关能力不作为新的 `Models` 一级主题平铺；统一业务入口继续留在 `models/workflows/`，而 transport 算法与 provider 主体归入 `docs/api/relaxtime/transport/`。
- `docs/api/relaxtime/*`：relaxtime 主链路相关的 API 文档（按功能拆分到 polarization/propagator/scattering/transport/workflow）。
- `docs/api/relaxtime/transport/`：transport 主题入口目录。
	- 主题入口：`docs/api/relaxtime/transport/README.md`
	- 用户入口：`docs/api/relaxtime/transport/Overview.md`
	- 职责核心：`docs/api/relaxtime/transport/CoreConcepts.md`
	- 导出 API 全集：`docs/api/relaxtime/transport/generated/Exports.md`
	- 细节页：`TransportCoefficients.md`、`RelaxationTime.md`、`AverageScatteringRate.md`
- `docs/api/pnjl/PNJL.md`：PNJL 兼容层聚合说明页；求解、扫描、相图、外磁场等现行主入口均已迁移到 `docs/api/models/*` 对应主题。
- `docs/api/integrals/*.md`：跨模块通用的数值积分/数学工具类 API 文档（与具体物理链路解耦）。
- `docs/api/relaxtime/`：relaxtime 领域主题总入口（见该目录下 README）。

## 命名约定

- 文档文件名尽量与 `src/` 中的 Julia 源文件同名（如 `TotalCrossSection.jl` → `TotalCrossSection.md`），便于从名字直接定位实现。
- 如果某 API 属于新架构 `models` 主题，优先放入对应主题子目录，例如 `docs/api/models/phase/`、`docs/api/models/workflows/`、`docs/api/models/scans/`、`docs/api/models/solver/`。
- 但并非所有 `Models` 主题都应直接作为 `docs/api/models/` 下的一级目录平铺；模型变体类与衍生量类主题应先判断是否需要额外分类层。
- 三层视图中的第二层统一按“职责核心层”理解；若主题偏算法判据，可使用 `Algorithms.md`，若主题偏职责边界与模块协作，可使用 `CoreConcepts.md`。
- 历史 `docs/api/pnjl/` 页面在完成迁移后应退出主导航，避免形成并列主入口。
- 新 `docs/api/models/` 主题应优先吸收旧路径页面中仍有价值的内容；不建议把“旧路径详见另页”长期保留为新主题的主要补充方式。

如果后续需要做目录重构，建议先完成新主题吸收与内部链接迁移，再统一移除薄兼容旧页，避免长期双入口并存。
