---
title: Transport 相关导出边界与主题归属治理任务单
archived: true
original: docs/dev/active/2026-03-10_Transport相关导出边界与主题归属治理任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Transport 相关导出边界与主题归属治理任务单

更新日期：2026-03-10

> 目标：在 `Models` 衍生量主题已经完成 susceptibility / derivatives 首轮治理之后，进一步梳理 transport 相关公开导出的主题边界、主次归属与目录方案，避免 `Models` 聚合入口、`relaxtime` 领域算法页与 derivatives 上游依赖继续交叉混写。

---

## 1. 背景与问题定义

在已完成的 `Models` 文档治理中：

- `docs/api/models/workflows/` 已经承接 `Models.solve_gap_and_transport` 与 `Models.solve_transport_from_equilibrium` 这类高层业务入口；
- `docs/api/models/derived/derivatives/` 已经承接 `bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p` 等体粘滞上游导数能力；
- `docs/api/relaxtime/transport/TransportCoefficients.md` 与 `docs/api/relaxtime/workflow/TransportWorkflow.md` 又分别承担了输运积分公式与 transport workflow 细节说明。

这导致 transport 主题出现三类职责混叠：

1. `Models` 对外聚合入口：用户从统一命名空间直接触发“平衡态 + 输运”工作流；
2. `relaxtime` 领域算法与 provider 契约：输运积分、弛豫时间、provider 行为与积分配置属于 transport 领域细节；
3. `derivatives` 上游依赖：体粘滞需要的热力学导数组合属于导数主题，不应反向吞并 transport 主主题。

当前最明显的问题是：

- `TransportProvider` 与 `prepare_transport_provider` 在 [docs/api/generated/models/ModelsExportIndex.md](../../api/generated/models/ModelsExportIndex.md) 中仍为 0 命中；
- `transport_provider` 虽有 1 次命中，但目前只是在 [docs/api/relaxtime/transport/TransportCoefficients.md](../../api/relaxtime/transport/TransportCoefficients.md) 中被顺带提及，尚未形成正式主题归属；
- `solve_gap_and_transport` 与 `solve_transport_from_equilibrium` 已经进入 `Models` workflow 文档，但其与 `relaxtime` 领域页之间的主次关系还没有被显式固定；
- `bulk_viscosity_coefficients` 与 `legacy_transport_c_p` 已经在 derivatives 中降级定位，但 transport 主题尚未明确“它们只是上游依赖，而不是 transport 首页主角”。

因此，本任务单不直接展开 transport 正文大规模改写，而是先完成：

1. transport 相关导出的主题边界梳理；
2. `Models` 与 `relaxtime` 的主次归属决策；
3. 后续目录与迁移顺序的初步方案固定。

---

## 2. 任务目标

- [x] 固定 transport 主题的主归属应偏向 `relaxtime` 领域，而不是在 `Models` 下新开一个并列一级 transport 主题。
- [x] 固定 `Models` 侧只保留统一工作流入口与桥接说明，不重复承担 transport 领域细节正文。
- [x] 明确 `transport_provider`、`TransportProvider`、`prepare_transport_provider` 的首选主题归属与次级入口位置。
- [x] 明确 `bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p` 与 transport 主题之间的边界。
- [x] 提出 transport 相关目录的初步扩展方案，为下一轮正文治理提供可执行蓝图。
- [x] 给出分阶段迁移建议与 DoD，限定本轮只做边界治理，不提前摊开全文实施。

任务目标落地说明（2026-03-10）：

- 已基于源码、现有 API 页与 generated index 梳理 transport 主题的三类职责边界；
- 已确认 transport 的领域主入口应放在 `docs/api/relaxtime/transport/` 及其相邻 workflow 页面；
- 已确认 `Models` 侧保留 `docs/api/models/workflows/TransportWorkflow.md` 作为业务入口页，不另设 `docs/api/models/transport/`；
- 已确认 `TransportProvider` 一簇虽然由 `Models` 导出，但语义上属于 transport 领域契约，应在 `relaxtime` transport 主题下取得正式主说明；
- 已固定后续目录与迁移阶段建议，供下一轮实施任务直接承接。

---

## 3. 范围与非目标

### 3.1 范围

- transport 相关 `Models` 公开导出的主题归属与边界；
- `docs/api/models/workflows/` 与 `docs/api/relaxtime/transport/` / `docs/api/relaxtime/workflow/` 的主次关系；
- `TransportProvider` 契约与 workflow 入口的分层方案；
- transport 与 derivatives 上游依赖的文档边界；
- 后续正文迁移的阶段计划与验收标准。

### 3.2 非目标

- 不在本任务中直接重写 [docs/api/relaxtime/transport/TransportCoefficients.md](../../api/relaxtime/transport/TransportCoefficients.md)；
- 不在本任务中直接新增 `docs/api/relaxtime/transport/README.md`、`Overview.md` 或 provider 细分页；
- 不在本任务中修改 `src/models/transport_provider.jl`、`src/models/workflows/TransportWorkflow.jl` 或 `src/relaxtime/*` 实现；
- 不在本任务中移动或改写 derivatives 已完成的页面；
- 不在本任务中直接闭合 `TransportProvider` / `prepare_transport_provider` 的 generated 覆盖缺口。

---

## 4. 现状盘点

### 4.1 当前 `Models` 相关公开导出

根据 `src/models/Models.jl`、`src/models/entrypoints.jl` 与 `src/models/transport_provider.jl`，transport 相关公开导出可分为两组：

#### A. `Models` 高层工作流入口

- `solve_gap_and_transport`
- `solve_transport_from_equilibrium`
- `transport_workflow_module`

这些导出的职责是：

- 暴露统一工作流入口；
- 将用户从 `Models` 命名空间转发到 transport workflow；
- 保持“先求平衡态，再输运”或“已知平衡态后处理”的业务视角。

#### B. `Models` 对外暴露的 transport provider 桥接契约

- `transport_provider`
- `TransportProvider`
- `prepare_transport_provider`

这些导出的职责是：

- 定义 models-side 可注入 provider；
- 把 `Models` 平衡态物理量组织成 transport 层可消费的最小契约；
- 作为 `Models` 与 `relaxtime` transport 积分器之间的耦合边界。

### 4.2 已有文档承接现状

当前 transport 主题相关文档主要分布为：

- `Models` 工作流入口页：
  - [docs/api/models/workflows/Overview.md](../../api/models/workflows/Overview.md)
  - [docs/api/models/workflows/TransportWorkflow.md](../../api/models/workflows/TransportWorkflow.md)
  - [docs/api/models/workflows/CoreConcepts.md](../../api/models/workflows/CoreConcepts.md)
  - [docs/api/models/workflows/ModuleAccessors.md](../../api/models/workflows/ModuleAccessors.md)
- `relaxtime` transport / workflow 领域页：
  - [docs/api/relaxtime/transport/TransportCoefficients.md](../../api/relaxtime/transport/TransportCoefficients.md)
  - [docs/api/relaxtime/transport/RelaxationTime.md](../../api/relaxtime/transport/RelaxationTime.md)
  - [docs/api/relaxtime/transport/AverageScatteringRate.md](../../api/relaxtime/transport/AverageScatteringRate.md)
  - [docs/api/relaxtime/workflow/TransportWorkflow.md](../../api/relaxtime/workflow/TransportWorkflow.md)
- `derived` 上游依赖页：
  - [docs/api/models/derived/derivatives/BulkViscosityDerivatives.md](../../api/models/derived/derivatives/BulkViscosityDerivatives.md)

### 4.3 generated 覆盖信号

根据 [docs/api/generated/models/ModelsExportIndex.md](../../api/generated/models/ModelsExportIndex.md)：

- `solve_gap_and_transport`：已有命中，主要落在 `models/workflows/*`；
- `solve_transport_from_equilibrium`：已有命中，主要落在 `models/workflows/*`；
- `transport_workflow_module`：已有命中，落在 `models/workflows/*`；
- `transport_provider`：仅 1 次命中，且主要是领域页中的顺带提及；
- `TransportProvider`：0 命中；
- `prepare_transport_provider`：0 命中。

这说明当前真正缺口不是 workflow 入口，而是 provider 契约在主题归属上尚未被正式承接。

### 4.4 实现边界事实

从源码可确认以下边界：

- [src/models/entrypoints.jl](../../../src/models/entrypoints.jl) 只做统一入口转发，不承载 transport 领域公式与 provider 细节；
- [src/models/workflows/TransportWorkflow.jl](../../../src/models/workflows/TransportWorkflow.jl) 负责串联平衡态、τ、bulk 与 transport 计算，是 workflow 协调层；
- [src/models/transport_provider.jl](../../../src/models/transport_provider.jl) 定义 `TransportProvider` 与 `prepare_transport_provider`，是 models-to-transport 契约边界；
- [docs/api/relaxtime/transport/TransportCoefficients.md](../../api/relaxtime/transport/TransportCoefficients.md) 已经详细说明 provider 字段契约、积分配置与 `prefer_energy_aniso` 行为；
- `bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p` 虽与 transport workflow 有调用关系，但它们的物理主题仍是 derivatives，而不是 transport provider / transport coefficient 本体。

---

## 5. 主题归属决策

### 5.1 总体原则

- transport 作为算法与领域主题，主归属固定在 `relaxtime`；
- `Models` 只保留“统一入口 / 聚合工作流 / 模块桥接”视角；
- 任何需要解释积分核、provider 契约、`prefer_energy_aniso` 优先级、积分配置的正文，应优先放在 `relaxtime`；
- 任何需要解释 `bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p` 的正文，应继续放在 derivatives，不回流到 transport 主页。

### 5.2 主次归属表

| 能力 / 导出 | 主归属 | 次级入口 | 说明 |
| --- | --- | --- | --- |
| `solve_gap_and_transport` | `docs/api/models/workflows/` | `docs/api/relaxtime/workflow/` | 从用户视角是 `Models` 聚合入口；领域 workflow 页负责细节 |
| `solve_transport_from_equilibrium` | `docs/api/models/workflows/` | `docs/api/relaxtime/workflow/` | 同上，属于后处理入口 |
| `transport_workflow_module` | `docs/api/models/workflows/ModuleAccessors.md` | 无需主领域扩展 | 属于模块访问器，不是 transport 领域首页能力 |
| `transport_provider` | `docs/api/relaxtime/transport/` | `docs/api/models/workflows/CoreConcepts.md` 或等效桥接页 | 虽由 `Models` 导出，但语义上是 transport provider 工厂 |
| `TransportProvider` | `docs/api/relaxtime/transport/` | `docs/api/models/workflows/CoreConcepts.md` 或等效桥接页 | 属于 provider 契约，不应只挂在 `Models` 导出表 |
| `prepare_transport_provider` | `docs/api/relaxtime/transport/` | `docs/api/models/workflows/CoreConcepts.md` 或等效桥接页 | 属于 provider 预备/缓存步骤，应与契约页同主题 |
| `bulk_viscosity_coefficients` | `docs/api/models/derived/derivatives/` | transport 页仅可链接引用 | 它是上游导数组合，不是 transport 主题本体 |
| `compute_B_bracket` | `docs/api/models/derived/derivatives/` | transport 页仅可链接引用 | 同上 |
| `legacy_transport_c_p` | `docs/api/models/derived/derivatives/` | transport 页不再升格 | 兼容导出，保持降级 |

### 5.3 不建议的方案

当前不建议：

- 在 `docs/api/models/` 下新开一级 `transport/` 主题；
- 把 provider 契约完全塞进 `docs/api/models/workflows/TransportWorkflow.md`；
- 为了闭合 `Models` exported coverage，把 `bulk_viscosity_coefficients` 重新并入 transport 主题；
- 让 `docs/api/relaxtime/workflow/TransportWorkflow.md` 继续承担 `Models` 用户入口的主导航角色。

这些方案都会重新制造边界混乱。

---

## 6. 初步目录方案

### 6.1 `relaxtime` 侧作为主领域主题

后续推荐把 transport 主主题固定为：

- `docs/api/relaxtime/transport/`

在不破坏现有页面的前提下，建议下一轮补齐：

- `docs/api/relaxtime/transport/README.md`
- `docs/api/relaxtime/transport/Overview.md`
- `docs/api/relaxtime/transport/CoreConcepts.md`

并保留现有细节页：

- `TransportCoefficients.md`
- `RelaxationTime.md`
- `AverageScatteringRate.md`

建议职责分工：

- `README.md`
  - transport 主题入口、阅读顺序、与 workflow / scattering / derivatives 的关系
- `Overview.md`
  - 面向用户的 transport 主题导航，先解释“算 transport 需要哪些层”
- `CoreConcepts.md`
  - 以“Provider Contract”为核心章节，正式承接 `transport_provider`、`TransportProvider`、`prepare_transport_provider`
- `TransportCoefficients.md`
  - 保持公式、积分核、`TransportIntegrationConfig`、provider 字段细节

### 6.2 `Models` 侧保留统一入口和桥接

后续 `Models` 侧不新建一级 transport 目录，继续以：

- `docs/api/models/workflows/TransportWorkflow.md`
- `docs/api/models/workflows/CoreConcepts.md`
- `docs/api/models/workflows/ModuleAccessors.md`

作为承接面。

建议后续微调方向：

- `TransportWorkflow.md`
  - 继续只讲 `solve_gap_and_transport` / `solve_transport_from_equilibrium` 的业务入口定位
- `CoreConcepts.md`
  - 增加一小节解释 provider 是“桥接契约”，正文主说明跳转到 `relaxtime/transport`
- `ModuleAccessors.md`
  - 继续承接 `transport_workflow_module` 等 accessor，不扩大到 provider 细节

如果后续发现 `CoreConcepts.md` 过载，再考虑新增：

- `docs/api/models/workflows/TransportBoundary.md`

但当前阶段不建议预先扩页。

### 6.3 derivatives 侧维持上游依赖主题

以下内容保持在已有 derivatives 主题：

- `bulk_viscosity_coefficients`
- `compute_B_bracket`
- `legacy_transport_c_p`

transport 正文只应：

- 说明这些量在 workflow 中被消费；
- 链接到 derivatives 页；
- 不再重复解释其公式与返回结构。

---

## 7. 后续导出覆盖策略

### 7.1 当前 coverage 状态

本轮原本只要求先固定边界，但在补齐 `relaxtime/transport` 主入口页与 provider 契约页之后，已经进一步完成了 coverage 验证：

- 已生成 `docs/api/relaxtime/transport/generated/Exports.md`；
- 已刷新 `docs/api/generated/models/ModelsExportIndex.md`；
- `transport_provider`、`TransportProvider`、`prepare_transport_provider` 已不再处于 0 命中状态。

### 7.2 下一轮实施时的优先闭合对象

下一轮正文治理应优先闭合：

- `transport_provider`
- `TransportProvider`
- `prepare_transport_provider`

闭合方式应优先通过 `relaxtime/transport` 主主题页面完成，而不是为了 coverage 再写一套 `Models` 平行说明。

### 7.3 generated 页策略结论

现已确认生成脚本支持 transport 主题的混合来源：

- `Models` 来源：`transport_provider`、`TransportProvider`、`prepare_transport_provider`
- `relaxtime` 来源：transport coefficient / relaxation time 相关公开接口

当前已采用以下方式落地：

1. 用人工维护页在 `relaxtime/transport` 闭合 provider 三件套的主说明；
2. 用多 `--module-file` 生成 `docs/api/relaxtime/transport/generated/Exports.md`；
3. 刷新 `ModelsExportIndex.md`，确认 provider 三件套在总索引中同步命中。

---

## 8. 分阶段实施建议

### 阶段 A：主入口页与目录补齐

- 在 `docs/api/relaxtime/transport/` 下补 `README.md` 与 `Overview.md`；
- 固定 transport 主题阅读顺序；
- 在 `docs/api/relaxtime/README.md` 中加入更清晰的 transport 导航。

实施进展（2026-03-10）：

- 已新增 `docs/api/relaxtime/transport/README.md`；
- 已新增 `docs/api/relaxtime/transport/Overview.md`；
- 已更新 `docs/api/relaxtime/README.md` 与 `docs/api/README.md`，把 transport 主主题导航接入总入口。

### 阶段 B：provider 契约落页

- 新增 `CoreConcepts.md`；
- 正式承接 `transport_provider`、`TransportProvider`、`prepare_transport_provider`；
- 从 `TransportCoefficients.md` 中抽离过长的 provider 契约说明，避免单页过载。

实施进展（2026-03-10）：

- 已新增 `docs/api/relaxtime/transport/CoreConcepts.md`；
- 已在 `docs/api/relaxtime/transport/TransportCoefficients.md` 增加到 `CoreConcepts.md` 的主跳转说明；
- provider 三件套的人工主说明已经有正式落点，但 generated coverage 仍待下一阶段验证。

### 阶段 C：`Models` 侧桥接收口

- 在 `docs/api/models/workflows/TransportWorkflow.md` 明确标注“领域细节请转 `relaxtime/transport`”；
- 在 `CoreConcepts.md` 补充 provider 桥接边界；
- 必要时补一个轻量桥接页，但不建立 `models/transport/` 一级主题。

实施进展（2026-03-10）：

- 已更新 `docs/api/models/workflows/TransportWorkflow.md`，加入 provider 契约跳转；
- 已更新 `docs/api/models/workflows/CoreConcepts.md`，明确 provider 只在 `Models` 侧承担桥接语义；
- 当前无需新增 `TransportBoundary.md`。

### 阶段 D：覆盖验证与旧说明收口

- 验证 `ModelsExportIndex.md` 中 provider 三件套是否已命中；
- 检查 `TransportCoefficients.md` 与 workflow 页是否仍有重复段落；
- 只保留“入口说明 -> 领域细节页”的单向导航，不再双向抢主说明权。

实施计划细化（2026-03-10）：

- 先用 `scripts/dev/generate_api_export_index.jl` 的多 `--module-file` 能力生成 `docs/api/relaxtime/transport/generated/Exports.md`；
- 用该 generated 页与 `ModelsExportIndex.md` 双重检查 provider 三件套是否已被人工页面承接；
- 若 coverage 已闭合，则阶段 D 可转入“重复说明收口”；若仍未闭合，再定点补 mention，而不是扩大主题边界。

实施进展（2026-03-10）：

- 已生成 `docs/api/relaxtime/transport/generated/Exports.md`；
- generated 页统计结果为 17/17 命中；
- 已刷新 `docs/api/generated/models/ModelsExportIndex.md`，其中：
  - `transport_provider` 为 6 命中；
  - `TransportProvider` 为 4 命中；
  - `prepare_transport_provider` 为 6 命中；
- 阶段 D 已进入“重复说明收口”而非“coverage 补洞”状态。

---

## 9. DoD

本轮边界治理完成的判定标准：

- [x] 已明确 transport 主题主归属为 `relaxtime`；
- [x] 已明确 `Models` 侧仅保留 workflow 入口与桥接角色；
- [x] 已明确 provider 三件套的首选主题归属；
- [x] 已明确 derivatives 与 transport 的边界；
- [x] 已给出后续目录与迁移阶段建议；
- [x] 已限制本轮不提前展开正文重写，防止范围失控。

下一轮正文实施的额外 DoD 预告：

- [x] `relaxtime/transport` 主入口页补齐；
- [x] provider 契约页面落地；
- [x] `Models` workflow 页完成桥接收口；
- [x] provider 三件套在 generated index 中实现有效命中。

本轮完成说明（2026-03-10）：

- 边界治理、主次归属、主题目录、generated coverage 与桥接收口均已完成；
- `relaxtime/transport` 职责核心页文件名已统一为 `CoreConcepts.md`，并将 “Provider Contract” 保留为页内主章节而不是文件名；
- 本任务单可以视为完成，不再作为后续去重与减载工作的承载文档。

---

## 10. 结论

本轮结论可以压缩为三句话：

- transport 的“算法与契约主题”归 `relaxtime`；
- transport 的“统一业务入口”归 `Models.workflows`；
- transport 所需的体粘滞导数组合继续归 `Models.derived.derivatives`。

后续实施应顺着这个主线推进，而不是再把三类内容折叠回同一页里。