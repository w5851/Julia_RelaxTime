---
title: Transport 主题重复说明收口与页面减载任务单
archived: true
original: docs/dev/active/2026-03-10_Transport主题重复说明收口与页面减载任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Transport 主题重复说明收口与页面减载任务单

更新日期：2026-03-10

> 目标：在 transport 主题边界、目录与 generated coverage 已完成的前提下，收口 `relaxtime/transport` 与 workflow 相关页面之间仍然存在的重复说明，降低单页负载，稳定“入口页 -> 职责核心页 -> 细节页”的阅读路径。

---

## 1. 背景

上一轮任务已经完成：

- transport 主归属固定到 `docs/api/relaxtime/transport/`；
- `Models` 侧只保留 workflow 入口与桥接；
- `CoreConcepts.md` 已作为 transport 主题的职责核心页落地；
- `generated/Exports.md` 与 `ModelsExportIndex.md` 中 provider 三件套 coverage 已闭合。

当前剩余问题不再是“主题归属”，而是“说明重复与页面过载”：

- `TransportCoefficients.md` 仍保留较长的 provider 相关说明；
- `docs/api/relaxtime/workflow/TransportWorkflow.md` 与 `docs/api/models/workflows/TransportWorkflow.md` 在输入分层与 workflow 说明上存在部分重叠；
- `RelaxationTime.md` 与 `AverageScatteringRate.md` 在默认行为、截断口径和缓存说明上仍有可压缩空间；
- 新的 transport 阅读路径虽然已建立，但“应该在哪一页停下继续读”还可以更明确。

---

## 2. 任务目标

- [x] 收缩 `TransportCoefficients.md` 中重复于 `CoreConcepts.md` 的 provider 说明。
- [x] 明确 `relaxtime/workflow/TransportWorkflow.md` 与 `models/workflows/TransportWorkflow.md` 的分工，减少叙述交叉。
- [x] 审查 `RelaxationTime.md` 与 `AverageScatteringRate.md` 中可合并或跳转替代的说明段落。
- [x] 保持 generated coverage 不回退。
- [x] 在不改动源码的前提下，降低 transport 主题的阅读成本。

---

## 3. 范围与非目标

### 3.1 范围

- `docs/api/relaxtime/transport/TransportCoefficients.md`
- `docs/api/relaxtime/workflow/TransportWorkflow.md`
- `docs/api/models/workflows/TransportWorkflow.md`
- `docs/api/relaxtime/transport/RelaxationTime.md`
- `docs/api/relaxtime/transport/AverageScatteringRate.md`

### 3.2 非目标

- 不改动 transport 相关 Julia 源码；
- 不重新设计 transport 主题目录；
- 不扩大到 scattering / propagator 其他子主题；
- 不重写 generated 页面。

---

## 4. 预期产出

- 更短、更聚焦的 `TransportCoefficients.md`
- 更清晰的 workflow 双页分工
- 更稳定的 transport 阅读顺序
- 保持 17/17 generated coverage 不回退

---

## 5. 建议执行顺序

1. 先精简 `TransportCoefficients.md` 的 provider 段落，只保留调用方必须知道的摘要与跳转。
2. 再处理两个 workflow 页面，明确一个讲业务入口，一个讲领域流程细节。
3. 最后处理 `RelaxationTime.md` 与 `AverageScatteringRate.md` 的重复默认行为说明。

实施进展（2026-03-10）：

- 已将 `TransportCoefficients.md` 中完整 provider 契约表收束为“provider 摘要 + 跳转到 `CoreConcepts.md`”；
- 已明确 `models/workflows/TransportWorkflow.md` 负责入口定位，而 `relaxtime/workflow/TransportWorkflow.md` 负责内部流程、输入合同与优先级；
- 已把 `RelaxationTime.md` 作为默认截断与 C++ 对齐的主说明页，并将 `AverageScatteringRate.md` 压缩为单过程积分与缓存细节页；
- 后续只需做回归验证与必要的微调，不再需要扩大 transport 文档边界。

结论：

- 本任务单目标已经完成，可作为已落地记录保留；
- 若继续推进 transport 文档，下一步不应再做结构治理，而应转向更小粒度的内容打磨，例如补更短示例、统一术语或增强 cross-link。