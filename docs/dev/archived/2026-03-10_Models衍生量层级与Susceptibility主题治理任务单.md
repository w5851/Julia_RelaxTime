---
title: Models 衍生量层级与 Susceptibility 主题治理任务单
archived: true
original: docs/dev/active/2026-03-10_Models衍生量层级与Susceptibility主题治理任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 衍生量层级与 Susceptibility 主题治理任务单

更新日期：2026-03-10

> 目标：在 `Models` 主题分类法已经确定的前提下，先为“衍生量”建立可复用的目录方案，并以 susceptibility 为首个落地对象，形成可执行的三层视图治理任务单。

---

## 1. 背景与问题定义

在 [2026-03-10_Models主题分类与目录分层治理任务单.md](2026-03-10_Models%E4%B8%BB%E9%A2%98%E5%88%86%E7%B1%BB%E4%B8%8E%E7%9B%AE%E5%BD%95%E5%88%86%E5%B1%82%E6%B2%BB%E7%90%86%E4%BB%BB%E5%8A%A1%E5%8D%95.md) 中，已经确认：

- susceptibility、cumulant、各类导数量不应直接作为 `docs/api/models/` 下的一级平铺主题；
- 这一簇能力应归入 `Models` 的“衍生量主题簇”；
- 当前建议先拆为 `susceptibility` 与 `derivatives` 两个子主题，而不是先做一个过大的总主题；
- `derived` 应作为与 `variants` 对应的分类总层。

当前这条主线仍处于明显的“导出存在、主题缺位”状态：

- `conserved_charge_susceptibility`、`chi_*`、`chi11_*`、`cumulant_*`、`baryon_*` 等导出在 [docs/api/generated/models/ModelsExportIndex.md](../../api/generated/models/ModelsExportIndex.md) 中大面积未命中；
- 旧文档资产主要集中在 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md)，更偏导数和体粘滞组合导数，而不是完整的 susceptibility / cumulant 主题；
- 这意味着如果不先建立 `derived` 总层与 susceptibility 子主题，就会继续出现“导出很多，但找不到主导航主题”的问题。

因此，本任务单的目标不是立即补一批零散说明，而是同时完成两件事：

1. 固定 `docs/api/models/derived/` 这一层的目录方案；
2. 以 susceptibility 为首个对象，建立后续三层视图文档治理的正式实施计划。

---

## 2. 任务目标

- [x] 固定 `docs/api/models/derived/` 作为 `Models` 衍生量主题总层的推荐路径。
- [x] 固定 susceptibility 的主题落点为 `docs/api/models/derived/susceptibility/`。
- [x] 设计 susceptibility 主题的三层视图结构与子页面拆分。
- [x] 明确旧 `docs/api/pnjl/` 导数页中哪些内容应被 susceptibility / derivatives 主题吸收。
- [x] 固定 susceptibility 主题的导出过滤集合与生成页策略。
- [x] 为后续真正落地正文页面提供可直接执行的任务顺序与 DoD。

任务目标落地说明（2026-03-10）：

- 已建立 [docs/api/models/derived/README.md](../../api/models/derived/README.md) 与 [docs/api/models/derived/susceptibility/README.md](../../api/models/derived/susceptibility/README.md)；
- 已建立 susceptibility 主题三层与细节页：
  - [docs/api/models/derived/susceptibility/Overview.md](../../api/models/derived/susceptibility/Overview.md)
  - [docs/api/models/derived/susceptibility/CoreConcepts.md](../../api/models/derived/susceptibility/CoreConcepts.md)
  - [docs/api/models/derived/susceptibility/ConservedChargeSusceptibilities.md](../../api/models/derived/susceptibility/ConservedChargeSusceptibilities.md)
  - [docs/api/models/derived/susceptibility/Cumulants.md](../../api/models/derived/susceptibility/Cumulants.md)
  - [docs/api/models/derived/susceptibility/generated/Exports.md](../../api/models/derived/susceptibility/generated/Exports.md)
- 总导航 [docs/api/README.md](../../api/README.md) 与脚本说明 [scripts/dev/README.md](../../../scripts/dev/README.md) 已同步更新；
- 旧页 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 已调整为 derivatives 过渡页，并增加 susceptibility 迁移指引；
- `Models` susceptibility 主题导出全集已生成，目标导出 24/24 均已在非 generated API 页面中命中。

---

## 3. 范围与非目标

### 3.1 范围

- `derived` 总层目录方案；
- susceptibility 主题的目标目录与页面结构；
- susceptibility 相关公开导出的过滤集合；
- 旧导数页与新主题之间的迁移/分流策略；
- 后续正文实施时的阶段拆分与验收标准。

### 3.2 非目标

- 不在本任务中直接落地 derivatives 子主题正文；
- 不在本任务中修改 `src/models/*` 实现代码；
- 不在本任务中处理 transport 相关主题。

---

## 4. 现状盘点

### 4.1 当前公开导出簇

根据 `src/models/Models.jl`，当前 susceptibility / cumulant / 导数相关公开导出至少包括：

- `conserved_charge_susceptibility`
- `chi_BQS`
- `chi_B`, `chi1_B`, `chi2_B`, `chi3_B`, `chi4_B`
- `chi_Q`, `chi1_Q`, `chi2_Q`, `chi3_Q`, `chi4_Q`
- `chi_S`, `chi1_S`, `chi2_S`, `chi3_S`, `chi4_S`
- `chi11_BQ`, `chi11_BS`, `chi11_QS`
- `cumulant_B`
- `cumulant_BQS`
- `baryon_Ssigma`
- `baryon_kappa_sigma2`
- `flavor_pressure_derivatives`
- `mass_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `thermo_derivatives`
- `compute_B_bracket`

### 4.2 当前旧文档资产

当前可直接复用的旧路径文档主要是：

- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md)

现有内容已覆盖：

- `thermo_derivatives`
- `dP_dT` / `dP_dmu`
- `bulk_derivative_coeffs`
- 组合导数与体粘滞使用场景

但仍存在明显缺口：

- conserved charge susceptibility 家族几乎没有系统性主题说明；
- `chi_*`、`chi11_*`、`cumulant_*` 与 `baryon_*` 系列尚未形成统一导航；
- 旧页更偏 PNJL 热力学偏导，而不是 `Models` 聚合导出的衍生量主题。

### 4.3 当前覆盖现状

根据 [docs/api/generated/models/ModelsExportIndex.md](../../api/generated/models/ModelsExportIndex.md)：

- `conserved_charge_susceptibility`、`chi_*`、`chi11_*`、`cumulant_*`、`baryon_*` 当前大量处于 0 命中；
- `thermo_derivatives`、`dP_dT`、`dP_dmu` 虽有旧页命中，但尚未进入 `Models` 新主题语境；
- 这表明 susceptibility 更适合作为 `derived` 总层下的首个子主题，而不是继续散落在旧导数页与 generated 索引之间。

---

## 5. 目录方案

### 5.1 derived 总层

推荐新增目录层：

- `docs/api/models/derived/`

其职责：

- 承接 `Models` 下不适合一级平铺、但又明确属于派生物理量或响应量族的主题；
- 避免继续把 susceptibility / derivatives / cumulants 误挂为 `docs/api/models/<theme>/` 一级目录；
- 为后续 response / fluctuation / derivative 相关主题预留统一层级。

### 5.2 susceptibility 主题落点

推荐主题目录：

- `docs/api/models/derived/susceptibility/`

说明：

- 该路径是当前推荐方案，但本任务单阶段只固定方案，不在本轮直接创建大量正文文件；
- 后续正文任务应以此路径为唯一主目录，而不是回退到 `docs/api/models/susceptibility/`。

---

## 6. susceptibility 主题结构设计

### 6.1 三层视图映射

| 视图 | 主页面 | 说明 |
| --- | --- | --- |
| 面向用户入口 | `derived/susceptibility/Overview.md` | 先说明 conserved charge susceptibility 家族、如何从 `Models` 入口理解 `chi_*` / `cumulant_*` / `baryon_*` |
| 职责核心 | `derived/susceptibility/CoreConcepts.md` | 解释 susceptibilities、cumulants、组合量与导数入口之间的边界 |
| 导出 API 全集 | `derived/susceptibility/generated/Exports.md` | 由脚本生成的 susceptibility 公开导出子集索引 |

### 6.2 建议页面结构

- `docs/api/models/derived/README.md`
- `docs/api/models/derived/susceptibility/README.md`
- `docs/api/models/derived/susceptibility/Overview.md`
- `docs/api/models/derived/susceptibility/CoreConcepts.md`
- `docs/api/models/derived/susceptibility/ConservedChargeSusceptibilities.md`
- `docs/api/models/derived/susceptibility/Cumulants.md`
- `docs/api/models/derived/susceptibility/generated/Exports.md`

### 6.3 页面职责建议

- `docs/api/models/derived/README.md`
  - `derived` 总层入口与后续主题索引
- `susceptibility/README.md`
  - susceptibility 主题主入口与阅读顺序
- `Overview.md`
  - 面向用户的 `Models` 衍生量入口说明，先讲族的结构而不是逐个函数清单
- `CoreConcepts.md`
  - susceptibilities、cumulants、`baryon_*` 组合量之间的关系
- `ConservedChargeSusceptibilities.md`
  - `conserved_charge_susceptibility`、`chi_BQS`、`chi_*`、`chi11_*`
- `Cumulants.md`
  - `cumulant_B`、`cumulant_BQS`、`baryon_Ssigma`、`baryon_kappa_sigma2`

说明：

- `thermo_derivatives`、`dP_dT`、`dP_dmu`、`bulk_derivative_coeffs` 等更适合作为后续 `derived/derivatives/` 子主题处理；
- 本 susceptibility 主题应优先处理 conserved charge / cumulant 主簇，而不是把所有导数能力混在一起。

---

## 7. 旧页面吸收与分流策略

### 7.1 需要吸收的内容

从 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 吸收：

- `thermo_derivatives` 的结果结构与组合导数语义
- `dP_dT` / `dP_dmu` 与组合导数的口径说明
- 体粘滞相关的导数耦合背景

### 7.2 分流原则

- 与 conserved charge susceptibility、cumulants、`baryon_*` 更直接相关的说明，应进入 `susceptibility` 子主题；
- 与 `thermo_derivatives`、`mass_derivatives`、`bulk_derivative_coeffs` 更直接相关的说明，应保留给后续 `derivatives` 子主题；
- 不建议把旧 `ThermoDerivatives.md` 整页简单搬迁到 susceptibility 主题中。

### 7.3 降级策略

后续当 `susceptibility` 与 `derivatives` 两个子主题都落地后：

- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 应降级为迁移说明或兼容说明页；
- 不再让旧页继续承担 `Models` 衍生量主题的主体说明职责。

---

## 8. 自动化策略

### 8.1 susceptibility 主题导出来源

- 首选来源：`src/models/Models.jl`

### 8.2 susceptibility 初始过滤集合

- `conserved_charge_susceptibility`
- `chi_BQS`
- `chi_B`, `chi1_B`, `chi2_B`, `chi3_B`, `chi4_B`
- `chi_Q`, `chi1_Q`, `chi2_Q`, `chi3_Q`, `chi4_Q`
- `chi_S`, `chi1_S`, `chi2_S`, `chi3_S`, `chi4_S`
- `chi11_BQ`, `chi11_BS`, `chi11_QS`
- `cumulant_B`
- `cumulant_BQS`
- `baryon_Ssigma`
- `baryon_kappa_sigma2`

### 8.3 暂不纳入 susceptibility 首轮过滤集合

以下导出更适合作为后续 `derivatives` 子主题处理：

- `flavor_pressure_derivatives`
- `mass_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `thermo_derivatives`
- `compute_B_bracket`

### 8.4 生成目标

- 输出位置：`docs/api/models/derived/susceptibility/generated/Exports.md`
- 仍使用 `scripts/dev/generate_api_export_index.jl`

---

## 9. 实施顺序

### 阶段 A：固定 derived 总层与 susceptibility 落点

- [x] 在任务层确认 `docs/api/models/derived/` 为推荐总层
- [x] 在任务层确认 `docs/api/models/derived/susceptibility/` 为唯一推荐落点
- [x] 避免后续再回到 `docs/api/models/susceptibility/` 的旧讨论

阶段 A 落地说明（2026-03-10）：

- [docs/api/models/derived/README.md](../../api/models/derived/README.md) 已建立为总层入口；
- [docs/api/README.md](../../api/README.md) 已将 `derived` 总层与 susceptibility 主题纳入主导航。

### 阶段 B：建立三层视图正文骨架

- [x] 新建 derived 总层入口与 susceptibility 主题目录
- [x] 新建 `README.md`
- [x] 新建 `Overview.md`
- [x] 新建 `CoreConcepts.md`
- [x] 建立 `generated/Exports.md`

阶段 B 落地说明（2026-03-10）：

- 已建立 `README.md`、`Overview.md`、`CoreConcepts.md` 以及两张细节页；
- 已通过 `scripts/dev/generate_api_export_index.jl` 生成 [docs/api/models/derived/susceptibility/generated/Exports.md](../../api/models/derived/susceptibility/generated/Exports.md)；
- 生成结果显示 24 个目标导出全部已在非 generated API 页面中命中。

### 阶段 C：吸收与分流旧页内容

- [x] 将旧 `ThermoDerivatives.md` 中与 susceptibility / cumulants 相关的高价值说明吸收入新主题
- [x] 明确哪些说明应留给后续 `derivatives` 子主题
- [x] 统一 `Models` 视角与旧 `PNJL.*` 兼容语义之间的边界

阶段 C 落地说明（2026-03-10）：

- 新主题已吸收守恒荷 susceptibility、cumulant 与 baryon 组合量的主说明，不再依赖旧 `PNJL.*` 页面做主体补充；
- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 已明确改为 derivatives 过渡页，并增加新主题跳转；
- `thermo_derivatives`、`mass_derivatives`、`bulk_derivative_coeffs`、`dP_dT`、`dP_dmu`、`compute_B_bracket` 已明确留给后续 `derived/derivatives/` 子主题。

### 阶段 D：旧页降级与导航收口

- [x] 后续在 susceptibility 与 derivatives 均落地后，再降级旧 `ThermoDerivatives.md`
- [x] `docs/api/README.md` 纳入 `derived` 总层导航
- [x] susceptibility 主题成为该能力簇的首个主入口

阶段 D 落地说明（2026-03-10）：

- [docs/api/README.md](../../api/README.md) 已纳入 `derived` 总层导航；
- susceptibility 主题当前唯一主入口为 [docs/api/models/derived/susceptibility/README.md](../../api/models/derived/susceptibility/README.md)；
- 旧 `ThermoDerivatives.md` 已从“主主题页”降为 derivatives 过渡页，待后续 derivatives 子主题落地后再进一步收口为迁移说明。

---

## 10. DoD

- [x] 已形成 `derived` 总层的正式目录方案
- [x] 已形成 susceptibility 主题的正式目标路径与页面结构
- [x] 已明确旧导数页应吸收与分流的范围
- [x] 已形成可直接执行的导出过滤集合
- [x] 后续正文实施可直接沿本任务单推进，而无需再次讨论目录落点

DoD 收口说明（2026-03-10）：

- `derived` 总层与 susceptibility 主题目录已经实际落地，而不是停留在方案层；
- susceptibility 主题的三层与细节页已经建立，自动生成的导出全集也已生成；
- 导出全集已确认目标导出 24/24 全部被新主题人工文档吸收；
- 旧 `ThermoDerivatives` 页面已退居 derivatives 过渡职责，不再承担守恒荷 susceptibility 主说明。

---

## 11. 后续直接衔接任务

本任务单完成后，下一步应直接进入：

1. `Models` derivatives 子主题任务单
2. `derived/derivatives/` 正文与生成页实施

说明：

- 完成 susceptibility 后，`derived` 总层就有了首个真实落地对象；
- 随后再做 derivatives 时，可以复用同一层级结构与旧页分流策略。