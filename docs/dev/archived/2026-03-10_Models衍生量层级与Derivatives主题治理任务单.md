---
title: Models 衍生量层级与 Derivatives 主题治理任务单
archived: true
original: docs/dev/active/2026-03-10_Models衍生量层级与Derivatives主题治理任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 衍生量层级与 Derivatives 主题治理任务单

更新日期：2026-03-10

> 目标：在 `Models` 主题分类法与 `derived` 总层已经确定、susceptibility 主题已完成首轮落地的前提下，为 derivatives 建立正式子主题任务单，承接热力学导数、质量导数与体粘滞组合导数相关公开能力。

---

## 1. 背景与问题定义

在 [2026-03-10_Models主题分类与目录分层治理任务单.md](2026-03-10_Models%E4%B8%BB%E9%A2%98%E5%88%86%E7%B1%BB%E4%B8%8E%E7%9B%AE%E5%BD%95%E5%88%86%E5%B1%82%E6%B2%BB%E7%90%86%E4%BB%BB%E5%8A%A1%E5%8D%95.md) 中，已经确认：

- susceptibility、cumulant、各类导数量不应直接作为 `docs/api/models/` 下的一级平铺主题；
- 这一簇能力应统一归入 `Models` 的 `derived` 总层；
- `susceptibility` 与 `derivatives` 应拆成两个子主题，而不是继续共用旧 `PNJL.*` 导数页。

当前，首个 `derived` 子主题 susceptibility 已经完成落地，导出覆盖 24/24 已闭合；但 derivatives 主线仍停留在旧路径：

- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md)

现状问题：

- 旧页虽然已经改成 derivatives bridge，但主说明仍未迁移到 `docs/api/models/derived/` 新架构；
- `mass_derivatives`、`thermo_derivatives`、`bulk_derivative_coeffs`、`dP_dT`、`dP_dmu`、`bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p` 仍未形成统一主题导航；
- [docs/api/generated/models/ModelsExportIndex.md](../../api/generated/models/ModelsExportIndex.md) 中，`mass_derivatives`、`compute_B_bracket`、`legacy_transport_c_p` 仍处于未命中状态，`dP_dT`、`dP_dmu`、`thermo_derivatives` 也仍主要只由旧页承接；
- derivatives 与 susceptibility 的边界已经在前一任务中明确，但还没有转化为正式的目录方案和实施任务。

因此，本任务单的目标不是再补一张旧页，而是同时完成两件事：

1. 固定 derivatives 在 `derived` 总层中的正式落点；
2. 为 derivatives 主题后续正文实施建立可直接执行的三层视图任务单。

---

## 2. 任务目标

- [x] 固定 derivatives 的主题落点为 `docs/api/models/derived/derivatives/`。
- [x] 设计 derivatives 主题的三层视图结构与子页面拆分。
- [x] 明确旧 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 中哪些内容应被新主题直接吸收。
- [x] 固定 derivatives 主题的导出过滤集合与生成页策略。
- [x] 明确 `legacy_transport_c_p` 的兼容层定位，不把它误写成首选用户入口。
- [x] 为后续正文实施提供可直接执行的阶段拆分与 DoD。

任务目标落地说明（2026-03-10）：

- 已建立 [docs/api/models/derived/derivatives/README.md](../../api/models/derived/derivatives/README.md)；
- 已建立 derivatives 主题三层与细节页：
  - [docs/api/models/derived/derivatives/Overview.md](../../api/models/derived/derivatives/Overview.md)
  - [docs/api/models/derived/derivatives/CoreConcepts.md](../../api/models/derived/derivatives/CoreConcepts.md)
  - [docs/api/models/derived/derivatives/MassDerivatives.md](../../api/models/derived/derivatives/MassDerivatives.md)
  - [docs/api/models/derived/derivatives/ThermoDerivatives.md](../../api/models/derived/derivatives/ThermoDerivatives.md)
  - [docs/api/models/derived/derivatives/BulkViscosityDerivatives.md](../../api/models/derived/derivatives/BulkViscosityDerivatives.md)
  - [docs/api/models/derived/derivatives/generated/Exports.md](../../api/models/derived/derivatives/generated/Exports.md)
- 总导航 [docs/api/README.md](../../api/README.md)、derived 总层入口 [docs/api/models/derived/README.md](../../api/models/derived/README.md) 与脚本说明 [scripts/dev/README.md](../../../scripts/dev/README.md) 已同步更新；
- 旧页 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 已进一步调整为迁移说明/兼容桥接页；
- `Models` derivatives 主题导出全集已生成，目标导出 8/8 均已在非 generated API 页面中命中。

---

## 3. 范围与非目标

### 3.1 范围

- `docs/api/models/derived/derivatives/` 的目录方案；
- derivatives 主题的页面结构；
- derivatives 相关公开导出的过滤集合；
- 旧 `ThermoDerivatives.md` 与新主题之间的迁移/降级策略；
- 与测试、baseline、integration 使用场景相匹配的主题边界说明。

### 3.2 非目标

- 不在本任务中修改 `src/models/derivatives/ThermoDerivatives.jl` 实现；
- 不在本任务中改动 `tests/` 或 baseline；
- 不在本任务中处理 transport provider 总体目录问题；
- 不在本任务中继续扩展 susceptibility 主题。

---

## 4. 现状盘点

### 4.1 当前公开导出簇

根据 `src/models/Models.jl` 与 `src/models/derivatives/ThermoDerivatives.jl`，当前 derivatives 主题的 `Models` 公开导出至少包括：

- `mass_derivatives`
- `thermo_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `bulk_viscosity_coefficients`
- `compute_B_bracket`
- `legacy_transport_c_p`

说明：

- `solve_equilibrium_mu`、`quasiparticle_energy`、`dE_dT`、`dE_dmu` 虽在旧页中出现，但当前不属于 `Models` 聚合导出表面，不应在新主题中被写成主公开 API。

### 4.2 当前旧文档资产

当前可直接复用的旧路径文档主要是：

- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md)

现有内容已覆盖：

- `thermo_derivatives`
- `dP_dT` / `dP_dmu`
- `bulk_derivative_coeffs`
- 体粘滞组合导数的使用场景

当前缺口：

- `mass_derivatives` 仍未形成正式主题说明；
- `compute_B_bracket` 与 `bulk_viscosity_coefficients` 尚未进入 `Models` 新主题导航；
- `legacy_transport_c_p` 仍没有被明确标记为兼容层导出；
- 旧页中仍混有非 `Models` 导出接口，容易模糊主题边界。

### 4.3 当前覆盖现状

根据 [docs/api/generated/models/ModelsExportIndex.md](../../api/generated/models/ModelsExportIndex.md)：

- `mass_derivatives`：0 命中
- `compute_B_bracket`：0 命中
- `legacy_transport_c_p`：0 命中
- `thermo_derivatives`：当前仅由旧 `ThermoDerivatives.md` 命中
- `dP_dT` / `dP_dmu`：当前仅由旧 `ThermoDerivatives.md` 命中
- `bulk_derivative_coeffs`：当前主要由旧页与 relaxtime workflow 页承接

这表明 derivatives 主题是当前 `derived` 总层下最自然的下一步治理对象。

### 4.4 当前实现与验证事实源

与后续 derivatives 文档应对齐的真实路径包括：

- 源码：
  - `src/models/derivatives/ThermoDerivatives.jl`
- 测试：
  - `tests/unit/pnjl/test_thermo_derivatives.jl`
  - `tests/integration/models/test_models_derivatives_dual_smoke.jl`
  - `tests/integration/relaxtime/test_bulk_viscosity_derivatives.jl`
  - `tests/regression/pnjl/test_pnjl_thermo_derivatives_regression.jl`
- baseline：
  - `tests/baselines/pnjl/baseline_pnjl_thermo_derivatives_fixedpoints_v1.csv`

---

## 5. 目录方案

### 5.1 derivatives 主题落点

推荐主题目录：

- `docs/api/models/derived/derivatives/`

说明：

- 该路径应作为 derivatives 子主题的唯一正式落点；
- 后续正文任务不应再回退到 `docs/api/models/derivatives/` 或继续扩写旧 `docs/api/pnjl/ThermoDerivatives.md`。

### 5.2 与 susceptibility 的边界

当前边界应保持为：

- `derived/susceptibility/`：守恒荷 susceptibility、cumulant、`baryon_*` 组合量
- `derived/derivatives/`：热力学导数、质量导数、体粘滞组合导数与 B bracket 主线

不建议在 derivatives 主题中重新承接 `chi_*`、`cumulant_*` 等说明。

---

## 6. derivatives 主题结构设计

### 6.1 三层视图映射

| 视图 | 主页面 | 说明 |
| --- | --- | --- |
| 面向用户入口 | `derived/derivatives/Overview.md` | 先说明 `thermo_derivatives`、`mass_derivatives`、`bulk_derivative_coeffs` 的典型使用顺序 |
| 职责核心 | `derived/derivatives/CoreConcepts.md` | 解释隐函数求导、链式法则、bulk viscosity 组合导数与兼容层边界 |
| 导出 API 全集 | `derived/derivatives/generated/Exports.md` | 由脚本生成的 derivatives 公开导出子集索引 |

### 6.2 建议页面结构

- `docs/api/models/derived/derivatives/README.md`
- `docs/api/models/derived/derivatives/Overview.md`
- `docs/api/models/derived/derivatives/CoreConcepts.md`
- `docs/api/models/derived/derivatives/MassDerivatives.md`
- `docs/api/models/derived/derivatives/ThermoDerivatives.md`
- `docs/api/models/derived/derivatives/BulkViscosityDerivatives.md`
- `docs/api/models/derived/derivatives/generated/Exports.md`

### 6.3 页面职责建议

- `README.md`
  - 主题主入口与阅读顺序
- `Overview.md`
  - 首选入口、最短使用路径、单位口径与推荐使用顺序
- `CoreConcepts.md`
  - 隐函数求导主线、AD 策略、链式法则、compat 边界
- `MassDerivatives.md`
  - `mass_derivatives` 与质量导数结果结构
- `ThermoDerivatives.md`
  - `thermo_derivatives`、`dP_dT`、`dP_dmu`、`bulk_derivative_coeffs`
- `BulkViscosityDerivatives.md`
  - `bulk_viscosity_coefficients`、`compute_B_bracket`、`legacy_transport_c_p`

说明：

- `legacy_transport_c_p` 应明确归类为兼容导出，而不是首页首选入口；
- 若后续页内篇幅过大，可再细分出 `Compatibility.md`，但当前阶段不必先行拆更多页。

---

## 7. 旧页面吸收与降级策略

### 7.1 需要吸收的内容

从 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 吸收：

- `thermo_derivatives` 的结果结构与一阶导数语义
- `dP_dT` / `dP_dmu` 的高阶导数接口口径
- `bulk_derivative_coeffs` 的业务场景
- 与体粘滞组合导数相关的最短示例和注意事项

### 7.2 需要重新组织而不是原样搬迁的内容

- 旧页中的 `solve_equilibrium_mu`、`quasiparticle_energy`、`dE_dT`、`dE_dmu` 说明，不应未经筛选直接进入 `Models` 新主题主页面；
- 若这些能力未来仍需保留，应明确标注其是否属于旧 `PNJL.*` 路径或非 `Models` 聚合导出接口。

### 7.3 降级策略

后续当 `derived/derivatives/` 正文落地后：

- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 应进一步降级为迁移说明或兼容说明页；
- 不再让旧页继续承担 derivatives 主题的主体说明职责。

---

## 8. 自动化策略

### 8.1 derivatives 主题导出来源

- 首选来源：`src/models/Models.jl`

### 8.2 derivatives 初始过滤集合

- `mass_derivatives`
- `thermo_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `bulk_viscosity_coefficients`
- `compute_B_bracket`
- `legacy_transport_c_p`

### 8.3 兼容层处理原则

- `legacy_transport_c_p` 应保留在 generated 导出全集中；
- 但在人工维护页中，应显式标注为 compatibility / legacy 入口，不与 `thermo_derivatives`、`bulk_viscosity_coefficients` 并列成首选入口。

### 8.4 生成目标

- 输出位置：`docs/api/models/derived/derivatives/generated/Exports.md`
- 仍使用 `scripts/dev/generate_api_export_index.jl`

---

## 9. 实施顺序

### 阶段 A：固定 derivatives 落点与边界

- [x] 在任务层确认 `docs/api/models/derived/derivatives/` 为唯一推荐落点
- [x] 明确 derivatives 与 susceptibility 的主题边界
- [x] 避免后续再回到旧 `PNJL.*` 导数页扩写的路径

阶段 A 落地说明（2026-03-10）：

- [docs/api/models/derived/README.md](../../api/models/derived/README.md) 已将 derivatives 纳入 `derived` 总层；
- [docs/api/README.md](../../api/README.md) 已将 `derived/derivatives/` 纳入主导航；
- derivatives 与 susceptibility 的边界已分别落到两套主题页中。

### 阶段 B：建立三层视图正文骨架

- [x] 新建 derivatives 主题目录与 `README.md`
- [x] 新建 `Overview.md`
- [x] 新建 `CoreConcepts.md`
- [x] 新建细节页骨架
- [x] 建立 `generated/Exports.md`

阶段 B 落地说明（2026-03-10）：

- 已建立 `README.md`、`Overview.md`、`CoreConcepts.md` 与三张细节页；
- 已通过 `scripts/dev/generate_api_export_index.jl` 生成 [docs/api/models/derived/derivatives/generated/Exports.md](../../api/models/derived/derivatives/generated/Exports.md)；
- 生成结果显示 8 个目标导出全部已在非 generated API 页面中命中。

### 阶段 C：吸收旧页内容并收紧公开表面

- [x] 将旧 `ThermoDerivatives.md` 中与 `Models` 导出一致的高价值说明吸收入新主题
- [x] 明确哪些说明属于非 `Models` 聚合导出接口，不进入首轮主说明
- [x] 显式标记 `legacy_transport_c_p` 的兼容层定位

阶段 C 落地说明（2026-03-10）：

- 新主题已吸收 `thermo_derivatives`、`mass_derivatives`、`bulk_derivative_coeffs`、`dP_dT`、`dP_dmu`、`bulk_viscosity_coefficients`、`compute_B_bracket` 的主说明；
- 旧页中 `solve_equilibrium_mu`、`quasiparticle_energy`、`dE_dT`、`dE_dmu` 等非 `Models` 聚合导出接口已在任务边界中显式排除，不进入首轮主说明；
- [docs/api/models/derived/derivatives/BulkViscosityDerivatives.md](../../api/models/derived/derivatives/BulkViscosityDerivatives.md) 已明确将 `legacy_transport_c_p` 标为 compatibility / legacy 导出。

### 阶段 D：旧页降级与导航收口

- [x] 在 derivatives 主题落地后，将旧 `ThermoDerivatives.md` 进一步降级为迁移说明
- [x] 在 `docs/api/README.md` 中继续沿 `derived` 总层承接 derivatives 导航
- [x] 让 derivatives 成为该导数能力簇的唯一主入口

阶段 D 落地说明（2026-03-10）：

- [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 已明确改为旧路径兼容说明与迁移桥接页；
- derivatives 当前唯一主入口为 [docs/api/models/derived/derivatives/README.md](../../api/models/derived/derivatives/README.md)；
- API 总导航与 derived 总层入口均已收口到新的主题结构。

---

## 10. DoD

- [x] 已形成 derivatives 主题的正式目标路径与页面结构
- [x] 已明确旧导数页应吸收与排除的范围
- [x] 已形成可直接执行的 derivatives 导出过滤集合
- [x] 已明确 `legacy_transport_c_p` 的兼容层定位
- [x] 后续正文实施可直接沿本任务单推进，而无需再次讨论目录落点或主题边界

DoD 收口说明（2026-03-10）：

- `derived/derivatives/` 主题目录已经实际落地，而不是停留在方案层；
- derivatives 主题的三层与细节页已经建立，自动生成的导出全集也已生成；
- 导出全集已确认目标导出 8/8 全部被新主题人工文档吸收；
- 旧 `ThermoDerivatives` 页面已退居兼容桥接职责，不再承担 `Models` derivatives 主题主说明。

---

## 11. 后续直接衔接任务

本任务单完成后，下一步应直接进入：

1. transport 相关导出边界梳理与主题归属任务
2. 旧 [docs/api/pnjl/ThermoDerivatives.md](../../api/pnjl/ThermoDerivatives.md) 的最终简化收口

说明：

- derivatives 完成后，`derived` 总层的两条主子线将同时具备真实落地对象；
- 之后再处理 transport 相关主题时，可以在 susceptibility / derivatives 已收口的前提下单独判断其边界。