---
title: Models 变体层级与 Magnetic 主题治理任务单
archived: true
original: docs/dev/active/2026-03-10_Models变体层级与Magnetic主题治理任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 变体层级与 Magnetic 主题治理任务单

更新日期：2026-03-10

> 目标：在 `Models` 主题分类法已经确定的前提下，先为“模型变体”建立可复用的目录方案，并以 magnetic 为首个落地对象，形成可执行的三层视图治理任务单。

---

## 1. 背景与问题定义

在 [2026-03-10_Models主题分类与目录分层治理任务单.md](2026-03-10_Models%E4%B8%BB%E9%A2%98%E5%88%86%E7%B1%BB%E4%B8%8E%E7%9B%AE%E5%BD%95%E5%88%86%E5%B1%82%E6%B2%BB%E7%90%86%E4%BB%BB%E5%8A%A1%E5%8D%95.md) 中，已经确认：

- `phase`、`workflows`、`scans`、`solver` 继续作为 `docs/api/models/` 下的一级平铺主题；
- magnetic 不应直接进入 `docs/api/models/magnetic/`；
- magnetic 应作为“模型变体主题”的首个落地对象；
- 后续需要先建立 `variants` 总层，再进入正文实施。

当前 magnetic 相关公开说明仍主要停留在旧路径：

- [docs/api/pnjl/MagneticPNJL.md](docs/api/pnjl/MagneticPNJL.md)
- [docs/api/pnjl/MagneticPNJLParameters.md](docs/api/pnjl/MagneticPNJLParameters.md)

现状问题：

- 旧页面仍以 `PNJL.*` 兼容语义为主，而不是从 `Models` 聚合导出视角组织；
- 旧页中的源码路径说明仍带有历史路径口径，尚未收敛到当前 `src/models/` 事实源；
- magnetic 相关导出是一整个高相关能力簇，但目前没有在 `docs/api/models/` 下形成主题导航；
- 如果直接新建 `docs/api/models/magnetic/`，会破坏刚刚确定的 variants 分层原则。

因此，本任务单的目标不是立即“补一页 magnetic 文档”，而是同时完成两件事：

1. 固定 `docs/api/models/variants/` 这一层的目录方案；
2. 以 magnetic 为首个对象，建立后续三层视图文档治理的正式实施计划。

---

## 2. 任务目标

- [x] 固定 `docs/api/models/variants/` 作为 `Models` 模型变体主题总层的推荐路径。
- [x] 固定 magnetic 的主题落点为 `docs/api/models/variants/magnetic/`。
- [x] 设计 magnetic 主题的三层视图结构与子页面拆分。
- [x] 明确旧 `docs/api/pnjl/` magnetic 页面中哪些内容应被新主题吸收。
- [x] 固定 magnetic 主题的导出过滤集合与生成页策略。
- [x] 为后续真正落地正文页面提供可直接执行的任务顺序与 DoD。

任务目标落地说明（2026-03-10）：

- 已建立 [docs/api/models/variants/README.md](../../api/models/variants/README.md) 与 [docs/api/models/variants/magnetic/README.md](../../api/models/variants/magnetic/README.md)；
- 已建立 magnetic 主题三层与细节页：
  - [docs/api/models/variants/magnetic/Overview.md](../../api/models/variants/magnetic/Overview.md)
  - [docs/api/models/variants/magnetic/CoreConcepts.md](../../api/models/variants/magnetic/CoreConcepts.md)
  - [docs/api/models/variants/magnetic/ModelAndConfig.md](../../api/models/variants/magnetic/ModelAndConfig.md)
  - [docs/api/models/variants/magnetic/LandauIntegrals.md](../../api/models/variants/magnetic/LandauIntegrals.md)
  - [docs/api/models/variants/magnetic/Thermodynamics.md](../../api/models/variants/magnetic/Thermodynamics.md)
  - [docs/api/models/variants/magnetic/generated/Exports.md](../../api/models/variants/magnetic/generated/Exports.md)
- 总导航 [docs/api/README.md](../../api/README.md) 与脚本说明 [scripts/dev/README.md](../../../scripts/dev/README.md) 已同步更新；
- 旧页 [docs/api/pnjl/MagneticPNJL.md](../../api/pnjl/MagneticPNJL.md) 与 [docs/api/pnjl/MagneticPNJLParameters.md](../../api/pnjl/MagneticPNJLParameters.md) 已降级为迁移说明。

---

## 3. 范围与非目标

### 3.1 范围

- `variants` 总层目录方案；
- magnetic 主题的目标目录与页面结构；
- magnetic 相关公开导出的过滤集合；
- 旧 magnetic 页面与新主题之间的迁移/降级策略；
- 后续正文实施时的阶段拆分与验收标准。

### 3.2 非目标

- 不在本任务中直接创建 `docs/api/models/variants/magnetic/` 正文页面；
- 不在本任务中直接生成 magnetic 导出索引页；
- 不在本任务中修改 `src/models/` 或 `src/pnjl/` 实现代码；
- 不在本任务中扩展到其它模型变体主题；
- 不在本任务中处理 susceptibility、transport 等其它未决主题。

---

## 4. 现状盘点

### 4.1 当前公开导出簇

根据 `src/models/Models.jl`，当前 magnetic 相关公开导出至少包括：

- `PNJLMagneticModel`
- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `omega0_flavor_landau`
- `omegat_flavor_landau`
- `density_flavor_landau`
- `MagneticIMCParams`
- `default_imc_params`
- `coupling_GB`
- `MagneticConfig`
- `default_magnetic_config`
- `calculate_magnetic_omega_components`
- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`
- `calculate_magnetic_rho`
- `calculate_magnetic_number_densities`
- `magnetic_nmax_convergence_report`

### 4.2 当前旧文档资产

现有主要页面：

- [docs/api/pnjl/MagneticPNJL.md](docs/api/pnjl/MagneticPNJL.md)
- [docs/api/pnjl/MagneticPNJLParameters.md](docs/api/pnjl/MagneticPNJLParameters.md)

现有内容已覆盖：

- 数据结构与单位口径
- Landau 能级与积分 API
- 磁场热力学 API
- `PNJLMagneticModel` 的最小说明
- 配置模板与回归脚本

但仍存在问题：

- 主视角仍偏 `PNJL.*` 与历史实现路径
- 用户入口、职责核心、导出全集尚未按新架构分层
- 参数说明页与 API 页分离，但尚未围绕 `Models` 聚合导出重组

### 4.3 当前实现与测试事实源

与后续文档应对齐的真实路径包括：

- 源码：
  - `src/models/pnjl_physics/PNJLMagneticModel.jl`
  - `src/models/pnjl_physics/core/MagneticIntegrals.jl`
  - `src/models/pnjl_physics/core/MagneticThermodynamics.jl`
- 配置：
  - `config/models/pnjl/magnetic_default.toml`
- 脚本：
  - `scripts/pnjl/run_magnetic_point.jl`
  - `scripts/pnjl/run_magnetic_eb_scan.jl`
  - `scripts/pnjl/run_magnetic_fixedpoint_regression.jl`
  - `scripts/pnjl/run_magnetic_stability_scan.jl`
- 测试：
  - `tests/unit/pnjl/test_pnjl_magnetic_model.jl`
  - `tests/unit/pnjl/test_magnetic_integrals.jl`
  - `tests/unit/pnjl/test_magnetic_thermodynamics.jl`
- baseline：
  - `tests/baselines/pnjl/baseline_pnjl_magnetic_fixedpoints_smoke_v1.csv`
  - `tests/baselines/pnjl/baseline_pnjl_magnetic_fixedpoints_nightly_v1.csv`

---

## 5. 目录方案

### 5.1 variants 总层

推荐新增目录层：

- `docs/api/models/variants/`

其职责：

- 承接 `Models` 下不适合一级平铺、但又明确属于模型家族/模型变体语义的主题；
- 避免继续把“变体主题”误放成 `docs/api/models/<theme>/` 一级目录；
- 为后续 NJL / PNJL / RPNJL / magnetic 等模型族整理预留统一层级。

### 5.2 magnetic 主题落点

推荐主题目录：

- `docs/api/models/variants/magnetic/`

说明：

- 该路径是当前推荐方案，但本任务单阶段只固定方案，不在本轮直接创建大量正文文件；
- 后续正文任务应以此路径为唯一主目录，而不是回退到 `docs/api/models/magnetic/`。

---

## 6. magnetic 主题结构设计

### 6.1 三层视图映射

| 视图 | 主页面 | 说明 |
| --- | --- | --- |
| 面向用户入口 | `variants/magnetic/Overview.md` | 先说明 `PNJLMagneticModel`、配置对象、常用计算入口、零磁场退化与最短示例 |
| 职责核心 | `variants/magnetic/CoreConcepts.md` | 解释 Landau 能级、IMC 耦合、`n_max` 收敛、磁场热力学组件之间的边界 |
| 导出 API 全集 | `variants/magnetic/generated/Exports.md` | 由脚本生成的 magnetic 公开导出子集索引 |

### 6.2 建议页面结构

- `docs/api/models/variants/magnetic/README.md`
- `docs/api/models/variants/magnetic/Overview.md`
- `docs/api/models/variants/magnetic/CoreConcepts.md`
- `docs/api/models/variants/magnetic/ModelAndConfig.md`
- `docs/api/models/variants/magnetic/LandauIntegrals.md`
- `docs/api/models/variants/magnetic/Thermodynamics.md`
- `docs/api/models/variants/magnetic/generated/Exports.md`

### 6.3 页面职责建议

- `README.md`
  - 主题主入口与阅读顺序
- `Overview.md`
  - `PNJLMagneticModel`、`MagneticConfig`、`default_magnetic_config`、最短使用路径
- `CoreConcepts.md`
  - 为什么 magnetic 属于模型变体、Landau 离散化、零磁场退化、`n_max` 收敛治理
- `ModelAndConfig.md`
  - `PNJLMagneticModel`、`MagneticIMCParams`、`MagneticConfig`、参数模板
- `LandauIntegrals.md`
  - `alpha_n`、`energy_landau`、`smooth_cutoff`、`resolve_nmax_from_cutoff`、各 flavor Landau 量
- `Thermodynamics.md`
  - `coupling_GB`、`calculate_magnetic_omega*`、`calculate_magnetic_pressure`、`calculate_magnetic_rho`、`calculate_magnetic_number_densities`、`magnetic_nmax_convergence_report`

---

## 7. 旧页面吸收与降级策略

### 7.1 需要吸收的内容

从 [docs/api/pnjl/MagneticPNJL.md](docs/api/pnjl/MagneticPNJL.md) 吸收：

- 核心 API 清单
- 单位口径
- `PNJLMagneticModel` 最小说明
- `eB -> 0` 退化行为
- 相关脚本与测试索引

从 [docs/api/pnjl/MagneticPNJLParameters.md](docs/api/pnjl/MagneticPNJLParameters.md) 吸收：

- `MagneticConfig` 与 `MagneticIMCParams` 字段口径
- 推荐参数范围
- `n_max` 收敛口径
- baseline 与配置模板引用

### 7.2 降级策略

后续正文主题落地完成后：

- `docs/api/pnjl/MagneticPNJL.md` 应降级为迁移说明/跳转页；
- `docs/api/pnjl/MagneticPNJLParameters.md` 应降级为迁移说明或参数兼容说明页；
- 不再让旧页继续承担新主题的主体说明职责。

---

## 8. 自动化策略

### 8.1 导出来源

- 首选来源：`src/models/Models.jl`

### 8.2 初始过滤集合

- `PNJLMagneticModel`
- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `omega0_flavor_landau`
- `omegat_flavor_landau`
- `density_flavor_landau`
- `MagneticIMCParams`
- `default_imc_params`
- `coupling_GB`
- `MagneticConfig`
- `default_magnetic_config`
- `calculate_magnetic_omega_components`
- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`
- `calculate_magnetic_rho`
- `calculate_magnetic_number_densities`
- `magnetic_nmax_convergence_report`

### 8.3 生成目标

- 输出位置：`docs/api/models/variants/magnetic/generated/Exports.md`
- 仍使用 `scripts/dev/generate_api_export_index.jl`

---

## 9. 实施顺序

### 阶段 A：固定 variants 总层与 magnetic 落点

- [x] 在任务层确认 `docs/api/models/variants/` 为推荐总层
- [x] 在任务层确认 `docs/api/models/variants/magnetic/` 为唯一推荐落点
- [x] 避免后续再回到 `docs/api/models/magnetic/` 的旧讨论

阶段 A 落地说明（2026-03-10）：

- [docs/api/models/variants/README.md](../../api/models/variants/README.md) 已建立为总层入口；
- [docs/api/README.md](../../api/README.md) 已将 `variants` 总层与 magnetic 主题纳入主导航。

### 阶段 B：建立三层视图正文骨架

- [x] 新建 magnetic 主题目录与 `README.md`
- [x] 新建 `Overview.md`
- [x] 新建 `CoreConcepts.md`
- [x] 建立 `generated/Exports.md`

阶段 B 落地说明（2026-03-10）：

- 已建立 `README.md`、`Overview.md`、`CoreConcepts.md` 与三张细节页；
- 已通过 `scripts/dev/generate_api_export_index.jl` 生成 [docs/api/models/variants/magnetic/generated/Exports.md](../../api/models/variants/magnetic/generated/Exports.md)；
- 生成结果显示 19 个目标导出全部已在非 generated API 页面中命中。

### 阶段 C：吸收旧页内容

- [x] 将旧 `MagneticPNJL.md` 的高价值 API 说明吸收入新主题
- [x] 将旧 `MagneticPNJLParameters.md` 的参数口径与 baseline 信息吸收入新主题
- [x] 统一源码、脚本、测试与 baseline 路径引用到当前真实路径

阶段 C 落地说明（2026-03-10）：

- 单位口径、配置对象、Landau 积分接口、热力学主接口、脚本与 baseline 指引已吸收进新主题；
- 新主题引用的源码、配置、脚本、测试与 baseline 路径均已对齐当前 `src/models/`、`config/models/`、`scripts/pnjl/` 与 `tests/` 事实源。

### 阶段 D：旧页降级与导航收口

- [x] 旧 `docs/api/pnjl/` magnetic 页面降级为迁移说明
- [x] `docs/api/README.md` 纳入 `variants` 总层导航
- [x] magnetic 主题成为该能力簇唯一主入口

阶段 D 落地说明（2026-03-10）：

- [docs/api/pnjl/MagneticPNJL.md](../../api/pnjl/MagneticPNJL.md) 与 [docs/api/pnjl/MagneticPNJLParameters.md](../../api/pnjl/MagneticPNJLParameters.md) 已改写为迁移说明；
- magnetic 主题当前唯一主入口为 [docs/api/models/variants/magnetic/README.md](../../api/models/variants/magnetic/README.md)。

---

## 10. DoD

- [x] 已形成 `variants` 总层的正式目录方案
- [x] 已形成 magnetic 主题的正式目标路径与页面结构
- [x] 已明确旧 magnetic 页面应吸收与降级的范围
- [x] 已形成可直接执行的导出过滤集合
- [x] 后续正文实施可直接沿本任务单推进，而无需再次讨论目录落点

DoD 收口说明（2026-03-10）：

- `variants` 总层与 magnetic 主题目录已经实际落地，而不是停留在方案层；
- magnetic 主题的三层与细节页已经建立，旧 `pnjl` 页面也已完成降级；
- 导出全集已成功生成，并确认目标导出 19/19 均已被新主题人工文档吸收。

---

## 11. 后续直接衔接任务

本任务单完成后，下一步应直接进入：

1. `Models` 衍生量总层与 susceptibility 子主题任务单
2. `Models` transport 相关导出边界梳理说明

说明：

- magnetic 已作为 `variants` 总层的首个真实落地对象完成；
- 之后再进入其它模型变体主题时，可以直接复用这一层级结构与治理模板。