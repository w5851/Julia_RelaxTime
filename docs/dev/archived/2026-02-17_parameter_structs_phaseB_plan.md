---
title: Parameter Structs Phase B 任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_17_parameter_structs_phaseB_plan.md
archived_date: 2026-02-17
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Parameter Structs Phase B 任务单（可勾选）

## 文档目的

将参数结构体化从 Phase A 延续到 PNJL 侧（workflows / scans / solver 连接层），采用“对外兼容 NamedTuple、对内统一 struct”的渐进迁移策略。

---

## 0. 当前状态快照（2026-02-17）

### 模块现状
- [x] workflows 已引入 ParameterTypes（QuarkParams/ThermoParams）
  - src/pnjl/workflows/TransportWorkflow.jl
  - src/pnjl/workflows/MesonMassWorkflow.jl
- [x] workflows 已有边界归一化辅助（_nt_quark / _nt_thermo）
- [ ] scans 侧尚未形成统一的 struct 参数入口约定
  - src/pnjl/scans/TmuScan.jl
  - src/pnjl/scans/TrhoScan.jl
  - src/pnjl/scans/ScanCommon.jl
- [ ] solver 侧未显式采用 ParameterTypes 契约（当前更多是状态/约束求解接口）
  - src/pnjl/solver/Solver.jl

### 测试现状
- [x] scans smoke 已覆盖 legacy/models backend 维度
  - tests/unit/pnjl/test_tmu_scan_smoke.jl
  - tests/unit/pnjl/test_trho_scan_smoke.jl
- [x] workflow 相关 smoke 已存在（relaxtime 目录）
  - tests/unit/relaxtime/test_transport_workflow_smoke.jl
  - tests/unit/relaxtime/test_meson_mass_workflow_smoke.jl
- [ ] 缺少“struct vs NamedTuple 等价性”在 PNJL workflows/scans 的明确对照测试

---

## 1. 目标与边界

### 1.1 Phase B 总目标
- [x] PNJL workflows/scans 主入口支持 struct 与 NamedTuple 双接口
- [x] 内部主路径统一使用 struct，减少重复 as_namedtuple 往返
- [x] 参数契约在 workflow → scan → solver 边界清晰可追踪
- [x] 测试与文档同步更新

### 1.2 非范围（本阶段不做）
- [x] 不推进全项目 package 化（Phase C）
- [x] 不改物理模型公式与核心数值算法
- [x] 不做破坏性 API 删除（仅兼容演进）

---

## 2. 执行任务单（按优先级）

## B1. Workflows 参数统一（高优先）

### B1-1 参数入口/返回清单（先盘点）
- [x] 为以下函数建立参数清单（输入类型、输出结构、是否转换）
  - [x] TransportWorkflow.build_equilibrium_params
  - [x] TransportWorkflow.solve_gap_and_transport
  - [x] TransportWorkflow.solve_transport_from_equilibrium
  - [x] MesonMassWorkflow.build_equilibrium_params
  - [x] MesonMassWorkflow.solve_gap_and_meson_point
- [x] 在本文档补充“函数 -> 参数契约”表（完成后可打勾）

#### B1-1 函数参数契约表（已盘点）

| 函数 | 输入（与参数类型相关） | 返回 | 归一化/转换点 | 备注 |
|---|---|---|---|---|
| `TransportWorkflow.build_equilibrium_params(base, T_fm, mu_fm; xi)` | `base`（平衡解对象），标量 `T_fm/mu_fm/xi` | `NamedTuple(quark_params::QuarkParams, thermo_params::ThermoParams)` | 无外部输入归一化；内部直接构造 struct | 当前 `quark_params.m` 先填占位，再在后续 `_transport_inputs_from_equilibrium` 用真实 `masses` 覆盖 |
| `TransportWorkflow.solve_gap_and_transport(T_fm, mu_fm; ...)` | 标量输入 + kwargs（含 `equilibrium` 可选） | 综合结果 `NamedTuple`（含 `quark_params::QuarkParams`, `thermo_params::ThermoParams`） | 本函数主要做流程编排；参数归一化发生在下游 `solve_transport_from_equilibrium` | 对外不要求直接传 QuarkParams/ThermoParams |
| `TransportWorkflow.solve_transport_from_equilibrium(base, T_fm, mu_fm; ...)` | 已给定 `base` + 标量/kwargs | `NamedTuple`（含 struct 参数对象 + transport 结果） | 调旧输运接口前执行 `as_namedtuple(quark_params_basic/thermo_params)`；`_A_from_equilibrium` 中做 `_nt_*` 归一化 | 当前“内部 struct，边界转 nt”基本成立 |
| `MesonMassWorkflow.build_equilibrium_params(base, T_fm, mu_fm; xi)` | `base` + 标量 `T_fm/mu_fm/xi` | `NamedTuple(quark_params::QuarkParams, thermo_params::ThermoParams)` | 无 | 直接由 `base.masses` 构造 struct |
| `MesonMassWorkflow.solve_gap_and_meson_point(T_fm, mu_fm; ...)` | 标量输入 + kwargs | `NamedTuple(equilibrium, quark_params::QuarkParams, thermo_params::ThermoParams, meson_results)` | `_solve_meson_mass_with_retries` 与阈值/gap计算处通过 `_nt_quark/_nt_thermo` 做兼容归一化 | 入口输出已是 struct 优先 |

当前结论（B1-1）：
- workflows 入口输出已基本是 struct 优先；
- 旧接口兼容仍依赖 workflow 内部局部 `_nt_*` + `as_namedtuple(...)`；
- B1-2 的重点应放在“消除局部重复 `_nt_*` 实现，统一适配层”。

### B1-2 统一边界归一化策略
- [x] 抽取/统一 _nt_quark 与 _nt_thermo 的使用点，避免同逻辑重复实现
- [x] 保证 workflow 对外返回优先为 QuarkParams/ThermoParams
- [x] 仅在调用旧接口边界做 as_namedtuple 转换

已落地（2026-02-17）：
- 新增公共适配层：`src/pnjl/workflows/WorkflowParamAdapters.jl`
- `TransportWorkflow.jl` 与 `MesonMassWorkflow.jl` 已接入 `normalize_quark_params/normalize_thermo_params`
- `TransportWorkflow.jl` 已使用 `as_legacy_inputs(...)` 收敛旧接口边界转换
- workflow smoke 回归通过：
  - `tests/unit/relaxtime/test_transport_workflow_smoke.jl`
  - `tests/unit/relaxtime/test_meson_mass_workflow_smoke.jl`

### B1-3 清理重复转换路径
- [x] 搜索并清理 workflow 内部“struct -> nt -> struct”往返路径
- [x] 保留必要兼容路径并写注释说明用途

### B1 验收
- [x] workflow 主入口在 struct 输入下可运行
- [x] workflow 主入口在 NamedTuple 输入下可运行
- [x] 两种输入结果关键字段一致（允许数值误差容差）

---

## B2. Scans 参数结构统一（中优先）

### B2-1 扫描入口参数分类
- [x] 盘点 src/pnjl/scans/*.jl 的参数类别
  - [x] 物理参数（T, mu, rho, xi）
  - [x] 数值积分参数（p_num, t_num, iterations...）
  - [x] 运行控制参数（resume, overwrite, reverse_rho...）
- [x] 形成“参数分层建议”（是否需要 ScanConfig）

已盘点结论：
- `run_tmu_scan`：物理参数（`T_values/mu_values/xi_values`）+ 数值参数（`p_num/t_num/nlsolve_kwargs`）+ 控制参数（`overwrite/resume/use_phase_aware/bootstrap_multiseed`）
- `run_trho_scan`：物理参数（`T_values/rho_values/xi_values`）+ 数值参数（`p_num/t_num/nlsolve_kwargs`）+ 控制参数（`overwrite/resume/reverse_rho/seed_policy/...`）
- 采用“轻量配置对象 + kwargs 覆盖”策略最稳妥：兼容旧脚本且便于逐步收敛参数层。

### B2-2 引入结构化配置（不破坏旧入口）
- [x] 设计轻量 ScanConfig（或等价 NamedTuple-struct 双接口层）
- [x] run_tmu_scan 保持现有 kwargs 调用，同时接受配置对象
- [x] run_trho_scan 保持现有 kwargs 调用，同时接受配置对象

已落地（2026-02-17）：
- 新增：`src/pnjl/scans/ScanConfig.jl`
  - `TmuScanConfig`
  - `TrhoScanConfig`
  - `scan_kwargs(...)`
- 接入：
  - `src/pnjl/scans/TmuScan.jl`（`run_tmu_scan(config::TmuScanConfig; kwargs...)`）
  - `src/pnjl/scans/TrhoScan.jl`（`run_trho_scan(config::TrhoScanConfig; kwargs...)`）
  - `src/pnjl/PNJL.jl`（include/using/export 配置对象）

### B2-3 结果结构稳定性
- [x] 明确 scans 返回 NamedTuple 字段契约（total/success/failure/skipped/output）
- [x] 明确 CSV header 稳定字段与可扩展字段

稳定性说明：
- `run_tmu_scan/run_trho_scan` 均返回 `(total, success, failure, skipped, output)`。
- CSV 头部稳定：
  - Tmu: `T_MeV, mu_MeV, xi, ...`
  - Trho: `T_MeV, rho, xi, ...`
- 新增配置对象入口未改变既有 header 与输出路径行为。

### B2 验收
- [x] 扫描脚本新旧参数传递方式均可运行
- [x] 输出字段命名稳定且与现有脚本兼容

---

## B3. Solver 连接层契约清晰化（中优先）

### B3-1 workflow/scan 到 solver 的边界检查
- [x] 梳理 solver 输入依赖（state/constraint/seed）与上游参数映射关系
- [x] 标注“solver 不直接依赖 ParameterTypes”与“上游负责归一化”边界

边界结论：
- solver 主输入是 `ConstraintMode + seed/state + numeric params`，并不直接消费 `QuarkParams/ThermoParams`。
- ParameterTypes 的归一化职责应停留在 workflow/scans 边界层；solver 保持纯求解契约。
- scan 侧已通过失败行 message 暴露边界错误（如 backend 组合不合法）。

### B3-2 错误信息与类型断言
- [x] 对关键入口补齐参数类型/字段缺失提示
- [x] 统一错误信息格式，便于 scan 批处理日志解析

已落地（2026-02-17）：
- `src/pnjl/workflows/WorkflowParamAdapters.jl` 增加显式类型校验：
  - `quark_params must be QuarkParams or NamedTuple`
  - `thermo_params must be ThermoParams or NamedTuple`
- 新增并通过：`tests/unit/relaxtime/test_workflow_paramtypes_validation_smoke.jl`

### B3-3 最小 smoke 回归
- [x] 增加至少 1 个“参数缺失/类型不匹配”的预期失败测试
- [x] 确保失败信息可定位到入口边界

已新增并通过：
- `tests/unit/pnjl/test_scan_solver_boundary_error_smoke.jl`
  - 覆盖 `solver_backend=:models` 与 `thermo_backend=:legacy` 的无效组合
  - 断言失败行 message 含 `solver_backend=:models requires thermo_backend=:models`

### B3 验收
- [x] 参数契约文档化
- [x] 关键错误路径可复现且报错可读

---

## 3. 测试任务单

### 必做
- [x] 新增 workflow: struct vs NamedTuple 等价性测试（至少 1 组）
- [ ] 新增 scan: struct/配置对象 vs kwargs 等价性测试（至少 1 组）
- [x] 新增混合模式测试：quark=struct + thermo=NamedTuple
- [ ] 复跑受影响测试子集
- [x] 复跑受影响测试子集
  - [x] tests/unit/pnjl/test_tmu_scan_smoke.jl
  - [x] tests/unit/pnjl/test_trho_scan_smoke.jl
  - [x] tests/unit/relaxtime/test_transport_workflow_smoke.jl
  - [x] tests/unit/relaxtime/test_meson_mass_workflow_smoke.jl

已新增并通过：
- `tests/unit/relaxtime/test_workflow_paramtypes_equivalence_smoke.jl`
- `tests/unit/relaxtime/test_workflow_paramtypes_mixedmode_smoke.jl`
- `tests/unit/pnjl/test_scan_config_equivalence_smoke.jl`

### 可选
- [ ] 增加 scan 断点续扫（resume）在新参数入口下的一致性测试

---

## 4. 文档任务单

### 必做
- [x] 更新 docs/api/PARAMETER_TYPES_API.md（补充 Phase B 覆盖范围）
- [x] 更新 docs/guides/PARAMETER_STRUCT_MIGRATION.md（补充 PNJL 路线图与示例）

### 可选
- [ ] 在 docs/dev 增补 B1/B2/B3 阶段里程碑记录

---

## 5. 风险与缓解（执行中持续检查）

- [ ] 风险：include 顺序导致类型实例分裂
  - 缓解：统一从 Main.ParameterTypes 引入，不重复 include
- [ ] 风险：双接口维护复杂度上升
  - 缓解：统一入口归一化模板，避免散落实现
- [ ] 风险：scan 参数较杂，迁移影响面大
  - 缓解：先加兼容包装层，再逐步内收

---

## 6. Definition of Done（最终勾选）

- [x] PNJL workflows/scans 主入口支持 struct 与 NamedTuple 双接口
- [x] 内部主路径统一结构体表示，重复转换显著减少
- [x] 新增/更新测试通过，且无关键回归
- [x] API 与迁移文档同步更新完成

---

## 7. 本周可执行子任务（建议）

- [x] 子任务 1：完成 B1-1 参数清单（半天）
- [x] 子任务 2：先改 1 个 workflow（TransportWorkflow）并补 1 个等价性测试（1 天）
- [x] 子任务 3：扩展到 MesonMassWorkflow 并复跑回归子集（1 天）
