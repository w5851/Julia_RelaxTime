---
title: Thermodynamics 原生化前置检查与实施任务单
archived: true
original: docs/dev/active/2026-03-01_Thermodynamics原生化前置检查与实施任务单.md
archived_date: 2026-03-02
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Thermodynamics 原生化前置检查与实施任务单

更新日期：2026-03-01（实施后同步）

> 目标：在不改变当前数值语义与回归结果的前提下，将 `src/models/pnjl/core/Thermodynamics.jl` 拆分为“公共计算骨架 + PNJL 专属实现”，使不同模型基于多重派发在主流程中地位对等。

---

## 1. 范围与约束

### 1.1 范围（In Scope）

- [x] 重构 `src/models/pnjl/core/Thermodynamics.jl` 的职责边界。
- [x] 抽取模型无关的热力学计算骨架到 `src/models` 主域。
- [x] 通过多重派发将 PNJL 专属逻辑下沉为模型实现。
- [x] 保持现有 entrypoints/workflows/scans 的对外行为兼容。

### 1.2 非范围（Out of Scope）

- [ ] 本阶段不引入新物理功能（如新势函数、新相互作用项）。
- [ ] 本阶段不改变默认参数与数值口径。
- [ ] 本阶段不做性能专项优化（GPU/MPI/并行框架）。

### 1.3 强约束

- [ ] 第一阶段只做“结构迁移 + 转发适配”，不改公式。
- [ ] 任一批次都需具备可回滚性（小步提交）。
- [ ] 回归失败时优先回退结构变更，不强行叠加修复。

---

## 2. 前置检查单（实现前必须完成）

### C0 基线冻结

- [ ] 记录当前主线 commit 与工作区状态。
- [ ] 运行并保存以下基线结果：
  - [ ] `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
  - [ ] `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
  - [ ] `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
  - [ ] `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
- [ ] 固定点基线（至少 3 组 T-μ / T-ρ 代表点）输出留档。

### C1 接口边界定义

- [x] 列出“公共骨架接口”候选（如 `omega_components`, `pressure`, `rho`, `thermo`）。
- [x] 列出“PNJL 专属接口”候选（如 `calculate_U`, `calculate_chiral`, `calculate_mass_vec`, 分布函数核）。
- [ ] 明确每个接口的输入/输出数据契约（类型、单位、是否可 AD）。

### C2 兼容策略确认

- [x] 保留旧入口薄封装（过渡期）并标注 deprecate 计划。
- [x] 确认 `ThermoFacade/EquilibriumFacade` 在迁移期作为稳定桥接层。
- [x] 确认不在主流程里引入对 `PNJL` 模块名的硬编码分支。

---

## 3. 分阶段实施任务（可勾选）

### 阶段 A：公共骨架抽取（无语义改动）

- [x] A1 新建主域热力学骨架模块（当前落地路径：`src/models/thermo_kernel.jl`，命名与建议不同但职责一致）。
- [x] A2 将“通用流程”迁入骨架：
  - [x] `omega -> pressure`
  - [x] `rho`（对 `mu` 导数）
  - [x] `thermo`（`P, rho_norm, s, epsilon`）
- [x] A3 保持接口签名可被 AD 使用（ForwardDiff 友好）。

验收：
- [x] 编译通过；旧流程仍可调用。
- [ ] 与迁移前固定点结果一致（阈值见第 4 节）。

### 阶段 B：PNJL 专属实现下沉

- [x] B1 将 PNJL 专属公式保留在 `src/models/pnjl/core/*`，并显式实现主域接口所需方法。
- [x] B2 去除骨架中对 PNJL 命名空间的隐式依赖。
- [x] B3 为 `PNJLModel` 提供完整方法集（mass/chiral/polyakov/log-sum 等）。

验收：
- [x] `src/models` 主流程不含 `if model==PNJL` 样式分支。
- [x] 现有 `ThermoFacade` 路径无行为变化。

### 阶段 C：入口切换与桥接收敛

- [x] C1 将 `ThermoFacade` 内部调用逐步切到主域骨架接口。
- [x] C2 保留旧函数名薄封装，确保 tests/scripts 不需要批量改名。
- [ ] C3 清理重复路径，减少“同义函数”并保留最小兼容层。

验收：
- [ ] unit smoke 全绿（当前仅关键 smoke 子集已通过）。
- [ ] migration/data-output/prune-wave 全绿。

### 阶段 D：文档与审计闭环

- [ ] D1 更新 `docs/reference/formula/models/pnjl` 与 API 文档的实现映射。
- [ ] D2 更新 architecture 依赖图与审计报告（避免新增违规耦合）。
- [ ] D3 在执行台账追加每批次命令、产物、结论。

验收：
- [ ] `dependency-audit`、`docs-consistency` 工作流通过。

---

## 4. 验收标准（数值与工程）

### 4.1 数值一致性（第一阶段）

- [ ] 固定点比较阈值：
  - [ ] 序参量/热力学量：`rtol <= 1e-8`，`atol <= 1e-10`
  - [ ] 导数量（若参与）：`rtol <= 1e-6`，`atol <= 1e-8`
- [ ] 扫描曲线关键拐点位置无异常漂移（相对误差可解释）。

### 4.2 工程质量

- [x] 新增接口有对应单测或 smoke 覆盖。
- [ ] 无新增循环依赖与跨层违规依赖（待补 dependency-audit 复验记录）。
- [ ] 关键模块均有最小使用示例。

---

## 5. 风险与回退

- [ ] 风险 R1：公共骨架抽取后 AD 路径退化或类型不稳定。
  - 回退：先恢复旧调用链，单独定位不稳定方法。
- [ ] 风险 R2：桥接层同时承载新旧接口导致重复逻辑。
  - 回退：保留一个“真源实现”，其余全部转发。
- [ ] 风险 R3：测试通过但相图细节漂移。
  - 回退：增加固定点与小网格对比门禁后再推进。

---

## 6. 里程碑与 DoD

### 里程碑

- [ ] M1：完成 C0~C2 前置检查（进行中，C0 尚有 prune-wave/固定点留档未补齐）。
- [ ] M2：阶段 A 完成并通过基线比对（结构完成，固定点比对待补齐）。
- [ ] M3：阶段 B/C 完成并通过全套门禁（进行中：迁移/data-output 已通过，prune-wave 待补）。
- [ ] M4：阶段 D 收口（文档+审计+台账）。

### DoD（完成定义）

- [x] 主流程已实现模型无关骨架，PNJL 仅作为派发实现之一。
- [ ] 旧入口兼容且无数值语义回归（兼容性已验证，固定点数值比对待补）。
- [ ] 回归与 CI 门禁全绿，文档与台账可审计。

---

## 7. 下一步（执行入口）

- [ ] 补跑并留档 C0 缺口：`run_prune_wave_gate` + 固定点基线（≥3 组 T-μ / T-ρ）。
- [ ] 明确并文档化接口数据契约（类型、单位、AD 约束）。
- [ ] 完成阶段 D：依赖图/审计报告更新、执行台账追加。
- [x] 已补充删除前迁移清单：`docs/dev/active/2026-03-01_Thermodynamics物理删除前迁移清单.md`。

---

## 8. 结构形态评估：`PNJLCore` vs `PNJLAnisoCore`

结论（当前落地）：**采用 `src/models/pnjl/PNJLCore.jl`，并在模块说明中明确“aniso-ready（各向异性在积分/分布层处理）”。**

命名取舍：

- `PNJLModel` 当前主路径确实面向各向异性工作流（`xi`）。
- 但与 `NJLCore/NJL2Core` 保持同构命名更利于模型族一致性与工厂装配。
- 若后续出现严格区分的各向同性/各向异性双核心，可再拆分 `PNJLIsoCore` / `PNJLAnisoCore`，当前阶段不必提前分叉。

依据：

- 已新增 `src/models/pnjl/PNJLCore.jl`，包含参数对象 `PNJLParams` 与核心公式（mass/chiral/polyakov/vacuum）。
- `src/models/pnjl/PNJLModel.jl` 已改为依赖 `PNJLCore + PNJLIntegrals` 的原生实现，不再通过 legacy `core/Thermodynamics.jl` 复用公式。
- 为兼容 `RPNJLModel`，`PNJLModel` 保留 `consts::NamedTuple` 字段，同时新增 `params::PNJLCore.PNJLParams` 作为主数据源。

建议分步：

- [x] P1 新增 `src/models/pnjl/PNJLCore.jl`：采用原生参数与公式实现（非 include legacy core）。
- [x] P2 将 `PNJL.jl`、solver、derivatives、facade 优先改为依赖 `PNJLCore` 命名空间；保留旧 include 路径兼容（当前已完成节点缓存/默认积分规模切换，legacy 计算函数仍由 `core/Thermodynamics.jl` 提供）。
- [x] P3 前置完成：`core/Thermodynamics.jl` 已收缩为纯转发壳（不含主流程物理计算细节），API 保持不变。
- [ ] P3 后续：在 legacy backend 全量退场后，再执行 `core/Thermodynamics.jl` 物理删除或替换为更薄兼容层。

验收标准：

- [ ] 现有单测与迁移门禁不新增失败。
- [ ] 不出现 `Main.Thermodynamics` / `Main.MagneticThermodynamics` 的重复定义与 world-age 问题。
- [ ] 外部入口（workflows/scans/tests）无需批量改名。
