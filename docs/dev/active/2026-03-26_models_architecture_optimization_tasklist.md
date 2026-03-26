# Models 架构优化设计任务单（最小风险版）

## 1. 背景与目标

- 背景：`src/models` 已形成“统一入口 + 多模型实现 + workflow/solver/phase”的主干架构，但仍存在少量高耦合点（`Main.*` 直接依赖、重复缓存、上层直连 `PNJLCore`）。
- 目标：在**不移动核心模型目录**（`njl/`、`pnjl_physics/`、`rpnjl/`）的前提下，降低耦合、提升可维护性，形成可分阶段落地的重构路径。
- 目标结果：
  - 依赖方向更清晰（上层通过 `Models` 稳定接口，而非 `Main.*`/内部核心模块散点访问）；
  - 求解器与条件模块的 model cache 逻辑统一；
  - 文档结构与代码结构一致（`variants` 主题补齐）。

## 2. 范围 / 非范围

### 2.1 本期范围（MVP）

- `src/models` 内部架构治理与文档设计。
- 输出分阶段任务清单、测试门禁与验收标准。
- 优先顺序：
  1. 收敛 `Main.*` 依赖；
  2. 统一 solver 层 model cache；
  3. 下沉 `PNJLCore` 直接依赖；
  4. 聚合 `Models.jl` 导出层；
  5. 工厂注册化；
  6. `variants` 文档补齐。

### 2.2 非范围（本期不做）

- 不将 `src/models/njl`、`src/models/pnjl_physics`、`src/models/rpnjl` 直接迁移到 `src/models/variants/`。
- 不改动物理公式/参数含义，不调整数值基线。
- 不进行跨域大重构（如 `src/relaxtime` 全局目录改造）。

## 3. 现状盘点（已实现 / 待实现）

### 3.1 已实现（现有基础）

- [x] 已建立统一入口：`src/models/Models.jl`、`src/models/entrypoints.jl`。
- [x] 已建立变体目录并落地 `rotation` / `gas_liquid`：`src/models/variants/`。
- [x] 工厂已统一模型构造入口：`src/models/factory.jl`。

### 3.2 待优化（结构性缺口）

- [ ] workflow 与部分模型仍存在 `Main.*` 直接访问（如 `src/models/workflows/TransportWorkflow.jl`、`src/models/workflows/MesonMassWorkflow.jl`、`src/models/pnjl_physics/PNJLMagneticModel.jl`）。
- [ ] `src/models/solver/Conditions.jl` 与 `src/models/solver/ImplicitSolver.jl` 存在重复 model cache 与 model 获取策略。
- [ ] 上层模块仍直接依赖 `PNJLCore` 默认网格常量（耦合实现细节）。
- [ ] `Models.jl` include/export 体量较大，聚合层次可进一步清晰化。
- [ ] `create_model` 仍是 `if-elseif` 分支，可扩展性一般。
- [ ] API 文档在 `variants` 层仅明确 magnetic，rotation/gas_liquid 主题不足。

## 4. 任务分解（可勾选）

## 阶段 A：依赖方向收敛（高优先级）

- [ ] A1. 新增/补齐 `Models` 内部稳定 accessor（例如 workflow 需要的节点、常量、模块访问器），避免 workflow 直接依赖 `Main.*`。
  - 验收：`TransportWorkflow` 与 `MesonMassWorkflow` 中关键 `Main.Models.PNJLCore/PNJLIntegrals` 直接访问显著减少。
- [ ] A2. 梳理 `PNJLMagneticModel` 对 `Main.MagneticThermodynamics` 的依赖边界，优先封装为 `Models` 内部门面函数。
  - 验收：模型侧调用点集中在有限门面函数，外部行为不变。

## 阶段 B：solver cache 统一（高优先级）

- [ ] B1. 抽取统一 model cache/getter（可放在 `src/models/factory.jl` 或新增轻量 `model_registry.jl`）。
  - 验收：`Conditions` 与 `ImplicitSolver` 不再各自维护重复 `_MODEL_CACHE`。
- [ ] B2. 统一并发访问策略（必要时加入锁或明确线程模型说明）。
  - 验收：接口行为一致，回归测试通过。

## 阶段 C：核心实现细节下沉（中优先级）

- [ ] C1. 为积分网格与默认节点定义稳定契约（经 `entrypoints` 或专用 accessor 暴露）。
  - 验收：workflow 不再直接取 `PNJLCore.DEFAULT_*`。
- [ ] C2. 将 `PNJLCore` 作为“内部实现依赖”逐步限制在模型/求解核心，减少向上游扩散。
  - 验收：新增调用方默认使用契约层接口。

## 阶段 D：聚合与工厂优化（中优先级）

- [ ] D1. 拆分 `Models.jl` 导出聚合段（如 core/solver/workflows/variants 分块），保持 include 顺序安全。
  - 验收：文件可读性提升，不引入 world-age/加载顺序回归。
- [ ] D2. 将 `create_model` 升级为注册表机制（保留旧符号兼容）。
  - 验收：新增模型不再需要修改长 `if-elseif` 主分支。

## 阶段 E：文档对齐（中优先级）

- [ ] E1. 在 `docs/api/models/variants/` 新增 `rotation` 与 `gas_liquid` 主题页（至少 README + Overview + generated exports）。
  - 验收：代码与 API 导航语义一致。
- [ ] E2. 更新 `docs/api/models/variants/README.md` 的主题列表。
  - 验收：`variants` 总览可直接发现全部主变体。

## 5. 测试与验收标准

### 5.1 最小回归集（每阶段至少执行）

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

### 5.2 治理脚本（文档/迁移门禁）

- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] `julia --project=. scripts/dev/check_pnjl_migration_guard.jl`

### 5.3 分阶段验收口径

- A 阶段：依赖收敛后，核心 workflow 行为与结果无语义回归。
- B 阶段：cache 行为统一且并发语义明确，无重复实例化异常。
- C 阶段：上层不再依赖 `PNJLCore` 常量细节。
- D 阶段：聚合与工厂改造不改变既有 public entrypoints。
- E 阶段：API 文档可完整反映 variants 主题结构。

## 6. 里程碑

- [ ] M1（短期，1-2 次 PR）：完成阶段 A + B。
- [ ] M2（中期，1-2 次 PR）：完成阶段 C + D。
- [ ] M3（收尾，1 次 PR）：完成阶段 E，并进行文档一致性检查。

## 7. DoD（Definition of Done）

- [ ] 所有本期范围任务勾选完成，且非范围未被误改。
- [ ] `Models` 入口兼容性保持（`create_model`、`entrypoints`、核心导出不破坏）。
- [ ] smoke 单元/集成/回归测试通过。
- [ ] 文档门禁脚本通过，`docs/api/models/variants` 导航可用。
- [ ] 迁移与兼容说明完整，便于后续归档。

## 8. 风险与回退方案

### 8.1 主要风险

- 依赖收敛可能触发 include 顺序与 world-age 问题。
- cache 统一若处理不当，可能引入并发可见性问题。
- 上层去耦可能短期增加适配层复杂度。

### 8.2 回退策略

- 每阶段独立 PR，禁止跨阶段大包提交。
- 为关键入口保留兼容别名（过渡期），并附 deprecation 注记。
- 若阶段内出现回归，按阶段粒度回退，不影响未变更模块。

## 9. 本文档与实施关系

- 本文档是实施前设计任务单，不代表代码已改动。
- 后续实施建议按阶段推进，并在每个阶段结束后更新本文件的勾选状态与验证记录。
