# ADR-0004: Solver 三层契约与向量内核统一

## 状态

提议

## 背景

当前 solver 相关能力在仓库内已形成较清晰分层基础，但仍存在口径分散与边界混用问题：

- mode 语义与求解内核仍有耦合，`FixedMu` 与非 `FixedMu` 治理面不完全一致；
- 导数能力虽已有独立实现基础（`implicit_gap.jl`），但对外入口仍与 solver 层有历史包装耦合；
- 设计评审中明确提出：对外接口可采用命名语义输入，但数值与微分内核应统一向量契约。

同时，已有 5 份 active 文档已达成一致方向：

- 目标态蓝图；
- ConstraintSpec 草案；
- 导数解耦草案；
- 三层最小接口清单；
- 项目级时序图。

这些内容属于架构级决策，需提升为 ADR，避免执行阶段反复漂移。

## 决策

我们决定采用以下统一架构约束：

1. 采用三层模型：
   - 语义外壳层（NamedTuple 面向业务）；
   - Schema 桥接层（NamedTuple <-> Vector）；
   - 向量内核层（solver/AD/隐式微分共享）。

2. 统一内核契约为：
   - `residual_vec!(F, x_vec, theta_vec, cfg)`；
   - solver 与导数引擎复用同一向量残差契约。

3. mode 降格为约束配置输入：
   - 通过 `ConstraintSpec` 统一表达约束问题；
   - mode 作为兼容构造器和标记，不再驱动主架构分支。

4. 统一策略入口：
   - 采用 `primary_strategy` 作为 method + multi-seed + fallback 的单一组合入口；
   - 不再把 multi-seed 作为并列散点能力。

5. fallback 插件化边界：
   - legacy/weighted fallback 作为显式可开关插件能力；
   - 默认不污染主链求解语义。

6. 维度契约最小化：
   - solver 内核仅维护 `x_dim/theta_dim`；
   - 物理语义（如 state/mu 切片）通过回调与 schema 表达，不进入内核维度字段。

## 理由

- 与数值库现实对齐：NLsolve、ForwardDiff、隐式微分链路均更稳定于向量内核。
- 降低维护复杂度：模型差异进入 schema 与 residual 语义层，不复制转换算法。
- 提升扩展性：新增模型/约束优先新增 schema/spec，不改内核主干。
- 提升可验证性：统一契约后可建立一致的 parity/regression 验收矩阵。

## 备选方案

### 方案 A：全链路 NamedTuple 内核化

优点：业务语义直观。

缺点：与现有数值/微分库主路径不匹配，性能与类型稳定风险更高。

未选择原因：不符合当前工程稳定性与 AD 可靠性要求。

### 方案 B：维持当前 mixed 模式（不冻结三层契约）

优点：短期改动最少。

缺点：口径持续漂移，执行阶段反复返工风险高。

未选择原因：与当前跨文档统一目标冲突。

## 后果

### 正面影响

- 求解与导数共享统一内核契约，边界清晰。
- mode 迁移为 spec-first 后，主链治理面可收敛一致。
- fallback 行为可审计、可开关、可逐步退役。

### 负面影响

- 需要引入 schema 注册与适配基础设施。
- 迁移期存在新旧接口并存维护成本。

### 风险与缓解

- 风险：边界转换实现不当导致性能回退。
  - 缓解：采用“零重复转换”准则，转换仅发生在边界。
- 风险：AD 路径被类型强转破坏。
  - 缓解：禁止导数路径 `Float64(...)` 强转，保留泛型元素类型。
- 风险：spec-first 与旧 mode 路径结果不一致。
  - 缓解：A/B parity 开关 + integration/regression 门禁。

## 实施指引

- 执行主路线参照：
  - `docs/dev/active/2026-04-06_solver_5docs_integrated_execution_roadmap.md`
- 实现基线参照：
  - `docs/dev/active/2026-04-06_solver_three_layer_minimal_interfaces.md`

## 相关决策

- [ADR-0002](0002-models-solver-contract.md)
- [ADR-0003](0003-pnjl-solver-decoupling-governance-and-ad-implicit-contract.md)

## 参考资料

- `docs/dev/active/2026-04-06_solver_target_architecture_blueprint.md`
- `docs/dev/active/2026-04-06_solver_constraintspec_interface_draft.md`
- `docs/dev/active/2026-04-06_solver_derivative_decoupling_interface_and_migration.md`
- `docs/dev/active/2026-04-06_solver_three_layer_minimal_interfaces.md`
- `docs/dev/active/2026-04-06_solver_project_sequence_diagram.md`

---

**日期**：2026-04-06
