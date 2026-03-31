# ADR-0003: PNJL 求解器解耦治理与 AD/隐式微分契约

## 状态

提议

## 背景

当前仓库的 PNJL 求解能力已经形成 `Models` 统一入口，但在约束求解与导数链路上仍存在多条并行实现路径：

- `src/models/constraint_solver.jl` 与 `src/models/solver/ImplicitSolver.jl` 存在职责重叠；
- `build_residual!` 闭包风格已可用，但“显式参数接口”契约尚未冻结；
- AD/隐式微分链路可工作，但配置与边界条件存在全局状态和可重入风险。

同时，`docs/dev/active/2026-03-31_PNJL求解器解耦框架约束与分层草案.md` 与后续评审计划已明确：该事项属于跨模块架构治理，不是局部实现细节。

## 决策

我们决定：

1. 将 PNJL 求解器解耦纳入正式 ADR 管理，作为架构治理决策；
2. 在 Wave-A 先冻结契约，再推进实现收敛；
3. 采用三层文档分工，避免多文档漂移：
   - `docs/decisions/`：记录“为什么这样定”（稳定决策）；
   - `docs/superpowers/specs|plans/`：记录“怎么落地”（执行设计与计划）；
   - `docs/dev/active/`：记录项目内任务勾选、验证证据与进度状态。

4. Wave-A 必须先固化以下契约：
   - ProblemSpec（mode/x_dim/theta_dim/residual!/forward/conditions/postprocess/governance）；
   - AD/ID 契约（Dual-safe 条件路径、forward/conditions 统一语义）；
   - Candidate Governance 契约（hard_rule/selector/context）；
   - 兼容层迁移台账（旧入口 -> 新入口 -> 删除门槛）。

## 理由

- 风险优先：当前主要风险是“口径漂移与治理失配”，不是底层数值能力不足。
- 兼容优先：先冻结契约可在不重写后端的情况下推进并轨。
- 证据优先：用最小验证矩阵（Jacobian smoke、维度一致性、选择稳定性、配置隔离）约束演进。

## 备选方案

### 方案 A：仅保留 active 任务单，不建 ADR

优点：流程简单。

缺点：缺少“长期稳定的决策记录”，后续跨会话讨论容易回到同一争议。

未选择原因：与本事项的架构级影响不匹配。

### 方案 B：直接代码并轨，不先冻结契约

优点：短期改动速度快。

缺点：实现先行会放大接口漂移和回归风险。

未选择原因：与当前“先固化约束再收敛实现”的治理目标冲突。

## 后果

### 正面影响

- 决策依据、执行计划与任务证据职责边界清晰。
- 后续新增 mode/参数扩展时有统一契约可遵循。
- AD/ID 风险点可被系统化测试门槛覆盖。

### 负面影响

- 文档管理层级增加，需要保持链接与状态同步。

### 风险与缓解

- 风险：三层文档信息不一致。
  - 缓解：active 主任务单作为执行单一事实源；ADR 只记决策，不记勾选进度；spec/plan 只记实现方案。

## 实施计划

1. 建立 `docs/superpowers/specs` 设计文档与 `docs/superpowers/plans` 实施计划。
2. 在 active 主任务单中加入 ADR/spec/plan 链接，作为执行入口索引。
3. Wave-A 实现按计划推进，并在 active 文档回填验证证据。

## 相关文档

- `docs/decisions/0002-models-solver-contract.md`
- `docs/dev/active/2026-03-31_PNJL求解器解耦框架约束与分层草案.md`
- `docs/dev/active/2026-03-31_PNJL求解器解耦框架兼容性评审与Wave-A实现计划.md`
- `docs/superpowers/specs/2026-03-31-pnjl-solver-decoupling-wave-a-design.md`
- `docs/superpowers/plans/2026-03-31-pnjl-solver-decoupling-wave-a-implementation-plan.md`

---

**日期**：2026-03-31
