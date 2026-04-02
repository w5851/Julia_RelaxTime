# Models Dimension-Agnostic Plan-B Architecture Design

## Context

当前 `models` 主流程在求解器、约束组装与扫描链路中仍存在大量固定维度耦合（典型为 `5/3/8` 切片与 `SVector{5}/SVector{3}` 固定结构）。

这类耦合在模型扩展时持续放大维护心智负担：

1. 每次引入新模型或新状态变量时，需要同步修改多个“隐式维度假设点”；
2. 代码审查与回归定位高度依赖“记住各处切片语义”；
3. 入口统一后，内部结构仍非统一，长期形成“外统一、内分叉”。

在当前阶段，项目主要矛盾是维护复杂度，而非性能瓶颈。
因此本轮架构目标选择 Plan-B：优先统一动态向量管线，降低系统复杂度与演进成本。

## Goals

1. 建立单一 schema-driven 动态状态架构，替代固定维度切片编程模式。
2. 在 `solver -> constraint -> scans -> entrypoints` 主链路中移除旧架构依赖点。
3. 将兼容层定位为“迁移过渡资产”，最终目标是删除旧架构与兼容层，而非长期共存。
4. 保持对外行为稳定（入口参数与产物契约可兼容），并通过分层回归验证数值可接受。

## Implemented Model Inventory (Current)

基于 `src/models/factory.jl` 的已注册模型族：

- `:NJL`
- `:NJL2`
- `:PNJL`
- `:PNJLMagnetic`
- `:RPNJL`
- `:Rotation`
- `:GasLiquid`

说明：上述 7 个是当前“已实现并注册”的物理模型入口；Plan-B 最终态应覆盖这 7 个模型，不允许只覆盖扫描常用子集。

## Non-Goals

1. 不在本期重写物理模型公式与求解算法族。
2. 不将本期目标扩展为性能优化专项。
3. 不新增与 Plan-B 架构无关的模型功能。

## Target Architecture (End State)

### A1. Canonical State Schema

- 引入 `ModelStateSchema` 作为唯一状态布局定义源。
- 所有状态向量读写通过 schema 的 `flatten/unflatten` 与字段映射完成。
- 禁止在主链路手写固定索引（如 `x[1:5]`, `x[6:8]`）。

### A2. Dynamic Residual Assembly

- `Conditions` 与 `ImplicitSolver` 使用 schema 组装 residual 与约束。
- 约束模式（`fixed_mu/fixed_rho/...`）仅声明“语义约束”，不声明固定长度切片。

### A3. Unified Result Contract

- `SolverResult` 升级为 schema-aware 结构（状态字段动态化）。
- 扫描与工作流消费侧不再假设固定维度。

### A4. Decommission Legacy Architecture

- 旧固定维度路径仅允许在迁移阶段临时保留。
- 进入收口阶段后，删除旧路径与兼容桥接，形成唯一新架构主线。

### A5. Full-Model Coverage Requirement

- Plan-B 的“完成”定义按模型族全覆盖评估，而不是按扫描白名单评估。
- 对每个模型至少满足：
  1) schema 可构造；
  2) solver/constraint 主路径可运行；
  3) workflow 或 scan 至少一条主链路可验证。

## Migration Strategy

采用“过渡可用、终态清除”的两段式治理：

1. **Transition Window（迁移窗口）**
   - 允许最小兼容桥接存在，用于降低一次性切换风险；
   - 每个桥接点必须有迁移 owner 与移除条件。
2. **Finalization Window（收口窗口）**
   - 当 smoke/regression/governance 达标后，统一删除旧架构与 compat；
   - 不接受“永久双轨”作为完成状态。

## Acceptance Criteria

1. `src/models` 主链路不再存在固定 `5/3/8` 索引切片实现。
2. `scan/workflow/entrypoints` 在支持模型族上通过 smoke。
3. regression 代表点在可接受阈值内稳定。
4. 文档与 API 契约只描述新架构，不再把旧架构作为常规路径。
5. 旧架构与兼容层文件/接口完成删除或明确封存。
6. `NJL/NJL2/PNJL/PNJLMagnetic/RPNJL/Rotation/GasLiquid` 全部纳入覆盖矩阵并通过对应验证门禁。

## Verification Strategy

1. unit smoke：schema、conditions、solver 维度无关行为；
2. integration smoke：跨模型扫描与 workflow 主路径；
3. regression smoke：关键点收敛、残差、热力学量边界；
4. governance：docs 一致性与活动文档治理检查。

## Risks and Mitigations

1. Risk: 动态化后引入行为回归。
   - Mitigation: 每个迁移 chunk 均先写失败测试再实现，保留小步回归门。
2. Risk: 迁移期双轨拉长导致技术债固化。
   - Mitigation: 在计划中单列“旧架构移除”任务为硬门禁。
3. Risk: 外部使用方依赖旧字段结构。
   - Mitigation: 在过渡期提供明确 deprecation 提示与替代字段文档，并设置移除版本点。

## Deliverables

1. 本设计文档（spec）
2. 与 spec 对齐的实施计划（plan）
3. 最终仅保留新架构主线的实现与验证证据
