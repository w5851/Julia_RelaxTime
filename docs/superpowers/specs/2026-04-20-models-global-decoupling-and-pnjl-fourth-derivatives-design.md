# models 全域解耦与 PNJL 四阶导能力吸收设计

- 日期：2026-04-20
- 状态：draft（已用于创建 issue）
- 关联 issue：
  - https://github.com/w5851/Julia_RelaxTime/issues/89
  - https://github.com/w5851/Julia_RelaxTime/issues/90

## 1. 背景与目标

本次需求包含两条相互耦合但应分阶段推进的主线：

1. `src/models` 全域能力解耦与重复实现治理。
2. 在统一边界上吸收 PNJL 四阶导数能力（参考外部仓，不迁移代码）。

采用“基础先行”策略：先完成 models 全域边界收敛（Issue #89），再实施四阶导能力扩展（Issue #90）。

## 2. 现状核查结论（models 全量结构）

### 2.1 结构分层（现状）

- 核心物理域：`src/models/pnjl_physics/`, `src/models/njl/`, `src/models/rpnjl/`, `src/models/variants/*/core`
- 求解域：`src/models/solver/{runtime,spec,orchestrator,governance,diff,compat,api}`
- 导数域：`src/models/derivatives/*`，以及 `src/models/solver/diff/*`、`src/models/solver/compat/ImplicitGapLegacy.jl`
- 流程域：`src/models/scans/*`, `src/models/workflow/*`, `src/models/workflows/*`, `src/models/phase/*`
- 入口与装配域：`src/models/Models.jl`, `src/models/entrypoints.jl`, `src/models/factory.jl`, `src/models/precompile/*`

### 2.2 重复与边界混叠（首轮）

- 导数能力分布在三处：
  - `src/models/derivatives/ThermoDerivatives.jl`
  - `src/models/solver/diff/{Targets.jl,PilotAdapters.jl}`
  - `src/models/solver/compat/ImplicitGapLegacy.jl`
- 内部 `Main.Models.*` 反向调用较多，显示子域间依赖边界仍偏“总线式”。
- workflow/scans/phase 与 solver/derivatives 的装配点分散，增加变更影响半径与认知成本。

## 3. 设计原则

1. 单一权威实现：每类能力只保留一个权威层。
2. 适配层显式化：保留兼容/桥接能力，但不得成为并行主实现。
3. 稳定入口不破坏：`Models` 与 `entrypoints` 合同保持稳定。
4. 语义不漂移：保持 mixed-meson governance 与 non-fixedmu 联合求解语义。
5. 参考不迁移：外部仓仅用于公式和行为对照。

## 4. 分阶段方案

### Phase A（Issue #89）：models 全域边界收敛与去重

#### A1 规划产物

- `models` 能力地图：模块 -> 职责 -> 上游/下游依赖 -> 权威文件
- 去重裁决表：每类能力的 authority/adapter/compat 划分
- 目录边界重构方案：逻辑边界优先、迁移顺序、回滚点
- 治理守护项：防止多点实现回流

#### A2 收敛目标

- 导数链路由单一权威层统一：
  - 计算权威层（derivatives）
  - solver/diff 仅做目标映射与 jacobian 组装
  - compat 层仅保留历史桥接
- workflow/scans/phase 通过稳定入口契约使用能力，不新增旁路。

#### A3 DoD

- 能力地图与裁决表文档化
- 每类核心能力唯一权威
- 目录按能力域可读
- 现有入口合同与回归测试通过

### Phase B（Issue #90）：PNJL 四阶导能力吸收（依赖 Phase A）

#### B1 能力边界

- 主线落地：PNJL 四阶导数与相关比值（如 `kappa4/kappa2`）
- 保留跨模型扩展点：Rotation/Gas_Liquid 先定义抽象接口，不强制同批实现

#### B2 参考源（仅行为/公式）

- `D:/Desktop/PNJL_Simulation/src/Function_PNJL.jl`
- `D:/Desktop/PNJL_Simulation/src/Gas_Liquid/Advanced_ForwardDiff.jl`
- `D:/Desktop/PNJL_Simulation/src/Rotation/Advanced_ForwardDiff.jl`

#### B3 验证策略

- 与现有低阶（1-3阶）语义、缩放、单位保持一致
- 固定点/采样点回归基线
- 与外部参考实现及可用文献结果对照

#### B4 DoD

- 四阶能力通过主线 API 可调用
- 缩放与单位定义明确且一致
- 回归基线接入治理
- 扩展接口保留且不破坏现有能力

## 5. 风险与缓解

- 风险：边界收敛期间行为偏移
  - 缓解：契约测试先行，迁移分批进行
- 风险：四阶 AD 嵌套导致数值/性能波动
  - 缓解：固定点回归 + 分层 profiling + 混合容差治理
- 风险：兼容层清理过早
  - 缓解：阶段性 deprecate 与迁移映射

## 6. 本轮决策记录

- issue 创建方式：直接创建真实 issue
- 拆分策略：方案 B（基础先行）
- 四阶 issue 粒度：PNJL 主线 + 跨模型可复用抽象扩展点
- 跨仓策略：参考不迁移

## 7. 后续衔接

- 本设计已落地为 #89 与 #90。
- 用户确认本 spec 后，再进入 implementation plan 阶段。
