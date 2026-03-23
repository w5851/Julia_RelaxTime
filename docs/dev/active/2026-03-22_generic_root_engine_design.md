# 通用连续性求解引擎设计（零知识实施版）

## 1. 背景与目标

当前仓库已经在两条主线上实现了“连续性求解 + fallback + 质量判据”能力：

- Meson 混合态工作流（`src/models/workflows/MesonMassWorkflow.jl`）
- Gap 隐式求解器（`src/models/solver/ImplicitSolver.jl`）

问题在于两条链路的核心语义相近，但实现割裂，导致：

- 质量标签和诊断结构不统一，跨模块难聚合；
- fallback 和多初值策略重复实现；
- 领域约束（如物理性判据）与求解引擎边界不清。

本设计的目标是先定义最小通用接口和迁移路径，为后续 PR 分阶段落地，不在本轮直接改生产求解逻辑。

## 2. 范围与非范围

### 2.1 范围

- 统一 continuation / quality gate / fallback / diagnostics 的抽象边界。
- 覆盖变量维度 `2D/5D/8D` 与扫描变量（`T`、`mu`、`rho`、混合参数）的接口兼容。
- 支持单分支、双分支、多分支的追踪键语义。

### 2.2 非范围

- 不承诺“一套无配置 API 覆盖所有未来模型”。
- 不在本轮改写全部历史求解路径。
- 不在本轮实现混合态两阶段求根重构。

## 3. 设计原则

- 最小可迁移：只抽取两条现有主线共同语义。
- 领域适配外置：物理判据通过 adapter 注入，不嵌入引擎。
- 诊断优先兼容：新诊断结构可被现有字段映射，不破坏旧调用。
- 行为等价优先：Phase A/B 迁移门禁是“现有验证栈不退化”。

## 4. 核心抽象草案

以下为接口草案（以 Julia 伪代码描述，字段名按当前仓库命名习惯）。

```julia
struct RootProblemSpec{F,C}
    residual!::F               # (Fvec, x, ctx) -> nothing
    ctx::C                     # 领域上下文（温度、化学势、参数、模型句柄等）
    x_dim::Int                 # 未知量维度（2/5/8）
    branch_kind::Symbol        # :single | :paired | :multi
end

struct RootPolicy
    primary_method::Symbol     # :newton
    fallback_method::Symbol    # :trust_region
    use_fallback::Bool
    use_multiseed::Bool
    residual_norm_max::Float64
    require_converged::Bool
    diagnostics_level::Symbol  # :off | :basic | :full
end

struct ContinuationState
    seed_by_branch::Dict{Symbol, Vector{Float64}}  # branch_key => seed
    last_point::NamedTuple                         # 如 (T=..., mu=..., xi=...)
end

struct RootAttempt
    method::Symbol              # :newton / :trust_region / :multiseed
    seed_source::Symbol         # :continuation | :default | :multiseed
    converged::Bool
    residual_norm::Float64
    quality_tag::Symbol         # :good | :fallback | :degraded | :bad
end

struct RootDiagnostics
    selected_attempt::Int
    attempts::Vector{RootAttempt}
    branch_key::Symbol
    notes::Vector{String}
end
```

## 5. 最小执行入口草案

```julia
solve_root_with_policy(spec, x0; policy, continuation_state=nothing,
                       domain_quality=nothing, branch_key=nothing)

solve_root_continuation(scan_points, spec_factory;
                        policy, tracker=nothing,
                        domain_quality=nothing, branch_key=nothing)
```

约束说明：

- `spec_factory(scan_point)` 负责构造每个扫描点的 `RootProblemSpec`；
- `domain_quality` 可选，返回 `(ok::Bool, score::Float64, reason::Symbol)`；
- `branch_key` 可选，返回分支身份键（如 `:eta_minus`、`:eta_plus`、`:default`）。

## 6. Adapter 接口草案

### 6.1 残差接口

- `residual!(F, x, ctx)`
  - 输入：状态 `x` 与上下文 `ctx`
  - 输出：原地写入残差 `F`

### 6.2 领域质量接口（可选）

- `domain_quality(x, meta)`
  - 典型用于物理性判据（有限性、正性、相区约束）

### 6.3 分支键接口（可选）

- `branch_key(x, ctx)`
  - 单分支可固定 `:default`
  - 混合态可返回 `:eta_minus/:eta_plus`、`:sigma_minus/:sigma_plus`

## 7. 统一质量门禁与 fallback 流程

统一判据建议：

1. `converged == true`（若 policy 要求）
2. `residual_norm <= residual_norm_max`
3. `domain_quality == ok`（若提供）

统一流程建议：

1. 用 continuation seed（若有）跑 `primary_method`
2. 失败则同 seed 跑 `fallback_method`
3. 仍失败且启用多初值，则跑多 seed 候选
4. 候选由统一 selector 选择（默认优先 quality，再按 residual，再按领域分数）

## 8. 统一诊断结构

建议质量标签枚举固定为：

- `:good`：一次主方法通过
- `:fallback`：经回退方法得到可接受解
- `:degraded`：仅满足部分条件（如收敛但残差偏大）
- `:bad`：无可接受解

`attempt trace` 统一包含：

- `method`
- `seed_source`
- `converged`
- `residual_norm`
- `quality_tag`

## 9. 迁移映射（与现有模块对齐）

### 9.1 Phase A：Meson 先迁

- 保留在 `MesonMassWorkflow.jl`：
  - 混合态标签语义与输出格式；
  - 阈值/gap 相关领域结果。
- 抽到 `GenericRootEngine`：
  - 尝试链（primary/fallback/multiseed）；
  - quality gate 与 attempt trace。

### 9.2 Phase B：Gap 接入

- 将 `ImplicitSolver` 中 `_nlsolve_with_tr_fallback` 的共性逻辑改为引擎调用。
- Gap 侧特有判据（物理性、thermo finite）保留在 adapter。

### 9.3 Phase C：诊断统一

- workflow 与 solver 共享同一 `RootDiagnostics` 字段语义。
- 对外保留兼容层，避免直接破坏既有调用方。

## 10. 验证与门禁建议

迁移阶段最小门禁：

- Meson 相关：
  - `tests/validation/relaxtime/test_mixed_meson_root_quality_continuation.jl`
  - legacy/literature 关键验证栈
- Gap 相关：
  - FixedMu / FixedRho / FixedAsymmetricRho 对应回归集

通过标准：

- 行为不退化（通过率、质量标签分布不显著恶化）；
- 诊断字段可对齐旧结果；
- 无新增非有限/负质量等物理非法解泄漏。

## 11. 风险与应对

- 抽象过度风险：先只覆盖 Meson+Gap 两条主线，不提前泛化。
- 相变附近回归风险：维护异常点清单作为迁移门禁。
- 性能回退风险：`diagnostics_level=:basic` 作为默认轻量路径。

## 12. 下一 PR 实施建议（可直接照此拆分）

1. 新增 `src/models/solver/GenericRootEngine.jl` 最小骨架与单元测试。
2. 只迁 Meson 共性求解链路，保持输出接口兼容。
3. 验证全绿后再迁 Gap，并补齐回归门禁。
4. 若导出字段成为稳定 API，再更新 `docs/api/`。
