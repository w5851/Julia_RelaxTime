# Solver ConstraintSpec 接口草案（软去 Mode）

更新日期：2026-04-06

## 1. 背景与目标

- 背景：当前不同 `ConstraintMode` 在实现上存在架构分支，主链与兼容链并存。
- 目标：将 mode 从“架构分支控制器”降级为“约束配置输入”，把差异统一表达为 `extra_constraints`。
- 原则：保持对外 API 稳定，采用软迁移，不一次性破坏旧调用。

## 2. 设计原则

- 统一求解内核：内核仅关心 residual、policy、candidate governance。
- 统一差异表达：`FixedMu/FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 差异仅在约束配置。
- 统一结果契约：`SolverResult + Diagnostics` 字段语义跨模式保持一致。
- 兼容优先：旧 mode 保留为 `ConstraintSpec` 构造器。

## 3. 核心数据结构（草案）

```julia
struct ConstraintDims
    x_dim::Int
    theta_dim::Int
end

struct ConstraintSpec{R,U,P,H,S,E}
    tag::Symbol
    dims::ConstraintDims
    base_residual!::R
    unpack_solution::U
    postprocess::P
    hard_rules::H
    selector::S
    extra_constraints::E
end
```

建议输入语义：

```julia
x_named::NamedTuple
theta_named::NamedTuple
```

执行层适配：

```julia
residual_vec!(F, x_vec, theta_vec, cfg)
```

说明：对外使用命名输入以承载语义，进入内核前一次映射到 `theta_vec`；数值内核统一使用向量以适配 `NLsolve`。

字段说明：

- `tag`：仅用于日志、诊断、兼容层标记，不应驱动主架构分支。
- `dims`：统一数值维度契约，只约束 `x` 与 `theta`，避免散落硬编码。
- `base_residual!`：基础方程（如平衡条件）。
- `extra_constraints`：模式差异唯一入口。

语义边界：

- solver 内核只关心 `x_dim/theta_dim` 与残差一致性。
- 物理语义（如 state/mu 切片）通过回调表达，不作为维度字段。
- 推荐回调：`unpack_solution(x, cfg) -> domain_state`、`pack_seed(domain_state, cfg) -> x0`。

## 4. 目标 API（Spec-First）

```julia
solve_constraint(model, spec::ConstraintSpec, theta; kwargs...) -> SolverResult
solve_multi(model, spec::ConstraintSpec, theta; kwargs...) -> SolverResult
```

兼容 API（迁移期保留）：

```julia
solve_constraint(model, mode::ConstraintMode, ...)
```

兼容语义：

- 旧 mode 接口内部只负责 `mode -> spec`。
- 实际求解一律进入 spec-first 主链。

## 5. Residual 组合协议

内核只认统一协议：

```julia
residual_vec!(F, x_vec, theta_vec, cfg) = [base_equations; extra_constraints...]
```

约束：

- mode 差异不得在主调度分支中硬编码。
- 新增模式时，优先新增约束配置而非新增调度分支。

## 6. Policy 解耦草案

```julia
struct PolicyProfile
    name::Symbol
    root_policy::RootPolicy
    primary_strategy
    residual_norm_max::Float64
    allow_legacy_fallback::Bool
end
```

建议：

- `ConstraintSpec` 不承载可变策略细节，只引用 policy profile。
- 阈值与 fallback 规则集中在 policy center 管理。
- `primary_strategy` 包含 method、multi-seed 与 fallback 组合，不再把 multi-seed 作为并列能力。

## 7. 兼容迁移路线（软去 Mode）

### Phase A：引入 Spec，不改变外部调用

- 新增 `build_constraint_spec(mode, model, theta; kwargs...)`。
- 旧入口继续可用。

### Phase B：内核切换为 Spec-First

- `solve_constraint(mode,...)` 入口内部立即转 spec。
- 主链只处理 `ConstraintSpec`。

### Phase C：文档与治理收口

- 文档改为 spec-first。
- mode 标记为 legacy-friendly wrapper。

## 8. 最小文件改动建议

- `src/models/solver/ProblemSpec.jl`
  - 新增 `ConstraintSpec` 定义与构建器。
- `src/models/solver/Solver.jl`
  - 新增 `solve_constraint(..., spec::ConstraintSpec, ...)`。
  - 旧 mode 入口改为构建 spec 后转调。
- `src/models/solver/Conditions.jl`
  - 收敛为基础残差 + 额外约束拼装。
- `src/models/constraint_solver.jl`
  - 聚焦治理能力（rules/selector/diagnostics），剥离 mode 业务分支。

## 9. 风险与护栏

- 风险：`FixedMu` 迁移后收敛路径变化。
- 护栏：A/B 双通路比对（mode path vs spec path），对比 `residual_norm`、`pressure`、`selected_reason`、`candidate_count`。
- 准入：仅在回归与集成用例一致时切主。

## 10. 决策待确认项

- `ConstraintSpec` 是否直接复用 `ProblemSpec`，或并行存在后再合并。
- `PolicyProfile` 是否允许 workflow 层覆盖，覆盖边界到哪里。
- diagnostics 最小稳定字段集合（跨 mode 强约束）。

## 11. 推荐实现模式：Schema 注册表（避免模型分支散落）

建议采用“通用转换函数 + schema 注册表”模式：

- 通用层只实现一套：`x_to_vec/vec_to_x/theta_to_vec/vec_to_theta`。
- 模型差异通过 `schema` 声明表达，不在转换函数里写模型分支。
- 新模型或新约束只需新增 schema 注册，不改通用适配逻辑。

示例：

```julia
struct VarSchema{NX,NT}
    x_keys::NTuple{NX,Symbol}
    theta_keys::NTuple{NT,Symbol}
end

struct SchemaRegistry
    table::Dict{Tuple{Symbol,Symbol},VarSchema}
end

function register_schema!(reg::SchemaRegistry, model_kind::Symbol, spec_tag::Symbol, schema::VarSchema)
    reg.table[(model_kind, spec_tag)] = schema
    return reg
end

function schema_for(reg::SchemaRegistry, model_kind::Symbol, spec_tag::Symbol)
    key = (model_kind, spec_tag)
    haskey(reg.table, key) || throw(ArgumentError("missing schema for $(key)"))
    return reg.table[key]
end
```

初始化示意：

```julia
reg = SchemaRegistry(Dict{Tuple{Symbol,Symbol},VarSchema}())

register_schema!(reg, :PNJL, :fixed_mu,
    VarSchema((:phi_u,:phi_d,:phi_s,:Phi,:PhiBar), (:T,:mu)))

register_schema!(reg, :PNJL, :fixed_rho,
    VarSchema((:phi_u,:phi_d,:phi_s,:Phi,:PhiBar,:mu_u,:mu_d,:mu_s), (:T,:rho_target)))
```

内核调用路径建议：

1. `schema = schema_for(reg, model_kind, spec.tag)`
2. 通过通用适配函数完成 named <-> vec 转换
3. 调用 `residual_named(...)` 计算语义残差
4. 回写到 `F` 并交给 `NLsolve`

该模式直接支持“mode 降格为配置输入”的目标态，并降低后续维护成本。
