# Solver 三层模型最小接口清单（实现基线）

更新日期：2026-04-06

## 1. 目的

- 固化“语义外壳 + Schema 桥接 + 向量内核”三层模型的最小接口。
- 作为后续任务拆解与实现验收的基线文档。

## 2. 三层模型

### 2.1 语义外壳层（NamedTuple 面向业务）

建议接口：

```julia
solve_named(model, spec, theta_named; kwargs...) -> SolveNamedResult
derive_named(model, spec, theta_named; kwargs...) -> DerivativeNamedResult
```

职责：

- 参数可读性与语义表达。
- 默认值补全与用户输入校验。
- 错误信息语义化（避免暴露内核细节）。

### 2.2 Schema 桥接层（NamedTuple <-> Vector）

建议接口：

```julia
schema_for(model_kind, spec_tag) -> VarSchema

named_to_vec(x_named, schema, :x)
vec_to_named(x_vec, schema, :x)

named_to_vec(theta_named, schema, :theta)
vec_to_named(theta_vec, schema, :theta)

validate_schema(schema, dims)
```

职责：

- 固定顺序映射（键顺序唯一且稳定）。
- NamedTuple 与向量之间无语义损失重排。
- 类型保持（禁止在桥接层降精度强转）。

### 2.3 向量内核层（Solver/AD 共用）

建议接口：

```julia
solve_vec(model, spec, theta_vec; primary_strategy, kwargs...) -> SolvePointResult
derive_vec(model, spec, theta_vec; x_star=nothing, method=:implicit_ad, kwargs...) -> DerivativeResult
residual_vec!(F, x_vec, theta_vec, cfg)
```

职责：

- 统一数值求解执行。
- 统一隐式微分执行。
- 统一诊断输出格式。

## 3. 最小统一契约对象

```julia
ConstraintDims(x_dim, theta_dim)

ConstraintSpec(
    tag,
    dims,
    residual_named!,
    unpack_solution,
    postprocess,
    hard_rules,
    selector,
    extra_constraints,
)

PolicyProfile(
    name,
    root_policy,
    primary_strategy,
    residual_norm_max,
    allow_legacy_fallback,
)
```

说明：

- `primary_strategy` 应统一包含 method、多 seed、fallback，不做并列分散配置。

## 4. 必做约束（防回归）

- 内核热路径禁止 `Float64(...)` 强转（保持 Dual 传播）。
- `schema` 顺序不可变，且由单一注册源提供。
- `primary_strategy` 作为统一策略入口，避免策略碎片化。
- 导数层仅复用 `residual_vec!` 契约，不读取 solver 迭代内部状态。

## 5. 结果对象（建议）

```julia
struct SolvePointResult
    x_star
    converged::Bool
    residual_norm
    diagnostics
end

struct DerivativeResult
    x_star
    dx_dtheta
    method::Symbol
    diagnostics
end
```

## 6. 验收口径（最小）

- 同一输入下：`solve_named` 与 `solve_vec` 的解和残差一致（容差内）。
- 同一输入下：`derive_named` 与 `derive_vec` 的一阶导一致（容差内）。
- 桥接层往返一致：`named -> vec -> named` 保持键值与顺序语义一致。
- AD 路径无类型破坏：Dual 能完整通过 residual 链路。

## 7. 适用边界

- 该清单用于 solver 核心与导数解耦落地。
- 不覆盖具体物理模型公式实现与 workflow 调度细节。

## 8. 实现准则：零重复转换（性能与可微分约束）

### 8.1 一次绑定，不在迭代里查表

- 进入求解前基于 schema 生成 `ResidualKernel`（预绑定索引和 `residual_core!`）。
- `residual!` 热路径仅做数组读写，不做键查找。

### 8.2 热路径禁止命名结构分配

- 允许：固定索引访问、`@views`、原地写 `F`。
- 禁止：`NamedTuple(...)`、`Dict`、`pairs`、字符串拼接等动态分配。

### 8.3 边界转换固定化

- `solve_named` 入口：`theta_named -> theta_vec`（一次）。
- `solve_named` 可选入口：`x_named -> x0_vec`（一次）。
- 内核迭代：纯 `vec`。
- 结果出口：`x_star_vec -> x_star_named`（一次）。

### 8.4 AD/Implicit 路径同契约

- `ForwardDiff` 与隐式微分统一调用 `residual_vec!`。
- 禁止在 Jacobian 回调中 `vec -> named -> vec` 往返转换。

### 8.5 类型纪律

- 适配函数使用 `AbstractVector{T}` 保持元素类型泛型。
- 禁止在导数路径中 `Vector{Float64}(...)` 或 `Float64(...)` 强转。
- 仅在最终用户输出层允许数值归一化到 `Float64`。

### 8.6 建议监控指标

- `residual!` 单点评估分配应接近 0（或固定极低）。
- AD 打开时无 Dual 类型破坏报错。
- 每个请求 `named->vec` 与 `vec->named` 各一次（可用 debug counter 验证）。

该准则目标是将语义转换前后移，保持数值循环纯向量化。
