# Models 参数容器契约（统一容器 + 语义分层）

## 背景

在 Wave-E 统一扫描入口与 compat 退场阶段，需要把多模型参数输入口径统一，同时保证：

- 不为参数化模式（如 `pnjl_aniso`）再引入并行脚本/并行模型范式；
- 保持 AD/ID 友好（类型稳定、可微路径可追溯）；
- 保持模式判定、参数校验和错误提示清晰。

## 结论

采用**统一参数容器（NamedTuple）+ 语义分层**：

- 统一容器：对外接口、脚本、配置、序列化使用单一容器。
- 语义分层：容器内部按职责分层，而非全量扁平字段。

推荐顶层结构：

```julia
params = (
    control = (
        mode = :fixed_mu,
        T = 150.0,
        mu = 0.0,
        rho = nothing,
    ),
    background = (
        xi = 0.0,
        eB = 0.0,
        omega = 0.0,
        profile = :default,
    ),
    model = (
        kind = :PNJL,
        profile = :aniso,
    ),
    numerics = (
        p_num = 24,
        t_num = 8,
    ),
)
```

## 参数类别说明

### `control`（控制/约束参数）

用于定义“求解什么问题”：

- `mode`
- `T`
- `mu`
- `rho`（或其他约束目标）

`control` 直接参与模式判定与问题构造（`ProblemSpec` / `ConstraintModes`）。

### `background`（背景/介质参数）

用于定义“在何种背景下求解该问题”：

- `xi`
- `eB`
- `omega`
- `profile`

`background` 主要改变残差核、热力学核或分布函数形态。

### `model`（模型元信息）

- `kind`：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`
- `profile`：用于参数化模式（如 `:aniso`）

约定：`pnjl_aniso` 作为 `model.kind=:PNJL` + `model.profile=:aniso` + `background.xi!=0` 的组合表达，不引入新的独立 `model_kind`。

### `numerics`（数值参数）

- 网格精度、迭代控制、容差等数值设置。

## AD/ID 兼容要求

1. AD 角度允许统一容器；关键是避免不可微/不稳定操作。
2. `conditions`/residual 路径禁止把 Dual 参数强制转换为 `Float64`。
3. 参数读取要类型稳定，避免动态 key 拼接。
4. 分支逻辑基于 `mode`/`model.kind` 的离散派发，连续参数保持可微路径。

## 与现有模式系统的关系

保留当前“强类型 mode 分发”主干：

- `FixedMu/FixedRho/...` 继续作为主分发边界；
- 统一容器通过适配器映射到现有 `mode + kwargs` 接口；
- 不引入“全量动态字典驱动求解”的反模式。

## 错误口径与校验

建议按层校验并报错：

1. `model` 层：模型是否已注册、是否支持该扫描子命令；
2. `control` 层：`mode` 与约束目标是否一致；
3. `background` 层：参数是否在物理/数值可接受区间；
4. `numerics` 层：网格参数/容差是否合法。

## 最小迁移策略（Wave-E）

1. 先在统一 CLI 层接入 `params` 容器适配器；
2. 扫描链路先支持 `PNJL/NJL/RPNJL`，再扩展到 `PNJLMagnetic/Rotation/GasLiquid`；
3. 保持旧脚本 hard-deprecated 状态到阈值满足后再 remove/archive。
