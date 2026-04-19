# #81 solver/diff 试点集成设计（统一服务方案）

## 背景与目标

Issue #81 要求在真实链路完成导数层试点接入，并给出 Phase-2（是否回收到 `ProblemSpec`）的决策依据。已确认采用两条代表链路：

1. 脚本链：`scripts/pnjl/calculate_derivatives.jl`
2. 分析链：`scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl`

本设计采用 **B 方案（统一服务）**：抽出共享导数适配服务，由两条链路共同调用，降低重复逻辑并提升后续扩展可维护性。

非目标：

- 不在本 issue 做全仓迁移。
- 不在未验证前深改 `ProblemSpec` 主契约。

## 设计概览

### 新增共享模块

新增 `src/models/solver/diff/PilotAdapters.jl`，作为 #81 试点期专用统一服务层。

模块职责：

- 将脚本/分析层输入映射到 `ThermoDiffContext`。
- 统一组装 `ParamSpec` 与目标 `DiffTarget` 集。
- 提供稳定求值接口，返回可直接写 CSV 或用于分析的命名结果。
- 统一错误语义与 xi/spec_override 解析规则。

### 对外接口（初版）

1. `build_pilot_diff_context(result; mode, model, theta, spec_override=nothing, jacobian_backend=nothing)`
   - 对 `build_thermo_diff_context` 的薄封装。
   - 统一别名归一（`mu_fm`/`μ_fm`）和 `xi` 回退规则。

2. `eval_pilot_derivatives(ctx; target_names, param_names)`
   - 入参：`target_names::Vector{Symbol}`、`param_names::Vector{Symbol}`。
   - 内部：`DiffTarget` 与 `ParamSpec` 解析 + `jacobian(ctx, targets, params)`。
   - 出参：`NamedTuple`，包含：
     - `targets`
     - `params`
     - `jacobian::Matrix{Float64}`
     - `by_name::NamedTuple`（例如 `dP_dT`、`dP_dmu` 的便捷索引）

3. `eval_pilot_core_thermo(ctx)`（可选）
   - 统一读取 `pressure/energy/entropy/rho_norm`，用于分析链一致性记录。

> 说明：接口命名和粒度保持“试点可演进”，避免过早冻结到最终公共 API。

## 两条链路接入方案

### 1) `calculate_derivatives.jl`（脚本链）

替换策略：

- 保持 CLI、参数网格、CSV 列名与输出路径不变。
- 在每个 `(T, μ)` 网格点：
  - 先 `solve(...)` 获得平衡态；
  - 构建 `ThermoDiffContext`（经 pilot 适配层）；
  - 通过 `eval_pilot_derivatives(...)` 获取导数矩阵与命名值；
  - 映射回原有写盘字段。

兼容性约束：

- `DERIVATIVES_HEADER` 与 `BULK_HEADER` 不改。
- 未收敛或导数不可用时，沿用当前失败计数与 verbose 告警行为。

### 2) `t190_sigma_chain_decomposition.jl`（分析链）

替换策略：

- 不改变主物理流程（散射振幅、截面、blocking 链保持原样）。
- 在状态构建阶段额外接入 pilot 适配层，记录统一来源的导数/热力学响应量，供一致性与可读性评估。
- 输出新增列仅限“评估辅助列”；不影响原有关键结论列。

兼容性约束：

- 旧有 summary/detail 文件主列含义不变。
- 新增列采用显式前缀（如 `pilot_`）避免歧义。

## 错误处理与回退策略

- 参数合法性：复用 `ParamSpec` 与上下文构建时已有 `ArgumentError` 语义。
- 可恢复失败：脚本网格点失败不终止全局，按现有策略记录失败计数。
- 非可恢复失败：契约破坏（shape 不匹配、非有限值）直接抛错并失败退出（确保尽早暴露问题）。

## 测试与验证计划

### 新增/扩展单元测试

新增：`tests/unit/models/solver/test_solver_diff_pilot_adapters.jl`

覆盖点：

1. 目标映射：`target_names` 能正确解析并与 `DiffTarget` 对齐。
2. shape 契约：`N targets x P params` 返回矩阵维度正确。
3. xi 与 `mu_fm/μ_fm` 归一规则。
4. 错误路径：未知 target、未知 param、非有限值、空输入。

### 回归与验收命令

按 issue #81 DoD 执行：

```bash
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
```

补充 smoke：

```bash
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
```

## 性能与稳定性评估输出

新增评估文档（建议放在 `docs/dev/active/` 对应 #81 条目）包含：

1. 数值一致性：接入前后关键导数差异（rtol/atol）。
2. 运行开销：网格脚本与分析脚本的 wall-time 对比。
3. 可维护性：改动行数、重复逻辑减少情况、阅读入口数量。
4. Phase-2 建议：
   - 回收到 `ProblemSpec` 的最小回调集合；
   - 保留在 diff 外层模块的目标集合；
   - 下一阶段迁移模板（可执行步骤）。

## 风险与缓解

1. 试点链路代表性不足
   - 缓解：脚本链 + 分析链双样本；记录未覆盖边界。

2. 统一服务抽象过早固化
   - 缓解：以 pilot 前缀接口控制承诺范围，Phase-2 再决定公共化。

3. 额外求解调用导致性能回退
   - 缓解：优先复用已求解状态；在适配层内限制重复求解。

## 里程碑与完成定义

M1. 统一服务模块落地 + 单测通过。  
M2. 两条链路接入完成（脚本链 + 分析链）。  
M3. issue 指定验证命令通过 + smoke 通过。  
M4. 产出评估记录与 Phase-2 可执行建议。

DoD 对齐：满足 issue #81 全部勾选项。
