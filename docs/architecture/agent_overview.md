# Agent 架构总览

供所有 agent（Codex、Claude Code、DeepSeek、Copilot）共用的模块与子系统结构说明。规则在 [AGENTS.md](../../AGENTS.md)，可执行命令在 [命令目录](../dev/agent_command_reference.md)，测试分层在 [testing_governance](../dev/testing_governance.md)。本页只回答「代码在哪里、按什么顺序串起来」。

## Include-driven 模块系统

没有 `src/PackageName.jl` 入口。模块直接在 `Main` 下 `include`，并自带重复加载守卫：

```julia
if !isdefined(Main, :ModuleName)
    include("src/path/to/module.jl")
end
```

`using` / `import` 放在模块顶部，编码与导入风格规则见 [AGENTS.md](../../AGENTS.md)。

## 模型体系（`src/models/`）

用户可见操作统一走 `src/models/entrypoints.jl`（`Models` 模块）；`src/pnjl/PNJL.jl` 是历史兼容背景，不作为入口。

类型层级（已按源码核对）：

```
AbstractQCDModel
├── AbstractNJLModel
│   ├── NJLModel
│   └── NJL2Model
├── AbstractPNJLModel
│   ├── PNJLModel
│   ├── PNJLMagneticModel
│   └── RPNJLModel
├── RotationModel
└── GasLiquidModel
```

每个具体模型按固定契约组织：`capabilities.jl`、`api.jl`、`adapters/entrypoint_adapter.jl`、`workflows/noop.jl`、`core/`。

## 求解器（`src/models/solver/`）

- 约束模式：`FixedMu`、`FixedRho`、`FixedEntropy`、`FixedSigma`、`FixedAsymmetricRho`。
- 调用链：`solve()` / `solve_constraint()` → `build_problem_spec()`（`ProblemSpec`）→ `ProblemSpec.forward_solve()` → `ConstraintSolver`（NLsolve 后端）或 `GapSolver`。
- `CandidateGovernance` 对多个 seed 策略的结果做评估并选出物理分支。
- 微分引擎（`src/models/solver/diff/`）：默认 `TaylorDiff`（`:auto` → `:taylordiff`）；gap 方程梯度用 `ForwardDiff`。`derivative_backend=:forwarddiff` 的 ImplicitFunction fallback **已下线**；底层 `solve_pnjl_with_flavor_mu_derivatives` 当前仅支持 `order=1`，更高阶张量留给 susceptibility 层（见 `src/models/solver/compat/ImplicitGapLegacy.jl`）。
- 结构不变量：非 fixedmu 模式使用 flattened joint-solve 语义；改动边界由 [AGENTS.md](../../AGENTS.md) 规定。

## 相图流水线（`src/models/phase/`）

- 主流程：`PhasePipeline.jl` / `ProductionPhasePipeline.jl` 串联各阶段。
- 阶段实现：`CEPDetector.jl`（含 spinodal 相关判定）、`CrossoverLine.jl`、`AdaptiveRhoRefinement.jl`、`RhoSupportRefinement.jl`、`PhaseGridConvergence.jl`。
- 产物：`PhaseArtifacts.jl` 写出 `phase_report.md` 与 `phase_summary.json`，读写细节在 `PhaseIO.jl`。

## RelaxTime / 输运（`src/relaxtime/`）

按依赖顺序：

`OneLoopIntegrals` → `PolarizationAniso`（含 `PolarizationCache`）→ `AFieldBuilder` → `EffectiveCouplings` → `MesonPropagator` → `TotalPropagator` → `MesonMass` → `MottTransition` → `DifferentialCrossSection` → `ScatteringAmplitude` → `AverageScatteringRate` → `RelaxationTime` → `TransportCoefficients`（η、ζ、σ）。

## 工作流引擎（`src/models/workflow_engine/`）

三类管线：

- **Workflow**：`solve_gap_and_transport`、`solve_gap_and_meson_point`
- **Scan**：`run_tmu_scan`、`run_trho_scan` 由 `ScanAdapter.jl` 绑定到 `scripts/models/run_unified_scan.jl`（该入口同时承接 `scan magnetic` 与 `workflow phase` 分支）
- **RelaxTime orchestrator**：统一的 `transport` / `cross-section` 编排

## 配置系统（`config/`）

TOML 加继承：模型参数 `config/models/<model>/default.toml`（可用 `unittest.toml` 覆盖）、共享物理常数 `config/physics/default.toml`、工作流档案 `config/workflows/`；由 `src/config/ConfigLoader.jl` 加载。
