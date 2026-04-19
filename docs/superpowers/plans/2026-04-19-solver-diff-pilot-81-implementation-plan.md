# Solver/Diff Pilot #81 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 `calculate_derivatives` 与 `t190_sigma_chain_decomposition` 两条链路落地统一导数服务接入，完成 #81 的可维护性/稳定性/性能评估与 Phase-2 决策产出。

**Architecture:** 在 `src/models/solver/diff/` 新增 `PilotAdapters.jl` 作为试点统一服务，集中处理 `ThermoDiffContext` 构建、目标解析、Jacobian 求值与命名结果映射。两条链路只保留“参数采集 + 调用服务 + 写盘/分析”，避免再次散落 ad-hoc 导数调用。通过新增单测和链路验证命令保证契约一致与回归可控。

**Tech Stack:** Julia 1.10+、Models solver/diff 模块、现有 tests 分层入口、GitHub issue #81 DoD 验证命令。

---

## File Structure

- Create: `src/models/solver/diff/PilotAdapters.jl`
  - 责任：#81 试点统一导数服务入口，封装 context 构建、target/param 解析、jacobian 求值、命名映射。
- Modify: `src/models/solver/topology.jl`
  - 责任：将 `PilotAdapters.jl` 纳入 solver/diff include 拓扑。
- Modify: `src/models/Models.jl`
  - 责任：导出试点服务 API（pilot 命名空间）。
- Create: `tests/unit/models/solver/test_solver_diff_pilot_adapters.jl`
  - 责任：覆盖 pilot 服务核心契约（映射、shape、xi/μ 别名、错误路径）。
- Modify: `scripts/pnjl/calculate_derivatives.jl`
  - 责任：脚本链接入统一服务，保持 CSV 头与 CLI 不变。
- Modify: `scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl`
  - 责任：分析链接入统一服务，新增 `pilot_` 辅助评估列。
- Create: `docs/dev/active/issue-81-solver-diff-pilot-evaluation.md`
  - 责任：记录接入前后稳定性/性能/可维护性对比与 Phase-2 建议。

---

### Task 1: Pilot 服务单测骨架（先红）

**Files:**
- Create: `tests/unit/models/solver/test_solver_diff_pilot_adapters.jl`

- [ ] **Step 1: Write the failing test file**

```julia
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver diff pilot adapters contract" begin
    T_fm = 100.0 / 197.327
    mu_fm = 0.0
    model = Models.create_model(:PNJL)
    mode = Models.FixedMu()
    result = Models.solve(model, mode, T_fm, mu_fm; p_num=8, t_num=4, residual_norm_max=1e-6)

    @testset "build_pilot_diff_context supports mu aliases" begin
        ctx = Models.build_pilot_diff_context(result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, μ_fm=mu_fm, xi=0.0),
            spec_override=(p_num=8, t_num=4, residual_norm_max=1e-6),
        )
        @test hasproperty(ctx.theta, :mu_fm)
        @test ctx.theta.mu_fm == mu_fm
    end

    @testset "eval_pilot_derivatives returns named matrix payload" begin
        ctx = Models.build_pilot_diff_context(result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
            spec_override=(p_num=8, t_num=4, residual_norm_max=1e-6),
        )
        payload = Models.eval_pilot_derivatives(ctx;
            target_names=[:pressure, :dP_dT],
            param_names=[:T_fm, :mu_fm],
        )
        @test size(payload.jacobian) == (2, 2)
        @test all(isfinite, payload.jacobian)
        @test hasproperty(payload.by_name, :pressure__d_dT_fm)
    end

    @testset "eval_pilot_derivatives rejects unknown target" begin
        ctx = Models.build_pilot_diff_context(result;
            mode=mode,
            model=model,
            theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
        )
        @test_throws ArgumentError Models.eval_pilot_derivatives(ctx;
            target_names=[:unknown_target],
            param_names=[:T_fm],
        )
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_diff_pilot_adapters.jl")'`  
Expected: FAIL，提示 `build_pilot_diff_context` / `eval_pilot_derivatives` 未定义。

- [ ] **Step 3: Commit**

```bash
git add tests/unit/models/solver/test_solver_diff_pilot_adapters.jl
git commit -m "test(models/solver): add failing pilot adapter contract tests"
```

---

### Task 2: 实现 PilotAdapters 最小可用接口（转绿）

**Files:**
- Create: `src/models/solver/diff/PilotAdapters.jl`
- Modify: `src/models/solver/topology.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/solver/test_solver_diff_pilot_adapters.jl`

- [ ] **Step 1: Write the minimal implementation**

```julia
const _PILOT_PARAM_ALIAS = Dict{Symbol,Symbol}(
    :μ_fm => :mu_fm,
)

@inline function _normalize_pilot_param_name(name::Symbol)
    return get(_PILOT_PARAM_ALIAS, name, name)
end

@inline function _pilot_theta(theta::NamedTuple)
    if hasproperty(theta, :mu_fm) && hasproperty(theta, :μ_fm)
        throw(ArgumentError("theta cannot contain both mu_fm and μ_fm"))
    end
    if hasproperty(theta, :μ_fm)
        return merge(Base.structdiff(theta, NamedTuple{(:μ_fm,)}), (mu_fm=getproperty(theta, :μ_fm),))
    end
    return theta
end

@inline function build_pilot_diff_context(result; mode, model, theta::NamedTuple, spec_override=nothing, jacobian_backend=nothing)
    θ = _pilot_theta(theta)
    if jacobian_backend === nothing
        return build_thermo_diff_context(result; mode=mode, model=model, theta=θ, spec_override=spec_override)
    end
    return build_thermo_diff_context(result; mode=mode, model=model, theta=θ, spec_override=spec_override, jacobian_backend=jacobian_backend)
end

@inline function _pilot_targets(target_names::AbstractVector{Symbol})
    isempty(target_names) && throw(ArgumentError("target_names must be non-empty"))
    return [diff_target(name) for name in target_names]
end

@inline function _pilot_params(param_names::AbstractVector{Symbol})
    isempty(param_names) && throw(ArgumentError("param_names must be non-empty"))
    normalized = Symbol[_normalize_pilot_param_name(name) for name in param_names]
    return ParamSpec(normalized)
end

@inline function _pilot_derivative_key(target::Symbol, param::Symbol)
    p = param === :T_fm ? :dT_fm : (param === :mu_fm ? :dmu_fm : Symbol(:d_, param))
    return Symbol(target, :__ , p)
end

@inline function eval_pilot_derivatives(ctx::ThermoDiffContext; target_names::Vector{Symbol}, param_names::Vector{Symbol})
    targets = _pilot_targets(target_names)
    params = _pilot_params(param_names)
    J = jacobian(ctx, targets, params)
    kv = Pair{Symbol,Float64}[]
    for (i, tname) in enumerate(target_names)
        for (j, pname) in enumerate(params.names)
            push!(kv, _pilot_derivative_key(tname, pname) => J[i, j])
        end
    end
    return (
        targets=target_names,
        params=params.names,
        jacobian=J,
        by_name=NamedTuple(kv),
    )
end
```

- [ ] **Step 2: Wire module include and exports**

在 `src/models/solver/topology.jl` 增加：

```julia
include(joinpath(@__DIR__, "diff", "PilotAdapters.jl"))
```

在 `src/models/Models.jl` 增加导出：

```julia
export build_pilot_diff_context, eval_pilot_derivatives
```

- [ ] **Step 3: Run pilot adapter test to verify it passes**

Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_diff_pilot_adapters.jl")'`  
Expected: PASS。

- [ ] **Step 4: Commit**

```bash
git add src/models/solver/diff/PilotAdapters.jl src/models/solver/topology.jl src/models/Models.jl tests/unit/models/solver/test_solver_diff_pilot_adapters.jl
git commit -m "feat(models/solver): add pilot diff adapters for issue 81"
```

---

### Task 3: 扩展单测覆盖错误路径与 shape 契约

**Files:**
- Modify: `tests/unit/models/solver/test_solver_diff_pilot_adapters.jl`

- [ ] **Step 1: Add failing edge-case tests**

```julia
@testset "pilot adapters edge cases" begin
    T_fm = 100.0 / 197.327
    mu_fm = 0.0
    model = Models.create_model(:PNJL)
    mode = Models.FixedMu()
    result = Models.solve(model, mode, T_fm, mu_fm; p_num=8, t_num=4, residual_norm_max=1e-6)

    ctx = Models.build_pilot_diff_context(result;
        mode=mode,
        model=model,
        theta=(T_fm=T_fm, mu_fm=mu_fm, xi=0.0),
    )

    @test_throws ArgumentError Models.eval_pilot_derivatives(ctx; target_names=Symbol[], param_names=[:T_fm])
    @test_throws ArgumentError Models.eval_pilot_derivatives(ctx; target_names=[:pressure], param_names=Symbol[])
    @test_throws ArgumentError Models.eval_pilot_derivatives(ctx; target_names=[:pressure], param_names=[:bad_param])

    payload = Models.eval_pilot_derivatives(ctx;
        target_names=[:pressure, :entropy, :rho_norm],
        param_names=[:T_fm],
    )
    @test size(payload.jacobian) == (3, 1)
end
```

- [ ] **Step 2: Run test to verify new tests fail first**

Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_diff_pilot_adapters.jl")'`  
Expected: FAIL（至少一个新增断言失败）。

- [ ] **Step 3: Implement minimal fixes if failing**

若失败源于 `NamedTuple(kv)` 构造问题，改为：

```julia
keys_tuple = Tuple(first.(kv))
vals_tuple = Tuple(last.(kv))
by_name = NamedTuple{keys_tuple}(vals_tuple)
```

若失败源于空输入未拦截，补 `isempty(...)` 判定抛 `ArgumentError`。

- [ ] **Step 4: Re-run test to verify pass**

Run: `julia --project=. -e 'include("tests/unit/models/solver/test_solver_diff_pilot_adapters.jl")'`  
Expected: PASS。

- [ ] **Step 5: Commit**

```bash
git add tests/unit/models/solver/test_solver_diff_pilot_adapters.jl src/models/solver/diff/PilotAdapters.jl
git commit -m "test(models/solver): cover pilot adapter edge and shape contracts"
```

---

### Task 4: 脚本链 `calculate_derivatives.jl` 接入统一服务

**Files:**
- Modify: `scripts/pnjl/calculate_derivatives.jl`
- Test: `tests/unit/models/solver/test_solver_diff_pilot_adapters.jl`

- [ ] **Step 1: Write failing integration-style assertion (script-level helper)**

在脚本中新增可测 helper（先写调用断言），例如：

```julia
function _pilot_eval_point(model, mode, T_fm, μ_fm, xi, p_num, t_num)
    result = Models.solve(model, mode, T_fm, μ_fm; xi=xi, p_num=p_num, t_num=t_num)
    ctx = Models.build_pilot_diff_context(result;
        mode=mode,
        model=model,
        theta=(T_fm=T_fm, mu_fm=μ_fm, xi=xi),
        spec_override=(p_num=p_num, t_num=t_num, residual_norm_max=1e-6),
    )
    return Models.eval_pilot_derivatives(ctx;
        target_names=[:pressure, :energy, :entropy, :rho_norm, :dP_dT, :dP_dmu],
        param_names=[:T_fm, :mu_fm],
    )
end
```

- [ ] **Step 2: Run focused unit tests before implementation wiring**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/solver/test_solver_diff_pilot_adapters.jl"; include("tests/unit/runtests.jl")'`  
Expected: PASS（确保服务可被脚本调用）。

- [ ] **Step 3: Replace ad-hoc derivative call path with pilot payload mapping**

在 `step1_mass_derivatives` 中将旧 `thermo_derivatives(...)` 调用改为 `_pilot_eval_point(...)` + jacobian 映射：

```julia
payload = _pilot_eval_point(model, mode, T_fm, μ_fm, config.xi, config.p_num, config.t_num)
J = payload.jacobian

dP_dT = J[1, 1]
dP_dmu = J[1, 2]
dEpsilon_dT = J[2, 1]
dEpsilon_dmu = J[2, 2]
dn_dT = J[4, 1]
dn_dmu = J[4, 2]
```

并保持原 CSV 头与字段顺序完全不变。

- [ ] **Step 4: Run script smoke on tiny grid**

Run: `julia --project=. scripts/pnjl/calculate_derivatives.jl --T_min=150 --T_max=150 --T_step=10 --mu_min=0 --mu_max=0 --mu_step=50 --p_num=8 --t_num=4 --output_dir=data/outputs/results/pnjl/derivatives_pilot_smoke`  
Expected: 生成两份 CSV 且字段头不变。

- [ ] **Step 5: Commit**

```bash
git add scripts/pnjl/calculate_derivatives.jl
git commit -m "refactor(scripts/pnjl): route derivative script through pilot diff adapters"
```

---

### Task 5: 分析链 `t190_sigma_chain_decomposition.jl` 接入统一服务

**Files:**
- Modify: `scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl`

- [ ] **Step 1: Add failing assertion for pilot helper output keys**

在分析脚本增加 helper：

```julia
function pilot_eval_state_derivatives(result, T_fm, mu_fm, xi)
    model = Models.create_model(:PNJL)
    mode = FixedMu()
    ctx = Models.build_pilot_diff_context(result;
        mode=mode,
        model=model,
        theta=(T_fm=T_fm, mu_fm=mu_fm, xi=xi),
        spec_override=(p_num=12, t_num=6, residual_norm_max=1e-6),
    )
    return Models.eval_pilot_derivatives(ctx;
        target_names=[:pressure, :entropy],
        param_names=[:T_fm],
    )
end
```

- [ ] **Step 2: Run script syntax check**

Run: `julia --project=. -e 'include("scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl")'`  
Expected: 若路径数据不存在可报运行时输出路径错误，但不应有语法错误。

- [ ] **Step 3: Wire pilot metrics into CSV as auxiliary columns**

在 detail/header 增加辅助列（示例）：

```julia
println(io, "...,pilot_dpressure_dT,pilot_dentropy_dT")
```

在 summary/header 增加辅助列（示例）：

```julia
println(io_sum, "...,pilot_dpressure_dT_A,pilot_dpressure_dT_B")
```

并确保旧主列顺序和语义不变，仅追加 `pilot_` 前缀列。

- [ ] **Step 4: Run script dry run with local output path**

将输出路径临时改为仓库内可写目录后运行一次：  
Run: `julia --project=. scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl`  
Expected: 输出 detail/summary 文件且包含 `pilot_` 新列。

- [ ] **Step 5: Commit**

```bash
git add scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl
git commit -m "refactor(analysis): add pilot diff adapter metrics to t190 decomposition chain"
```

---

### Task 6: 完成 #81 评估文档与 Phase-2 决策建议

**Files:**
- Create: `docs/dev/active/issue-81-solver-diff-pilot-evaluation.md`

- [ ] **Step 1: Write evaluation doc template with concrete sections**

```markdown
# Issue #81 Solver/Diff Pilot Evaluation

## Scope
- Script chain: scripts/pnjl/calculate_derivatives.jl
- Analysis chain: scripts/analysis/relaxtime/t190_sigma_chain_decomposition.jl

## Numerical Consistency
- Baseline command:
- Pilot command:
- Key metric table (rtol/atol):

## Runtime Cost
- Baseline wall-time:
- Pilot wall-time:
- Relative overhead (%):

## Maintainability
- Files touched:
- Repeated derivative logic removed:
- New unified entrypoints:

## Phase-2 Decision
- Candidate callbacks to fold into ProblemSpec:
- Targets that should remain in diff layer:
- Migration template for next chains:
```

- [ ] **Step 2: Fill doc with measured values from actual runs**

运行后将实际数值写入，不保留空字段。

- [ ] **Step 3: Commit**

```bash
git add docs/dev/active/issue-81-solver-diff-pilot-evaluation.md
git commit -m "docs(dev): record issue 81 pilot evaluation and phase-2 decision items"
```

---

### Task 7: 运行 issue #81 验证命令与补充 smoke

**Files:**
- Modify (if needed): `docs/dev/active/issue-81-solver-diff-pilot-evaluation.md`

- [ ] **Step 1: Run required integration core**

Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'`  
Expected: PASS。

- [ ] **Step 2: Run required validation smoke**

Run: `julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'`  
Expected: PASS。

- [ ] **Step 3: Run supplemental unit smoke**

Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`  
Expected: PASS。

- [ ] **Step 4: Append exact pass/fail evidence into evaluation doc**

把命令、时间戳、结果摘要写入 `issue-81-solver-diff-pilot-evaluation.md`。

- [ ] **Step 5: Commit**

```bash
git add docs/dev/active/issue-81-solver-diff-pilot-evaluation.md
git commit -m "test(validation): capture issue 81 pilot verification evidence"
```

---

### Task 8: API 文档与 PR 收口

**Files:**
- Modify: `docs/api/models/solver/README.md`
- Modify: `docs/api/models/solver/Overview.md`
- Modify: `docs/api/models/solver/generated/Exports.md` (regenerated)

- [ ] **Step 1: Add pilot adapter API references in solver docs**

在 `README.md` / `Overview.md` 增补：

```markdown
- `build_pilot_diff_context(result; ...)`: #81 试点统一上下文构建入口。
- `eval_pilot_derivatives(ctx; target_names, param_names)`: #81 试点统一 Jacobian 求值入口。
```

- [ ] **Step 2: Regenerate exported API section (if repo script available)**

Run: `julia --project=. scripts/dev/check_docs_consistency.jl`  
Expected: PASS；若检查提示导出清单需要更新，按仓内脚本更新 `generated/Exports.md`。

- [ ] **Step 3: Final local sanity before PR**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/solver/test_solver_diff_pilot_adapters.jl"; include("tests/unit/runtests.jl")'`  
Expected: PASS。

- [ ] **Step 4: Commit**

```bash
git add docs/api/models/solver/README.md docs/api/models/solver/Overview.md docs/api/models/solver/generated/Exports.md
git commit -m "docs(api): document issue 81 pilot diff adapter entrypoints"
```

- [ ] **Step 5: PR closure workflow**

```bash
git push -u origin feat/solver-diff-pilot-81
gh pr create --repo w5851/Julia_RelaxTime --title "refactor(models/solver): pilot unified diff adapters in script and analysis chains (#81)" --body "$(cat <<'EOF'
## Summary
- add unified pilot diff adapter service under solver/diff
- migrate one script chain and one analysis chain to shared Jacobian path
- provide #81 evaluation record with phase-2 actionable decisions

## Verification
- julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="core"; include("tests/integration/runtests.jl")'
- julia --project=. -e 'ENV["VALIDATION_PROFILE"]="smoke"; include("tests/validation/runtests.jl")'
- julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
EOF
)"
```

评审后按既定流程：修复评论 -> 重跑验证 -> squash merge -> 删除远端分支 -> 本地切回主干并删除特性分支。

---

## Plan Self-Review

- Spec coverage: 已覆盖两条试点链路接入、统一服务抽象、测试验证、性能稳定性评估、Phase-2 决策输出。
- Placeholder scan: 无 `TBD/TODO/implement later` 占位语句。
- Type consistency: 计划中统一使用 `build_pilot_diff_context` / `eval_pilot_derivatives`、`target_names::Vector{Symbol}`、`param_names::Vector{Symbol}`。
