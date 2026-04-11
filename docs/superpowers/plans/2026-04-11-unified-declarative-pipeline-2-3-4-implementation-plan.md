# Unified Declarative Pipeline (2+3+4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 `Models` workflow/scan 与 relaxtime orchestrator 建立统一声明式编排骨架，在不破坏现有公开入口的前提下输出统一 manifest 与回归锚点。

**Architecture:** 保留 `src/models/workflow` 现有执行内核（types/runner/hash/manifest），新增 adapter + catalog + io 扩展层，将 2/3/4 三类链路映射到统一 stage 语义。脚本入口仅做 CLI 薄适配，业务编排集中在 `Models` 侧。迁移采用四波次并保留回退开关。

**Tech Stack:** Julia 1.10+, `src/models/workflow/*`, `src/models/entrypoints.jl`, `src/models/Models.jl`, `scripts/models/run_unified_scan.jl`, `scripts/relaxtime/run_relaxtime_orchestrator.jl`, tests under `tests/unit|integration|regression`.

---

## File Structure Map

- Create: `src/models/workflow/adapters/CommonAdapterUtils.jl` (输入规范化、alias 映射、错误分类)
- Create: `src/models/workflow/adapters/WorkflowAdapter.jl` (workflow 家族适配)
- Create: `src/models/workflow/adapters/ScanAdapter.jl` (scan 家族适配)
- Create: `src/models/workflow/adapters/RelaxtimeOrchestratorAdapter.jl` (relaxtime orchestrator 家族适配)
- Create: `src/models/workflow/catalog/WorkflowCatalog.jl` (workflow stage 声明)
- Create: `src/models/workflow/catalog/ScanCatalog.jl` (scan stage 声明)
- Create: `src/models/workflow/catalog/RelaxtimeOrchestratorCatalog.jl` (orchestrator stage 声明)
- Create: `src/models/workflow/io/ManifestExtensions.jl` (pipeline_family/baseline_suite/physics_profile/adapter_version)
- Modify: `src/models/workflow/PipelineRunner.jl` (manifest 扩展字段写出)
- Modify: `src/models/entrypoints.jl` (入口接线到 adapter)
- Modify: `src/models/Models.jl` (include/export 新模块)
- Modify: `scripts/models/run_unified_scan.jl` (薄入口化)
- Modify: `scripts/relaxtime/run_relaxtime_orchestrator.jl` (薄入口化)
- Create: `tests/unit/models/workflow/test_adapter_normalization.jl`
- Create: `tests/unit/models/workflow/test_catalog_contracts_234.jl`
- Create: `tests/unit/models/workflow/test_manifest_extensions_234.jl`
- Create: `tests/integration/models/test_scan_pipeline_runner_smoke.jl`
- Create: `tests/integration/relaxtime/test_orchestrator_pipeline_runner_smoke.jl`
- Create: `tests/regression/models/test_workflow_scan_orchestrator_consistency_234.jl`

### Task 1: 建立 Adapter 公共层与错误分类

**Files:**
- Create: `src/models/workflow/adapters/CommonAdapterUtils.jl`
- Test: `tests/unit/models/workflow/test_adapter_normalization.jl`

- [ ] **Step 1: 写失败测试（归一化与错误分类）**

```julia
@testset "adapter utils normalization" begin
    normalized = Models.normalize_adapter_kwargs((T_mev=150.0, xi=0.1), Dict("T_MeV" => :T_mev))
    @test normalized.T_mev == 150.0
    @test normalized.xi == 0.1

    err = Models.classify_pipeline_error(ArgumentError("bad input"))
    @test err == :input_validation_error
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_adapter_normalization.jl")'`
Expected: FAIL，提示 `normalize_adapter_kwargs` 或 `classify_pipeline_error` 未定义。

- [ ] **Step 3: 实现最小 Adapter 公共工具**

```julia
module CommonAdapterUtils

export normalize_adapter_kwargs, classify_pipeline_error

function normalize_adapter_kwargs(kwargs::NamedTuple, alias_map::Dict{String,Symbol})
    out = Dict{Symbol,Any}(pairs(kwargs))
    for (k, sym) in alias_map
        if haskey(out, Symbol(k)) && !haskey(out, sym)
            out[sym] = out[Symbol(k)]
        end
    end
    return (; (k => out[k] for k in sort!(collect(keys(out)); by=String))...)
end

function classify_pipeline_error(err)::Symbol
    err isa ArgumentError && return :input_validation_error
    err isa IOError && return :artifact_io_error
    msg = lowercase(sprint(showerror, err))
    occursin("converge", msg) && return :numerical_convergence_error
    return :unexpected_error
end

end
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_adapter_normalization.jl")'`
Expected: PASS，测试集全部通过。

- [ ] **Step 5: 提交**

```bash
git add src/models/workflow/adapters/CommonAdapterUtils.jl tests/unit/models/workflow/test_adapter_normalization.jl
git commit -m "refactor(models/workflow): add adapter normalization and error taxonomy helpers"
```

### Task 2: 建立三类 Catalog（workflow/scan/orchestrator）

**Files:**
- Create: `src/models/workflow/catalog/WorkflowCatalog.jl`
- Create: `src/models/workflow/catalog/ScanCatalog.jl`
- Create: `src/models/workflow/catalog/RelaxtimeOrchestratorCatalog.jl`
- Test: `tests/unit/models/workflow/test_catalog_contracts_234.jl`

- [ ] **Step 1: 写失败测试（stage 合同完整性）**

```julia
@testset "catalog contracts 234" begin
    wf = Models.workflow_pipeline_stage_ids()
    sc = Models.scan_pipeline_stage_ids()
    ro = Models.relaxtime_orchestrator_stage_ids()
    expected = (:prepare_inputs, :solve_core, :postprocess, :export_artifacts, :emit_repro_manifest)
    @test Tuple(wf) == expected
    @test Tuple(sc) == expected
    @test Tuple(ro) == expected
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_catalog_contracts_234.jl")'`
Expected: FAIL，catalog API 尚未导出。

- [ ] **Step 3: 实现最小 Catalog 声明**

```julia
const PIPELINE_STAGE_SKELETON = (
    :prepare_inputs,
    :solve_core,
    :postprocess,
    :export_artifacts,
    :emit_repro_manifest,
)

workflow_pipeline_stage_ids() = PIPELINE_STAGE_SKELETON
scan_pipeline_stage_ids() = PIPELINE_STAGE_SKELETON
relaxtime_orchestrator_stage_ids() = PIPELINE_STAGE_SKELETON
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_catalog_contracts_234.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add src/models/workflow/catalog tests/unit/models/workflow/test_catalog_contracts_234.jl
git commit -m "refactor(models/workflow): add unified stage catalogs for workflow scan orchestrator"
```

### Task 3: Manifest 扩展字段接入 Runner

**Files:**
- Create: `src/models/workflow/io/ManifestExtensions.jl`
- Modify: `src/models/workflow/PipelineRunner.jl`
- Test: `tests/unit/models/workflow/test_manifest_extensions_234.jl`

- [ ] **Step 1: 写失败测试（manifest 扩展字段存在）**

```julia
@testset "manifest extensions" begin
    meta = Models.build_manifest_extensions((; pipeline_family=:workflow, baseline_suite=:smoke, physics_profile="default", adapter_version="v1"))
    @test meta["pipeline_family"] == "workflow"
    @test meta["baseline_suite"] == "smoke"
    @test meta["physics_profile"] == "default"
    @test meta["adapter_version"] == "v1"
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_manifest_extensions_234.jl")'`
Expected: FAIL，`build_manifest_extensions` 未定义。

- [ ] **Step 3: 实现 io 扩展与 runner 接线**

```julia
function build_manifest_extensions(meta::NamedTuple)
    return Dict(
        "pipeline_family" => String(get(meta, :pipeline_family, :unknown)),
        "baseline_suite" => String(get(meta, :baseline_suite, :none)),
        "physics_profile" => String(get(meta, :physics_profile, "default")),
        "adapter_version" => String(get(meta, :adapter_version, "v1")),
    )
end

# In _write_manifest(...) of PipelineRunner.jl
extra = haskey(ctx.state, :manifest_extensions) ? ctx.state[:manifest_extensions] : Dict{String,Any}()
merge!(pipeline_meta, extra)
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/unit/models/workflow/test_manifest_extensions_234.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add src/models/workflow/io/ManifestExtensions.jl src/models/workflow/PipelineRunner.jl tests/unit/models/workflow/test_manifest_extensions_234.jl
git commit -m "refactor(models/workflow): extend manifest metadata for unified pipeline families"
```

### Task 4: 接入 WorkflowAdapter 到 Models 入口

**Files:**
- Create: `src/models/workflow/adapters/WorkflowAdapter.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/integration/models/test_scan_pipeline_runner_smoke.jl`

- [ ] **Step 1: 写失败集成测试（workflow 路由可产出统一 manifest）**

```julia
@testset "workflow adapter smoke" begin
    tmp = mktempdir()
    res = Models.run_workflow_pipeline(:transport; T_fm=0.15, mu_fm=0.0, xi=0.0, output_dir=tmp)
    @test haskey(res, :run_manifest)
    @test isfile(res.run_manifest)
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'include("tests/integration/models/test_scan_pipeline_runner_smoke.jl")'`
Expected: FAIL，`run_workflow_pipeline` 未实现。

- [ ] **Step 3: 实现 workflow adapter 最小接线**

```julia
function run_workflow_pipeline(kind::Symbol; kwargs...)
    spec, stages, ctx = build_workflow_pipeline(kind; kwargs...)
    result = run_pipeline(spec, stages, ctx; manifest_path=joinpath(kwargs.output_dir, "run_manifest.json"))
    result.success || error("workflow pipeline failed at $(result.failed_stage)")
    return (run_manifest=result.manifest_path,)
end
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/integration/models/test_scan_pipeline_runner_smoke.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add src/models/workflow/adapters/WorkflowAdapter.jl src/models/entrypoints.jl src/models/Models.jl tests/integration/models/test_scan_pipeline_runner_smoke.jl
git commit -m "refactor(models): route workflow entrypoints through declarative adapter runner"
```

### Task 5: 接入 ScanAdapter 与 unified scan 脚本薄入口化

**Files:**
- Create: `src/models/workflow/adapters/ScanAdapter.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `scripts/models/run_unified_scan.jl`
- Test: `tests/integration/models/test_scan_pipeline_runner_smoke.jl`

- [ ] **Step 1: 写失败测试（scan 路由与 manifest 字段）**

```julia
@testset "scan adapter smoke" begin
    tmp = mktempdir()
    stats = Models.run_scan_pipeline(:tmu; model_kind=:PNJL, T_values=[150.0], mu_values=[0.0], xi_values=[0.0], output_path=joinpath(tmp, "out.csv"))
    @test hasproperty(stats, :manifest_path)
    @test isfile(stats.manifest_path)
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'include("tests/integration/models/test_scan_pipeline_runner_smoke.jl")'`
Expected: FAIL，`run_scan_pipeline` 未定义。

- [ ] **Step 3: 实现 ScanAdapter 与脚本路由替换**

```julia
# entrypoints.jl
function run_scan_pipeline(kind::Symbol; kwargs...)
    return ScanAdapter.run_scan_pipeline(kind; kwargs...)
end

# run_unified_scan.jl
if args[1] == "scan"
    return Main.Models.run_scan_pipeline(Symbol(args[2]); kwargs...)
end
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/integration/models/test_scan_pipeline_runner_smoke.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add src/models/workflow/adapters/ScanAdapter.jl src/models/entrypoints.jl scripts/models/run_unified_scan.jl tests/integration/models/test_scan_pipeline_runner_smoke.jl
git commit -m "refactor(models/scans): unify scan orchestration via declarative adapter"
```

### Task 6: 接入 RelaxtimeOrchestratorAdapter 与脚本薄入口化

**Files:**
- Create: `src/models/workflow/adapters/RelaxtimeOrchestratorAdapter.jl`
- Modify: `src/models/entrypoints.jl`
- Modify: `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- Test: `tests/integration/relaxtime/test_orchestrator_pipeline_runner_smoke.jl`

- [ ] **Step 1: 写失败测试（orchestrator 声明式路由）**

```julia
@testset "relaxtime orchestrator adapter smoke" begin
    tmp = mktempdir()
    out = Models.run_relaxtime_orchestrator_pipeline(:cross_section; outdir=tmp)
    @test haskey(out, :manifest_path)
    @test isfile(out.manifest_path)
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_orchestrator_pipeline_runner_smoke.jl")'`
Expected: FAIL，orchestrator adapter API 未实现。

- [ ] **Step 3: 实现 adapter 并替换脚本主逻辑调用**

```julia
# entrypoints.jl
function run_relaxtime_orchestrator_pipeline(cmd::Symbol; kwargs...)
    return RelaxtimeOrchestratorAdapter.run_orchestrator_pipeline(cmd; kwargs...)
end

# run_relaxtime_orchestrator.jl main path
out = Main.Models.run_relaxtime_orchestrator_pipeline(Symbol(cmd); opts...)
println("[orchestrator] manifest=$(out.manifest_path)")
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_orchestrator_pipeline_runner_smoke.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add src/models/workflow/adapters/RelaxtimeOrchestratorAdapter.jl src/models/entrypoints.jl scripts/relaxtime/run_relaxtime_orchestrator.jl tests/integration/relaxtime/test_orchestrator_pipeline_runner_smoke.jl
git commit -m "refactor(relaxtime): route orchestrator through models declarative adapter"
```

### Task 7: 回归一致性与闸门收口

**Files:**
- Create: `tests/regression/models/test_workflow_scan_orchestrator_consistency_234.jl`
- Modify: `docs/api/models/workflows/Overview.md`
- Modify: `docs/api/relaxtime/workflow/TransportWorkflow.md`

- [ ] **Step 1: 写失败回归测试（旧入口 vs 新入口一致性）**

```julia
@testset "pipeline 234 consistency" begin
    legacy = Models.solve_gap_and_transport(0.15, 0.0; xi=0.0, compute_tau=false, tau=(u=1.0,d=1.0,s=1.0,ubar=1.0,dbar=1.0,sbar=1.0), compute_bulk=false)
    via_pipeline = Models.run_workflow_pipeline(:transport; T_fm=0.15, mu_fm=0.0, xi=0.0)
    @test isapprox(legacy.transport.eta, via_pipeline.transport.eta; rtol=1e-6, atol=1e-8)
end
```

- [ ] **Step 2: 运行回归确认失败**

Run: `julia --project=. -e 'include("tests/regression/models/test_workflow_scan_orchestrator_consistency_234.jl")'`
Expected: FAIL，直到新路径输出合同与旧路径对齐。

- [ ] **Step 3: 补齐对齐逻辑与文档更新**

```julia
# Ensure pipeline output carries legacy-comparable transport fields
produced = Dict{Symbol,Any}(
    :transport => direct.transport,
    :equilibrium => direct.equilibrium,
)
```

- [ ] **Step 4: 运行回归与最小闸门命令确认通过**

Run: `julia --project=. -e 'include("tests/regression/models/test_workflow_scan_orchestrator_consistency_234.jl")'`
Expected: PASS。

Run: `julia --project=. scripts/dev/check_models_entry_contract.jl`
Expected: PASS，未破坏 Models 公开入口合同。

- [ ] **Step 5: 提交**

```bash
git add tests/regression/models/test_workflow_scan_orchestrator_consistency_234.jl docs/api/models/workflows/Overview.md docs/api/relaxtime/workflow/TransportWorkflow.md
git commit -m "refactor(models/workflow): close 2+3+4 declarative migration with regression parity"
```

## Final Verification Batch

- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_models_entry_contract.jl`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`

## Commit Style Note

Use existing repository style sampled from recent history:

- `refactor(models/workflow): ...`
- `refactor(models): ...`
- `refactor(relaxtime): ...`
- `docs(dev): ...` (if docs-only)

Do not invent new prefix format.
