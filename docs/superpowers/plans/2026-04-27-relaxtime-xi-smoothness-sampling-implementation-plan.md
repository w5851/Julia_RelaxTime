# Relaxtime Xi Smoothness Sampling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在固定参数预算下可复现地产生 24 组 `(T, μ)` 样本，批量生成 `xi` 输运曲线，并输出“自动首筛 + 人工复核”平滑性结论。

**Architecture:** 采用三段式流水线：`采样清单生成 -> 批量计算/绘图 -> 平滑性评分与复核合并`。每个阶段产物均落盘为带元数据的 CSV/图像 sidecar，并由 run manifest 串联。实现以脚本层为主，不改动核心物理求解内核。

**Tech Stack:** Julia 1.10+ (`scripts/relaxtime/*.jl`), Python plotting (`scripts/plot_scan_csv.py` + `scripts/common/provenance_image.py`), Julia `Test` (unit + integration smoke), existing `run_gap_transport_scan.jl` workflow.

---

## Preconditions

- 在独立 worktree 中执行该计划。
- 使用现有参考数据：
  - `data/reference/pnjl/boundary.csv`
  - `data/reference/pnjl/crossover.csv`
  - `data/reference/pnjl/cep.csv`
- 本计划不重生成相图数据，不修改 `src/models` 与 `src/relaxtime` 物理核心行为。

## File Map

- Create: `scripts/relaxtime/xi_smoothness_sampling_lib.jl`
  - 抽样核心逻辑（随机层 + 相变线邻域层）与清单写出。
- Create: `scripts/relaxtime/generate_xi_smoothness_params.jl`
  - 抽样 CLI，输出 `params_24_seed<seed>.csv`。
- Create: `scripts/relaxtime/run_xi_smoothness_batch.jl`
  - 读取参数清单，逐样本调用 `run_gap_transport_scan.jl`，落 run manifest。
- Create: `scripts/relaxtime/plot_xi_smoothness_batch.jl`
  - 对每个样本调用 `scripts/plot_scan_csv.py` 产图并触发图像 provenance sidecar。
- Create: `scripts/relaxtime/evaluate_xi_smoothness.jl`
  - 计算 `S2 / S1jump / N_spike`，输出 `scores/flags/review_queue`。
- Create: `scripts/relaxtime/merge_xi_smoothness_review.jl`
  - 合并人工复核标签，生成 `smoothness_final.csv`。
- Create: `tests/unit/relaxtime/test_xi_smoothness_sampling_lib.jl`
  - 抽样域约束、可复现性、字段契约单元测试。
- Create: `tests/unit/relaxtime/test_xi_smoothness_metrics.jl`
  - 指标计算与双门槛分类单元测试。
- Create: `tests/integration/relaxtime/test_xi_smoothness_batch_runner_smoke.jl`
  - 批量 runner `--dry-run` + manifest 契约。
- Create: `tests/integration/relaxtime/test_xi_smoothness_plot_smoke.jl`
  - 图像输出与 `.provenance.json` sidecar 契约。
- Create: `tests/integration/relaxtime/test_xi_smoothness_evaluation_smoke.jl`
  - `scores/flags/review_queue/final` 文件契约。
- Create: `tests/fixtures/relaxtime/xi_smoothness/*.csv`
  - 轻量 fixture（避免重算物理流程）。
- Modify: `tests/unit/runtests.jl`
  - 纳入新增 unit tests。
- Modify: `tests/integration/runtests.jl`
  - 纳入新增 integration smoke tests。

## Task 1: Build Sampling Library + CLI (TDD)

**Files:**
- Create: `scripts/relaxtime/xi_smoothness_sampling_lib.jl`
- Create: `scripts/relaxtime/generate_xi_smoothness_params.jl`
- Create: `tests/unit/relaxtime/test_xi_smoothness_sampling_lib.jl`
- Modify: `tests/unit/runtests.jl`

- [ ] **Step 1: Write failing unit tests for sampling contract**

```julia
@testset "xi smoothness sampling contract" begin
    rows = XiSmoothnessSampling.sample_params(
        24;
        seed=20260427,
        random_count=12,
        near_count=12,
        T_range=(50.0, 270.0),
        muq_range=(0.0, 360.0),
    )
    @test length(rows) == 24
    @test count(r -> r.source == "random_uniform", rows) == 12
    @test count(r -> r.source == "near_phase_line", rows) == 12
    @test all(50.0 <= r.T_MeV <= 270.0 for r in rows)
    @test all(0.0 <= r.muq_MeV <= 360.0 for r in rows)
end
```

- [ ] **Step 2: Run test and confirm failure**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_xi_smoothness_sampling_lib.jl")'`

Expected: FAIL with missing module/file error.

- [ ] **Step 3: Implement minimal sampling library**

```julia
module XiSmoothnessSampling

using Random

Base.@kwdef struct SampleRow
    sample_id::String
    source::String
    T_MeV::Float64
    muq_MeV::Float64
    muB_MeV::Float64
    anchor_type::String
    anchor_T_MeV::Float64
    anchor_muq_MeV::Float64
    delta_T::Float64
    delta_muq::Float64
    rng_seed::Int
end

function sample_params(total::Int; seed::Int, random_count::Int, near_count::Int, T_range, muq_range)
    rng = MersenneTwister(seed)
    # full implementation in this file
end

end
```

- [ ] **Step 4: Implement CLI writer with metadata header**

```julia
#!/usr/bin/env julia
include(joinpath(@__DIR__, "xi_smoothness_sampling_lib.jl"))
using .XiSmoothnessSampling

# parse args, call sample_params(...), write:
# data/outputs/results/relaxtime/xi_smoothness_sampling/sampling/params_24_seed<seed>.csv
```

- [ ] **Step 5: Re-run unit test and confirm pass**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_xi_smoothness_sampling_lib.jl"; include("tests/unit/runtests.jl")'`

Expected: PASS。

- [ ] **Step 6: Commit sampling layer**

```bash
git add scripts/relaxtime/xi_smoothness_sampling_lib.jl scripts/relaxtime/generate_xi_smoothness_params.jl tests/unit/relaxtime/test_xi_smoothness_sampling_lib.jl tests/unit/runtests.jl
git commit -m "feat(relaxtime): add reproducible xi smoothness parameter sampling"
```

## Task 2: Add Batch Runner With Resume/Manifest (TDD)

**Files:**
- Create: `scripts/relaxtime/run_xi_smoothness_batch.jl`
- Create: `tests/integration/relaxtime/test_xi_smoothness_batch_runner_smoke.jl`
- Modify: `tests/integration/runtests.jl`

- [ ] **Step 1: Write failing integration smoke for dry-run manifest**

```julia
@testset "xi smoothness batch dry-run" begin
    outdir = mktempdir()
    params_csv = joinpath(outdir, "params.csv")
    # write tiny 2-row params fixture
    run(`julia --project=. scripts/relaxtime/run_xi_smoothness_batch.jl --params $params_csv --out-root $outdir --dry-run`)
    @test isfile(joinpath(outdir, "run_manifest.json"))
end
```

- [ ] **Step 2: Run test and confirm failure**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_xi_smoothness_batch_runner_smoke.jl")'`

Expected: FAIL with script not found.

- [ ] **Step 3: Implement batch runner with per-sample command assembly**

```julia
cmd = `julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin $T --tmax $T --tstep 1 --mubmin $muB --mubmax $muB --mubstep 1 --xi-list $xi_list --mode finite_15 --tau-interpolation-mode hybrid_threshold --output $out_csv --failed-points-output $failed_csv --resume`
```

- [ ] **Step 4: Add `--dry-run` and `--resume` behavior**

```julia
if opts.dry_run
    status = "skipped"
else
    proc = run(ignorestatus(cmd))
    status = proc.exitcode == 0 ? "success" : "failed"
end
```

- [ ] **Step 5: Write run manifest JSON with required fields**

```julia
manifest_row = Dict(
    "sample_id" => sample_id,
    "status" => status,
    "elapsed_seconds" => elapsed,
    "input_sampling_csv" => opts.params_csv,
    "result_csv" => out_csv,
    "failed_points_path" => failed_csv,
)
```

- [ ] **Step 6: Re-run integration smoke and confirm pass**

Run: `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_xi_smoothness_batch_runner_smoke.jl"; include("tests/integration/runtests.jl")'`

Expected: PASS。

- [ ] **Step 7: Commit batch runner**

```bash
git add scripts/relaxtime/run_xi_smoothness_batch.jl tests/integration/relaxtime/test_xi_smoothness_batch_runner_smoke.jl tests/integration/runtests.jl
git commit -m "feat(relaxtime): add xi smoothness batch runner with manifest"
```

## Task 3: Add Plot Batch Script With Provenance Sidecars (TDD)

**Files:**
- Create: `scripts/relaxtime/plot_xi_smoothness_batch.jl`
- Create: `tests/integration/relaxtime/test_xi_smoothness_plot_smoke.jl`
- Create: `tests/fixtures/relaxtime/xi_smoothness/sample_transport_vs_xi.csv`

- [ ] **Step 1: Write failing integration smoke for png + sidecar generation**

```julia
@test isfile(joinpath(fig_dir, "eta_over_s_vs_xi.png"))
@test isfile(joinpath(fig_dir, "eta_over_s_vs_xi.png.provenance.json"))
```

- [ ] **Step 2: Run test and confirm failure**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_xi_smoothness_plot_smoke.jl")'`

Expected: FAIL with missing plotting batch script.

- [ ] **Step 3: Implement plotting batch script by invoking existing plot tool**

```julia
plot_cmd = `python scripts/plot_scan_csv.py --mode lines --csv $csv_path --x xi --ys tau_u,tau_s,eta_over_s,zeta_over_s,sigma_over_T --out-dir $fig_dir --check`
run(plot_cmd)
```

- [ ] **Step 4: Ensure command string is preserved for provenance**

```julia
println("[plot] command=", plot_cmd)
```

- [ ] **Step 5: Re-run integration smoke and confirm pass**

Run: `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_xi_smoothness_plot_smoke.jl"; include("tests/integration/runtests.jl")'`

Expected: PASS，并且 sidecar JSON 含 `script_path=\"scripts/plot_scan_csv.py\"`。

- [ ] **Step 6: Commit plot batch layer**

```bash
git add scripts/relaxtime/plot_xi_smoothness_batch.jl tests/integration/relaxtime/test_xi_smoothness_plot_smoke.jl tests/fixtures/relaxtime/xi_smoothness/sample_transport_vs_xi.csv
git commit -m "feat(relaxtime): add xi smoothness plotting batch with provenance"
```

## Task 4: Implement Smoothness Metrics + Flags + Review Queue (TDD)

**Files:**
- Create: `scripts/relaxtime/evaluate_xi_smoothness.jl`
- Create: `tests/unit/relaxtime/test_xi_smoothness_metrics.jl`
- Create: `tests/integration/relaxtime/test_xi_smoothness_evaluation_smoke.jl`
- Create: `tests/fixtures/relaxtime/xi_smoothness/eval_input/*.csv`

- [ ] **Step 1: Write failing unit tests for `S2 / S1jump / N_spike`**

```julia
@test XiSmoothnessEval.metric_s2([1.0, 1.0, 1.0]) == 0.0
@test XiSmoothnessEval.metric_s2([1.0, 2.0, 1.0]) > 0.0
@test XiSmoothnessEval.metric_spikes([0.0, 0.2, 0.0], 0.1) == 1
```

- [ ] **Step 2: Write failing unit tests for dual-threshold classification**

```julia
@test XiSmoothnessEval.classify(0.01, 0.01, 0; loose=0.05, strict=0.2) == "smooth"
@test XiSmoothnessEval.classify(0.30, 0.25, 2; loose=0.05, strict=0.2) == "not_smooth"
@test XiSmoothnessEval.classify(0.10, 0.08, 1; loose=0.05, strict=0.2) == "suspect"
```

- [ ] **Step 3: Run unit tests and confirm failure**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_xi_smoothness_metrics.jl")'`

Expected: FAIL with missing evaluator module.

- [ ] **Step 4: Implement evaluator script and CSV outputs**

```julia
# write:
# smoothness_scores.csv
# smoothness_flags.csv
# review_queue.csv
```

- [ ] **Step 5: Add integration smoke for output contract**

```julia
@test isfile(joinpath(outdir, "smoothness_scores.csv"))
@test isfile(joinpath(outdir, "smoothness_flags.csv"))
@test isfile(joinpath(outdir, "review_queue.csv"))
```

- [ ] **Step 6: Re-run unit + integration tests and confirm pass**

Run:

- `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_xi_smoothness_metrics.jl"; include("tests/unit/runtests.jl")'`
- `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_xi_smoothness_evaluation_smoke.jl"; include("tests/integration/runtests.jl")'`

Expected: PASS。

- [ ] **Step 7: Commit evaluator layer**

```bash
git add scripts/relaxtime/evaluate_xi_smoothness.jl tests/unit/relaxtime/test_xi_smoothness_metrics.jl tests/integration/relaxtime/test_xi_smoothness_evaluation_smoke.jl tests/fixtures/relaxtime/xi_smoothness/eval_input
git commit -m "feat(relaxtime): add xi smoothness scoring and review queue"
```

## Task 5: Implement Manual Review Merge Script (TDD)

**Files:**
- Create: `scripts/relaxtime/merge_xi_smoothness_review.jl`
- Create: `tests/unit/relaxtime/test_xi_smoothness_review_merge.jl`
- Modify: `tests/unit/runtests.jl`

- [ ] **Step 1: Write failing unit test for merge contract**

```julia
@testset "review merge" begin
    # flags.csv + manual_review.csv -> smoothness_final.csv
    @test isfile(final_csv)
    @test occursin("confirm_smooth", read(final_csv, String))
end
```

- [ ] **Step 2: Run test and confirm failure**

Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_xi_smoothness_review_merge.jl")'`

Expected: FAIL with missing script.

- [ ] **Step 3: Implement merge script**

```julia
# manual labels: confirm_smooth / confirm_not_smooth / undecided
# output: smoothness_final.csv
```

- [ ] **Step 4: Re-run unit test and confirm pass**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_xi_smoothness_review_merge.jl"; include("tests/unit/runtests.jl")'`

Expected: PASS。

- [ ] **Step 5: Commit merge step**

```bash
git add scripts/relaxtime/merge_xi_smoothness_review.jl tests/unit/relaxtime/test_xi_smoothness_review_merge.jl tests/unit/runtests.jl
git commit -m "feat(relaxtime): add manual review merge for xi smoothness labels"
```

## Task 6: Pipeline Documentation + Operational Commands

**Files:**
- Modify: `docs/superpowers/specs/2026-04-27-relaxtime-xi-smoothness-sampling-design.md`
- Create: `docs/dev/active/2026-04-27-xi-smoothness-runbook.md`

- [ ] **Step 1: Add executable runbook commands (no placeholders)**

```bash
julia --project=. scripts/relaxtime/generate_xi_smoothness_params.jl --seed 20260427 --out data/outputs/results/relaxtime/xi_smoothness_sampling/sampling/params_24_seed20260427.csv
julia --project=. scripts/relaxtime/run_xi_smoothness_batch.jl --params data/outputs/results/relaxtime/xi_smoothness_sampling/sampling/params_24_seed20260427.csv --xi-min -0.5 --xi-max 0.5 --xi-step 0.02
julia --project=. scripts/relaxtime/plot_xi_smoothness_batch.jl --manifest data/outputs/results/relaxtime/xi_smoothness_sampling/sampling/run_manifest.json
julia --project=. scripts/relaxtime/evaluate_xi_smoothness.jl --manifest data/outputs/results/relaxtime/xi_smoothness_sampling/sampling/run_manifest.json
```

- [ ] **Step 2: Add manual review template example**

```csv
sample_id,field,manual_label,reviewer,note
S013,tau_s,confirm_not_smooth,wmzx,"kink near xi=-0.34"
```

- [ ] **Step 3: Commit docs/runbook**

```bash
git add docs/dev/active/2026-04-27-xi-smoothness-runbook.md docs/superpowers/specs/2026-04-27-relaxtime-xi-smoothness-sampling-design.md
git commit -m "docs(relaxtime): add xi smoothness sampling runbook"
```

## Task 7: Final Verification Gate

**Files:**
- Modify: none

- [ ] **Step 1: Run targeted unit profile**

Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; ENV["UNIT_FILES"]="relaxtime/test_xi_smoothness_sampling_lib.jl,relaxtime/test_xi_smoothness_metrics.jl,relaxtime/test_xi_smoothness_review_merge.jl"; include("tests/unit/runtests.jl")'`

Expected: PASS。

- [ ] **Step 2: Run targeted integration profile**

Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; ENV["INTEGRATION_FILES"]="relaxtime/test_xi_smoothness_batch_runner_smoke.jl,relaxtime/test_xi_smoothness_plot_smoke.jl,relaxtime/test_xi_smoothness_evaluation_smoke.jl"; include("tests/integration/runtests.jl")'`

Expected: PASS。

- [ ] **Step 3: Run governance checks for scripts/docs consistency**

Run:

- `julia --project=. scripts/dev/check_script_entrypoints.jl`
- `julia --project=. scripts/dev/check_docs_consistency.jl`
- `julia --project=. scripts/dev/check_relaxtime_script_governance.jl`

Expected: PASS。

- [ ] **Step 4: Final commit if verification passed**

```bash
git add .
git commit -m "feat(relaxtime): deliver xi smoothness sampling and evaluation pipeline"
```

## Spec Coverage Self-Check

- 抽样 24 组（12+12）: Task 1。
- 参数域与固定 seed 可复现: Task 1。
- 每样本独立 `run_gap_transport_scan` + `xi=-0.5:0.02:0.5`: Task 2。
- 图像输出与元数据 sidecar: Task 3。
- 自动判定 `S2/S1jump/N_spike` + 双门槛: Task 4。
- 人工复核合并 `smoothness_final.csv`: Task 5。
- 运行治理（failed sidecar/resume/manifest）: Task 2。

无缺失项。

## Placeholder Scan

- 无 `TODO/TBD/implement later`。
- 每个代码变更步骤均给出文件与命令示例。
- 每个验证步骤均给出明确命令与预期结果。

## Type/Name Consistency Check

- 统一命名：`sample_id/source/T_MeV/muq_MeV/muB_MeV`。
- 统一输出名：`smoothness_scores.csv / smoothness_flags.csv / review_queue.csv / smoothness_final.csv`。
- 统一标签：`smooth/suspect/not_smooth` 与 `confirm_smooth/confirm_not_smooth/undecided`。
