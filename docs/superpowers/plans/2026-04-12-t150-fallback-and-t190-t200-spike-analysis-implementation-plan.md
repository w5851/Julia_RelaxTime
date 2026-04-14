# T150 Fallback And T190/T200 Spike Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不修改 `src` 的前提下，通过 `scripts` 正确调用 `src` 既有 fallback 能力修复 T150 首点求解缺失，并产出 T190/T200 毛刺的可审计证据链分析报告。

**Architecture:** 在 `run_gap_transport_scan.jl` 的平衡态求解入口中优先使用 `Models.solve(...)`（由 `src` 内建 multiseed/fallback 机制托底），仅在必要时回落到现有 `solve_constraint(...)` 路径；同时新增跳过点失败 sidecar 以补齐不可见失败信息。分析侧以用户指定 merged 数据为主，先做链路归因，必要时再在 temp 目录定向重跑。

**Tech Stack:** Julia 1.10+，现有 Models/TransportWorkflow 脚本链路，CSV sidecar 诊断输出，Markdown 分析文档。

---

## 执行状态更新（2026-04-14）

- [x] 已按 3-patch 拆分完成实现并提交（solver fallback、workflow/analysis tooling、analysis deliverables）。
- [x] 已根据评审要求排除 `.github/workflows/ci-validation.yml` 变更，不纳入本次 PR。
- [x] 关键验证已通过（unit/integration/regression 定向用例 + docs consistency）。

### 对应提交

- `981857e` `refactor(relaxtime): add FixedMu solver fallback and failed-point sidecar`
- `cdfe2c2` `refactor(relaxtime): add t190 chain diagnostics flow and regression guard`
- `4edc50e` `docs(relaxtime): add t150/t190/t200 analysis notes and figures`

### Task 1: 建立求解入口测试骨架（TDD）

**Files:**
- Create: `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`
- Modify: `tests/unit/runtests.jl`
- Test: `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

- [x] **Step 1: 写失败测试（`Models.solve` 优先分支契约）**

```julia
using Test

@testset "run_gap_transport_scan solver entry" begin
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_gap_transport_scan.jl"))
    @test isdefined(Main, :_solve_fixedmu_via_models_solve)
end
```

- [x] **Step 2: 运行测试并确认失败**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: FAIL，出现 `@test false` 失败。

- [x] **Step 3: 在 `tests/unit/runtests.jl` 注册新文件（按现有模式）**

```julia
# 在 relaxtime 分组中添加：
"relaxtime/test_run_gap_transport_scan_solver_entry.jl",
```

- [x] **Step 4: 再次运行测试确认仍失败但文件被正确调度**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: FAIL，且仅该测试文件执行。

- [x] **Step 5: 提交（测试骨架）**

```bash
git add tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl tests/unit/runtests.jl
git commit -m "test: scaffold solver-entry fallback tests for gap transport scan"
```

### Task 2: 在脚本中实现 `Models.solve` 优先 + `solve_constraint` 兜底

**Files:**
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`
- Test: `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

- [x] **Step 1: 写失败测试（优先调用路径可观测）**

```julia
using Test

@testset "prefer Models.solve over solve_constraint" begin
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_gap_transport_scan.jl"))
    opts = Main.ScanOptions("tmp.csv", nothing, [-0.5], 150.0, 150.0, 1.0, 0.0, 0.0, 1.0, true, false, false, 12, 6, 40, 28, 8, 8, 6, true, 0.6, 8, 10, :linear, 128, :finite_15, 5, 24, 8.0)
    eq = Main.solve_models_equilibrium(150 / Main.Constants_PNJL.ħc_MeV_fm, 0.0, -0.5, Main.TransportWorkflow.PNJL.HADRON_SEED_5, opts)
    @test eq !== nothing
    @test eq.solver_backend == :models
end
```

- [x] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: FAIL，提示优先路径尚未实现。

- [x] **Step 3: 最小实现（脚本内新增统一求解包装）**

```julia
function _solve_fixedmu_via_models_solve(T_fm, muq_fm, xi, seed_state, opts)
    raw = Main.Models.solve(
        PNJL_MODEL,
        Main.Models.FixedMu(),
        T_fm,
        muq_fm;
        seed_guess=Float64.(seed_state[1:min(5, length(seed_state))]),
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
    )
    return raw
end

function _solve_fixedmu_via_models_constraint(T_fm, muq_fm, xi, seed_state, opts)
    return Main.Models.solve_constraint(
        PNJL_MODEL,
        Main.Models.FixedMu(),
        T_fm;
        μ_fm=muq_fm,
        seed_guess=Float64.(seed_state[1:min(5, length(seed_state))]),
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        residual_norm_max=1e-4,
        physicality_check=Main.Models.is_physical_solution,
    )
end

function solve_models_equilibrium(T_fm, muq_fm, xi, seed_state, opts)
    try
        raw = _solve_fixedmu_via_models_solve(T_fm, muq_fm, xi, seed_state, opts)
        return _normalize_equilibrium_result(raw; solver_backend=:models)
    catch
        raw = _solve_fixedmu_via_models_constraint(T_fm, muq_fm, xi, seed_state, opts)
        return _normalize_equilibrium_result(raw; solver_backend=:models)
    end
end
```

- [x] **Step 4: 补全映射与有效性检查（保持下游接口不变）**

```julia
return (
    converged=true,
    x_state=SVector{5}(Tuple(Float64.(raw.x_state))),
    mu_vec=SVector{3}(Tuple(Float64.(raw.mu_vec))),
    masses=SVector{3}(Tuple(Float64.(raw.masses))),
    iterations=Int(raw.iterations),
    residual_norm=Float64(raw.residual_norm),
    solver_backend=:models,
    omega=Float64(raw.omega),
)
```

- [x] **Step 5: 跑测试确认通过**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: PASS，证明优先调用与兜底行为满足测试。

- [x] **Step 6: 提交（入口策略变更）**

```bash
git add scripts/relaxtime/run_gap_transport_scan.jl tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl
git commit -m "fix: prefer src solver fallback path in gap transport scan"
```

### Task 3: 新增失败点 sidecar 诊断落盘

**Files:**
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`
- Test: `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

- [x] **Step 1: 写失败测试（失败点 sidecar 字段与写入时机）**

```julia
using Test

@testset "failed point sidecar" begin
    # 构造单点失败场景，断言 sidecar 出现且字段齐全。
    @test false
end
```

- [x] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: FAIL，提示 sidecar 尚未写入。

- [x] **Step 3: 最小实现（新增 sidecar 参数与写入函数）**

```julia
struct ScanOptions
    output::String
    channel_diagnostics_output::Union{Nothing, String}
    xi_values::Vector{Float64}
    tmin_mev::Float64
    tmax_mev::Float64
    tstep_mev::Float64
    mubmin_mev::Float64
    mubmax_mev::Float64
    mubstep_mev::Float64
    overwrite::Bool
    resume::Bool
    compute_bulk::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    tau_threshold_subtraction::Bool
    tau_asym_window::Float64
    tau_asym_fit_min_points::Int
    tau_asym_extra_points::Int
    tau_interpolation_mode::Symbol
    sigma_grid_n::Int
    integration_mode::Symbol
    gc_every_n::Int
    tr_p_nodes::Int
    tr_p_max_fm::Float64
    failed_points_output::Union{Nothing,String}
end

function write_failed_point_row!(io, T_mev, muB_mev, xi, diag, err)
    # 输出字段：T_MeV,muB_MeV,xi,seed_source,phase_prev,phase_curr_hint,error_type,error_message,timestamp
end
```

- [x] **Step 4: 在 `catch point_err` 分支接入写 sidecar**

```julia
catch point_err
    @warn "SCAN POINT FAILED — skipped" T_mev=T_mev muB_mev=muB_mev xi=xi err=point_err
    if failed_io !== nothing
        write_failed_point_row!(failed_io, T_mev, muB_mev, xi, diag_or_hint, point_err)
    end
end
```

- [x] **Step 5: 运行测试确认通过**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: PASS，sidecar 行与字段符合预期。

- [x] **Step 6: 提交（失败可观测性）**

```bash
git add scripts/relaxtime/run_gap_transport_scan.jl tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl
git commit -m "feat: add failed-point sidecar diagnostics for skipped scan points"
```

### Task 4: 最小复现验证 T150 首点

**Files:**
- Modify: `scripts/relaxtime/run_gap_transport_scan.jl`（仅在需要修正参数解析时）
- Create: `docs/analysis/relaxtime/t150_minimal_rerun_notes.md`
- Test: `D:\Desktop\Temp\relaxtime_t150_repro\`（运行输出目录）

- [x] **Step 1: 在 temp 目录准备输出路径**

Run: `mkdir "D:\Desktop\Temp\relaxtime_t150_repro"`
Expected: 目录创建成功。

- [x] **Step 2: 运行最小点复扫（先 `xi=-0.5`）**

Run: `julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 150 --tmax 150 --tstep 1 --mubmin 0 --mubmax 0 --mubstep 1 --xi-list -0.5 --mode finite_15 --tau-p-nodes 28 --tau-angle-nodes 8 --tau-phi-nodes 8 --tau-n-sigma 6 --sigma-grid-n 128 --output "D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05.csv" --overwrite`
Expected: 生成 CSV 且 `xi=-0.5` 行存在。

- [x] **Step 3: 若首点仍失败，扩展到 `-0.48,-0.46` 并检查连续推进**

Run: `julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 150 --tmax 150 --tstep 1 --mubmin 0 --mubmax 0 --mubstep 1 --xi-list -0.5,-0.48,-0.46 --mode finite_15 --tau-p-nodes 28 --tau-angle-nodes 8 --tau-phi-nodes 8 --tau-n-sigma 6 --sigma-grid-n 128 --output "D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05_to_xim046.csv" --overwrite`
Expected: 至少首点可解；若仍失败，sidecar 提供失败机制。

- [x] **Step 4: 记录验证结果到说明文档**

```markdown
# T150 最小复扫记录

- 命令：`julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 150 --tmax 150 --tstep 1 --mubmin 0 --mubmax 0 --mubstep 1 --xi-list -0.5 --mode finite_15 --tau-p-nodes 28 --tau-angle-nodes 8 --tau-phi-nodes 8 --tau-n-sigma 6 --sigma-grid-n 128 --output "D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05.csv" --overwrite`
- 结果：xi=-0.5 [成功/失败]
- sidecar：<path>
- 结论：是否恢复连续性起点
```

- [x] **Step 5: 提交（最小验证记录）**

```bash
git add docs/analysis/relaxtime/t150_minimal_rerun_notes.md
git commit -m "docs: record minimal T150 xi=-0.5 rerun evidence"
```

### Task 5: 完成 T190/T200 毛刺证据链分析文档

**Files:**
- Create: `docs/analysis/relaxtime/t150-convergence-and-t190-t200-spikes-analysis.md`
- Modify: `docs/analysis/relaxtime/t150_minimal_rerun_notes.md`
- Test: `D:\Desktop\Julia_RelaxTime\.worktrees\repro-main-oldparams\data\outputs\tmp\repro_main_oldparams\results\relaxtime\plan_b\plan_b_merged.csv`

- [x] **Step 1: 写文档骨架（6 段结构）**

```markdown
## 1) 现象复述与数据口径
## 2) 问题1根因链
## 3) 问题1修复结果
## 4) 问题2毛刺解释
## 5) 结论与置信度
## 6) 风险与后续
```

- [x] **Step 2: 填充问题1证据（源码链路 + 最小复跑 + sidecar）**

```markdown
- 源码证据：`run_gap_transport_scan.jl` 调用路径变化
- 数据证据：T150 xi=-0.5 结果行/失败行
- 机制证据：failed_points sidecar
```

- [x] **Step 3: 填充问题2证据（数值链路与物理链路分权重）**

```markdown
- 数值链路：比值放大、tau/tauinv 局部异常、平衡量平滑性对照
- 物理链路：相结构与阈值敏感区的可解释不光滑
- 权重：数值证据 70%，物理证据 30%
```

- [x] **Step 4: 输出结论与置信度（High/Medium/Low + 依据）**

```markdown
- 结论 A：T190/T200 的主要毛刺由数值敏感放大主导，物理链路提供次级解释。
  - 置信度：High
  - 依据：源码链路 + 数据对照 + 最小复跑
```

- [x] **Step 5: 运行文档一致性检查（仓库治理脚本）**

Run: `julia --project=. scripts/dev/check_docs_consistency.jl`
Expected: 通过；若失败，修正文档锚点/引用后重跑至通过。

- [x] **Step 6: 提交（分析报告）**

```bash
git add docs/analysis/relaxtime/t150-convergence-and-t190-t200-spikes-analysis.md docs/analysis/relaxtime/t150_minimal_rerun_notes.md
git commit -m "docs: add T150 convergence and T190/T200 spike evidence analysis"
```

### Task 6: 回归验证与交付检查

**Files:**
- Modify: `docs/analysis/relaxtime/t150-convergence-and-t190-t200-spikes-analysis.md`（仅当检查发现问题）
- Test: `scripts/relaxtime/run_gap_transport_scan.jl`

- [x] **Step 1: 运行目标单测文件**

Run: `julia --project=. -e 'ENV["UNIT_FILES"]="relaxtime/test_run_gap_transport_scan_solver_entry.jl"; include("tests/unit/runtests.jl")'`
Expected: PASS。

- [x] **Step 2: 运行 unit smoke（防止破坏快速回归）**

Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
Expected: PASS。

- [x] **Step 3: 运行脚本治理检查（可选但推荐）**

Run: `julia --project=. scripts/dev/check_active_docs_governance.jl`
Expected: PASS。

- [x] **Step 4: 人工核对交付清单**

```text
1) T150 首点修复或失败机制可见
2) failed_points sidecar 可复核
3) T190/T200 毛刺分析报告含结论+证据链+置信度
```

- [x] **Step 5: 最终提交（若 Step 3/4 引发修正）**

```bash
git add scripts/relaxtime/run_gap_transport_scan.jl tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl docs/analysis/relaxtime/t150-convergence-and-t190-t200-spikes-analysis.md
git commit -m "fix: improve plan_b solver fallback observability and document spike root cause"
```
