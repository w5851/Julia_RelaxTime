# Mott Phase mode_ab Merge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new `mode_ab` plotting output that produces two comparison figures (K/pi groups) using only `xi=-0.3,0,0.3`, while preserving existing `mode_a` and `mode_b` behavior.

**Architecture:** Keep the current split architecture: Julia orchestrates, Python renders. Extend `scripts/relaxtime/run_mott_phase_plot_modes.jl` with a focused `_mode_ab` path that reuses existing `_run_python` behavior and copy/rename conventions. Protect existing outputs by adding smoke assertions before implementation and keeping old assertions intact.

**Tech Stack:** Julia 1.10+ scripts, existing Python matplotlib pipeline (`scripts/plot_scan_csv.py`), Julia `Test` integration smoke runner.

---

## File Responsibility Map

- Modify: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`
  - Responsibility: enforce expected output contract for plot mode script; add new `mode_ab` artifact assertions.
- Modify: `scripts/relaxtime/run_mott_phase_plot_modes.jl`
  - Responsibility: orchestrate plotting modes and file export. Add `_mode_ab` without changing `_mode_a`/`_mode_b` contracts.
- Verify (generated outputs): `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/`
  - Responsibility: manual artifact verification on real dataset run.

### Task 1: Lock expected behavior with failing smoke assertions

**Files:**
- Modify: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`
- Test: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`

- [ ] **Step 1: Add failing assertions for new `mode_ab` outputs**

Insert these assertions after existing `mode_b` assertions:

```julia
    @test isfile(joinpath(out_dir, "mode_ab", "mott_mode_ab__M_K__xi3.png"))
    @test isfile(joinpath(out_dir, "mode_ab", "mott_mode_ab__M_pi__xi3.png"))
```

- [ ] **Step 2: Run the single integration smoke file and verify failure**

Run:

```bash
julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_mott_phase_plot_modes_smoke.jl"; include("tests/integration/runtests.jl")'
```

Expected: FAIL because `mode_ab` files are not produced yet.

- [ ] **Step 3: Commit test-only change**

Run:

```bash
git add tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl
git commit -m "test: define mode_ab smoke output contract"
```

### Task 2: Implement `mode_ab` plotting in existing Julia orchestrator

**Files:**
- Modify: `scripts/relaxtime/run_mott_phase_plot_modes.jl`
- Test: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`

- [ ] **Step 1: Add constants for xi targets and grouped observables**

Near `OBSERVABLES`, add:

```julia
const MODE_AB_XIS = (-0.3, 0.0, 0.3)
const MODE_AB_GROUPS = [
    ("M_K", ["M_K", "Gamma_K", "M_u_plus_M_s"]),
    ("M_pi", ["M_pi", "Gamma_pi", "M_u_plus_M_d"]),
]
```

- [ ] **Step 2: Add xi-filter helper for stable float matching**

Add helper:

```julia
@inline function _is_mode_ab_xi(x::Real)
    xf = Float64(x)
    any(t -> isapprox(xf, t; atol=1e-12, rtol=0.0), MODE_AB_XIS)
end
```

- [ ] **Step 3: Add `_mode_ab` implementation**

Add function below `_mode_b`:

```julia
function _mode_ab(input_csv::String, out_root::String)
    mode_ab_dir = joinpath(out_root, "mode_ab")
    mkpath(mode_ab_dir)

    for (label, ys_list) in MODE_AB_GROUPS
        ys = join(ys_list, ',')
        tmp_out = mktempdir()
        args = String[
            PLOT_SCRIPT,
            "--mode", "lines",
            "--csv", input_csv,
            "--x", "T_MeV",
            "--ys", ys,
            "--group", "xi",
            "--out-dir", tmp_out,
        ]
        _run_python(args)

        selected_png = nothing
        for xi in MODE_AB_XIS
            src = joinpath(tmp_out, "$(ys_list[1])_vs_T_MeV.png")
            if isfile(src)
                selected_png = src
                break
            end
            split_file = joinpath(tmp_out, "xi=$(xi)", "multi_y_$(join(ys_list, '_'))_vs_T_MeV.png")
            if isfile(split_file)
                selected_png = split_file
                break
            end
        end

        if selected_png === nothing
            pngs = filter(f -> endswith(lowercase(f), ".png"), readdir(tmp_out; join=true))
            isempty(pngs) && continue
            selected_png = pngs[1]
        end

        dst = joinpath(mode_ab_dir, "mott_mode_ab__$(label)__xi3.png")
        cp(String(selected_png), dst; force=true)
    end
end
```

- [ ] **Step 4: Ensure xi filtering is actually applied before plotting**

In `_mode_ab`, replace direct `input_csv` usage with a filtered temporary CSV generation flow:

```julia
function _write_mode_ab_filtered_csv(input_csv::String)
    out_csv = joinpath(mktempdir(), "mode_ab_filtered.csv")
    open(input_csv, "r") do src
        open(out_csv, "w") do dst
            header_seen = false
            xi_idx = 0
            for line in eachline(src)
                s = strip(line)
                if isempty(s)
                    continue
                elseif startswith(s, "#")
                    println(dst, line)
                    continue
                elseif !header_seen
                    cols = [strip(x) for x in split(s, ',')]
                    xi_idx = findfirst(==("xi"), cols)
                    xi_idx === nothing && throw(ArgumentError("input csv missing xi column"))
                    println(dst, line)
                    header_seen = true
                    continue
                end

                parts = split(s, ',')
                xi_idx > length(parts) && continue
                xv = tryparse(Float64, strip(parts[xi_idx]))
                xv === nothing && continue
                _is_mode_ab_xi(xv) || continue
                println(dst, line)
            end
        end
    end
    return out_csv
end
```

Then in `_mode_ab`, call:

```julia
    filtered_csv = _write_mode_ab_filtered_csv(input_csv)
```

and pass `filtered_csv` to Python args.

- [ ] **Step 5: Wire `main()` to run `_mode_ab` after existing modes**

In `main()`, keep existing calls and append:

```julia
    _mode_ab(input_csv, out_dir)
```

- [ ] **Step 6: Run single smoke integration file and verify pass**

Run:

```bash
julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_mott_phase_plot_modes_smoke.jl"; include("tests/integration/runtests.jl")'
```

Expected: PASS, including old `mode_a`/`mode_b` assertions and new `mode_ab` assertions.

- [ ] **Step 7: Commit implementation**

Run:

```bash
git add scripts/relaxtime/run_mott_phase_plot_modes.jl tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl
git commit -m "feat: add mode_ab merged mott plotting outputs"
```

### Task 3: Real dataset artifact verification in main workspace

**Files:**
- Modify: none
- Verify output: `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/`

- [ ] **Step 1: Run plot script on existing derived CSV (no data regeneration)**

Run:

```bash
julia --project=. scripts/relaxtime/run_mott_phase_plot_modes.jl --in data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv --out-dir data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine
```

Expected: script completes and prints output root path.

- [ ] **Step 2: Verify required files exist**

Check these files exist:

```text
data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_K__xi3.png
data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_pi__xi3.png
```

- [ ] **Step 3: Verify no regressions in old outputs**

Confirm representative old files still exist:

```text
data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_a/mott_mode_a__xi0.png
data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_b/mott_mode_b__M_K.png
```

- [ ] **Step 4: Commit any generated figure updates only if repo policy requires tracking them**

Run:

```bash
git status
```

If figure artifacts are tracked and changed, include them in a separate commit with message:

```bash
git add data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_K__xi3.png data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_pi__xi3.png
git commit -m "docs: refresh mott mode_ab reference figures"
```

If artifacts are untracked/ignored by policy, skip commit and record verification in PR notes.

## Final verification checklist

- [ ] `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl` passes.
- [ ] `mode_ab` folder contains exactly the two agreed merged outputs.
- [ ] Existing `mode_a` and `mode_b` behavior remains available.
- [ ] No scan rerun or data regeneration performed.
