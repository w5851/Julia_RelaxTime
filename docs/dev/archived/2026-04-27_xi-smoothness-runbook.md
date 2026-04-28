---
title: xi smoothness Plan-C 运行手册（Task6）
archived: true
original: docs/dev/active/2026-04-27_xi-smoothness-runbook.md
archived_date: 2026-04-27
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# xi smoothness Plan-C 运行手册（Task6）

## 1. 目标与约定

- 本手册用于一次性串联 xi smoothness 采样、批量计算、批量绘图、自动评估与人工复核合并。
- 运行目录：仓库根目录（`D:\Desktop\Julia_RelaxTime\.worktrees\xi-smoothness-planc`）。
- 默认输出根目录：`data/outputs/results/relaxtime/plan_c`。

## 2. 可直接执行命令

### 2.1 生成参数清单

```bash
julia --project=. scripts/relaxtime/generate_xi_smoothness_params.jl \
  --seed 42 \
  --random-count 12 \
  --near-count 12 \
  --output data/outputs/results/relaxtime/plan_c/sampling/params_24_seed42.csv
```

### 2.2 批量计算

```bash
julia --project=. scripts/relaxtime/run_xi_smoothness_batch.jl \
  --params data/outputs/results/relaxtime/plan_c/sampling/params_24_seed42.csv \
  --out-root data/outputs/results/relaxtime/plan_c/sampling \
  --xi-min -0.5 \
  --xi-max 0.5 \
  --xi-step 0.02 \
  --resume
```

### 2.3 批量绘图

```bash
julia --project=. scripts/relaxtime/plot_xi_smoothness_batch.jl \
  --manifest data/outputs/results/relaxtime/plan_c/sampling/run_manifest.json \
  --fig-root data/outputs/figures/relaxtime/plan_c
```

### 2.4 自动评估

```bash
julia --project=. scripts/relaxtime/evaluate_xi_smoothness.jl \
  --manifest data/outputs/results/relaxtime/plan_c/sampling/run_manifest.json \
  --out-root data/outputs/results/relaxtime/plan_c/smoothness
```

### 2.5 人工复核合并

```bash
julia --project=. scripts/relaxtime/merge_xi_smoothness_review.jl \
  --flags data/outputs/results/relaxtime/plan_c/smoothness/smoothness_flags.csv \
  --manual-review data/outputs/results/relaxtime/plan_c/smoothness/manual_review.csv \
  --out data/outputs/results/relaxtime/plan_c/smoothness/smoothness_final.csv
```

## 3. manual_review.csv 模板

字段约定：`sample_id,field,manual_label,reason`。

- `manual_label` 允许值：`confirm_smooth`、`confirm_not_smooth`、`undecided`。
- `reason` 中包含逗号时请使用 CSV 引号。

```csv
sample_id,field,manual_label,reason
S001,tau_u,confirm_not_smooth,"spike near xi=0.14, recheck local grid"
S001,eta_over_s,undecided,"needs cross-check with figure set A"
S012,sigma_over_T,confirm_smooth,curve monotonic and stable
```

## 4. 产物检查点

- 参数清单：`data/outputs/results/relaxtime/plan_c/sampling/params_24_seed42.csv`
- 运行清单：`data/outputs/results/relaxtime/plan_c/sampling/run_manifest.json`
- 自动评估：
  - `data/outputs/results/relaxtime/plan_c/smoothness/smoothness_scores.csv`
  - `data/outputs/results/relaxtime/plan_c/smoothness/smoothness_flags.csv`
  - `data/outputs/results/relaxtime/plan_c/smoothness/review_queue.csv`
- 人工合并输出：`data/outputs/results/relaxtime/plan_c/smoothness/smoothness_final.csv`
