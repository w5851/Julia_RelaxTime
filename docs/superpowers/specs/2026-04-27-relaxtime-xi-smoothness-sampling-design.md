# Relaxtime xi Smoothness Sampling Design

## Scope

- 目标：为 Plan-C 的 xi smoothness 采样与评估流程提供稳定、可复现、可批量执行的脚本链路。
- 范围：参数采样、批量计算、批量绘图、自动评分、人工复核合并。
- 非范围：不在本设计中扩展新的物理模型或修改核心求解器实现。

## Workflow Contract

- 输入主线：采样参数 CSV。
- 中间产物：批量运行 manifest、单样本结果 CSV、批量图像。
- 评估产物：`smoothness_scores.csv`、`smoothness_flags.csv`、`review_queue.csv`。
- 人工决议产物：`smoothness_final.csv`。

## Implementation status

- 已落地：Plan-C 脚本链路已具备端到端执行能力（sampling -> batch -> plot -> evaluate -> merge）。
- 已落地：人工复核接口采用 `manual_review.csv` 覆写自动标签，支持 `reason` 字段引用/转义。
- 已落地：输出结构已拆分为“自动评估结果”与“人工合并结果”，便于审计追踪。

## 脚本清单

- `scripts/relaxtime/generate_xi_smoothness_params.jl`：生成采样参数清单（CSV）。
- `scripts/relaxtime/xi_smoothness_sampling_lib.jl`：采样逻辑库（随机点 + 近相线点）。
- `scripts/relaxtime/run_xi_smoothness_batch.jl`：按参数清单执行批量扫描并生成 `run_manifest.json`。
- `scripts/relaxtime/plot_xi_smoothness_batch.jl`：按 manifest 批量生成图像。
- `scripts/relaxtime/evaluate_xi_smoothness.jl`：自动计算平滑性指标并给出标签。
- `scripts/relaxtime/merge_xi_smoothness_review.jl`：合并自动标签与人工复核结论。

## 运行入口

标准入口命令如下（在仓库根目录执行）：

```bash
julia --project=. scripts/relaxtime/generate_xi_smoothness_params.jl --output data/outputs/results/relaxtime/plan_c/sampling/params_24_seed42.csv
julia --project=. scripts/relaxtime/run_xi_smoothness_batch.jl --params data/outputs/results/relaxtime/plan_c/sampling/params_24_seed42.csv --out-root data/outputs/results/relaxtime/plan_c/sampling --resume
julia --project=. scripts/relaxtime/plot_xi_smoothness_batch.jl --manifest data/outputs/results/relaxtime/plan_c/sampling/run_manifest.json --fig-root data/outputs/figures/relaxtime/plan_c
julia --project=. scripts/relaxtime/evaluate_xi_smoothness.jl --manifest data/outputs/results/relaxtime/plan_c/sampling/run_manifest.json --out-root data/outputs/results/relaxtime/plan_c/smoothness
julia --project=. scripts/relaxtime/merge_xi_smoothness_review.jl --flags data/outputs/results/relaxtime/plan_c/smoothness/smoothness_flags.csv --manual-review data/outputs/results/relaxtime/plan_c/smoothness/manual_review.csv --out data/outputs/results/relaxtime/plan_c/smoothness/smoothness_final.csv
```

## Notes

- 运行手册详见：`docs/dev/active/2026-04-27-xi-smoothness-runbook.md`。
- 本设计文档仅维护入口与治理契约，不重复维护逐参数说明。
