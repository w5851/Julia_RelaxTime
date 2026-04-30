---
name: transport-regression-keeper
description: 管理 Julia_RelaxTime 中 transport / relaxtime 相关改动的回归验证、基线漂移判断与证据汇报。适用于截面、散射率、弛豫时间、输运工作流、固定点基线、smoke/core/full 分层回归。关键词：transport, relaxtime, regression, baseline, drift, tests/regression
---

# transport-regression-keeper

## Purpose

将 transport / relaxtime 改动从“代码改了”收束到“数值影响已定位、测试层已选定、漂移证据已汇报”。

本 Skill 是 `baseline-regression-governance` 在本仓库 `relaxtime` 主链路上的窄化版本，优先服务局部改动的快速判断与最小验证闭环。

## Apply When

- 修改 `src/relaxtime/`、`src/models/` 中会影响 transport 输出的实现。
- 调整截面、散射率、弛豫时间、输运系数工作流。
- 需要判断该补 `unit`、`integration`、`regression` 还是 `validation`。
- 需要解释“数值变了”到底是修复、可接受漂移，还是回归退化。

## Do Not Use

- 任务核心是建立全仓通用 baseline 治理规则时：转 `baseline-regression-governance`。
- 任务核心是性能 profiling 或 benchmark 设计时：优先仓库脚本与 `benchmark/` 流程。

## Hard Rules

- 不得为了让测试通过而先放宽容差。
- 不得无说明更新 baseline 文件。
- 先选测试层，再选命令；不要默认全量跑。
- 若改动触及统一入口或求解链路，必须检查 `Models` 契约是否仍成立。

## Repository Defaults

- 代码入口优先看：`src/models/entrypoints.jl`、`src/relaxtime/`
- 测试层：
  - `tests/unit/relaxtime/`
  - `tests/integration/relaxtime/`
  - `tests/regression/relaxtime/`
  - `tests/validation/relaxtime/`
- 常见基线位置：`tests/baselines/relaxtime/`
- 常见命令：
  - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
  - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
  - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

## Standard Workflow

1. Map impact
   - 记录改动模块、影响输出量、可能波及的脚本或入口。
   - 区分“局部公式变动”“工作流拼装变动”“统一入口变动”。
2. Choose the narrowest sufficient test layer
   - 公式或 helper 级：先 `unit`
   - 跨模块链路：补 `integration`
   - 固定点数值承诺：补 `regression`
   - 对外物理参考或文献映射：补 `validation`
3. Run focused verification
   - 优先单文件或 selector 跑法。
   - 只有在局部结果不足以判断风险时，才扩大到 profile 级测试。
4. Judge drift
   - 对比 baseline、既有测试预期或固定点输出。
   - 先判断漂移来源，再决定是否允许更新基线。
5. Summarize evidence
   - 说明改动影响面、执行命令、变化点位、容差判断、残余风险。

## Decision Heuristics

- 改动只影响局部数值 helper，但调用链广：先跑对应 `unit`，再选 1 个最短 `integration` smoke。
- 改动影响 CSV 输出、固定点表或回归断言：默认补或运行 `regression`。
- 改动同时影响 `legacy` 与 `Models` 桥接：检查现有 bridge smoke，不要只看单边实现。

## Output Contract

汇报至少包含：

- 改动模块
- 选择的测试层与理由
- 实际执行命令
- 是否观察到数值漂移
- 漂移是否在现有容差/物理预期内
- 是否需要更新 baseline、文档或后续 validation
