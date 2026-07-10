---
name: transport-regression-keeper
description: 管理 Julia_RelaxTime transport/relaxtime 改动的测试层选择、固定点回归、数值漂移判断和证据汇报。用于截面、散射率、弛豫时间或输运工作流的具体改动；全仓 baseline 制度建设使用 baseline-regression-governance。
---

# Transport Regression Keeper

## Hard rules

- 不得为了通过测试而先放宽容差。
- 不得无证据刷新 baseline。
- 先判断影响与测试层，再选择命令。
- 触及统一入口或求解链路时，检查 `Models` 契约和非 fixed-μ 联合求解语义。

## Test-layer decision

- 局部公式、验证器或 helper：`tests/unit/relaxtime/`
- 跨模块 transport 工作流：`tests/integration/relaxtime/`
- 固定点数值承诺和 baseline：`tests/regression/relaxtime/`
- 外部物理参考或文献映射：`tests/validation/relaxtime/`

基线存放在 `tests/baselines/relaxtime/`。

## Workflow

1. 记录改动模块、影响输出量和调用入口。
2. 区分局部公式、工作流拼装和统一入口变更。
3. 选择最窄且足够的测试层与单文件 selector。
4. 对比 baseline 或既有固定点，定位变化点和变化量。
5. 诊断算法、物理口径、配置、精度或实现错误来源。
6. 判断漂移是预期修复、可接受误差还是回归退化。
7. 只有证据充分时才更新 baseline，并按通用治理记录差异。

## Evidence report

汇报改动模块、测试层及理由、命令、漂移点位和幅度、容差与物理判断、baseline 决定、未变化契约和残余 validation 风险。

局部 helper 调用链很广时，先跑 unit，再选择一个最短 integration smoke；CSV 或固定点断言受影响时默认运行 regression。
