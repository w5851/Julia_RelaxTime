# PNJL Maxwell tolerance-contract feasibility v1

verdict: `feasible_candidate`。本目录使用 deep-oracle run `30676440627` 的五个固定低温
anchor 做 solver-free post-processing replay；没有调用 equilibrium/Newton solver，也没有修改
production、reference 或历史 evidence。

扫描的内部 Maxwell bisection tolerance 为 `0.0001, 5.0e-5, 1.0e-5, 5.0e-6`。外层
coarse/fine geometry gate 固定为 position `0.025 MeV`、density `0.0025`、
area `5.0e-5`，因此本审计明确区分“求根停止条件”和“跨 rho 层收敛证书”。
严格候选要求 `5.0e-6` 下两层 residual 均不超过该值、topology 不切换、
geometry gate 通过且二分未耗尽 `60` 次迭代；`1e-5→5e-6` 只作为跨容差稳定性检查。

本 feasibility 不是 production gate，也不等于放宽 Maxwell area 容差。若 verdict 不是
`feasible_candidate`，不得创建 production tolerance-contract PR；应保留逐点曲线和失败原因
供物理/数值审核。
