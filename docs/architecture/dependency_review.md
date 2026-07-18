# 依赖审计报告
生成时间：2026-07-18T17:52:10.377

来源：docs/architecture/dependencies.mmd

## 跨目录依赖清单（汇总）

- relaxtime -> integration: 2

## 跨目录依赖明细

- relaxtime/OneLoopIntegrals.jl -> integration/IntervalQuadratureStrategies.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegralsAniso.jl -> integration/IntervalQuadratureStrategies.jl (relaxtime -> integration)

## 违规点（基于依赖矩阵）

- 未发现违规

## 调整建议

- 当前依赖符合矩阵，建议继续保持。

