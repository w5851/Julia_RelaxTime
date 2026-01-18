# 依赖审计报告
生成时间：2026-01-18T18:32:30.213

来源：docs/architecture/dependencies.mmd

## 跨目录依赖清单（汇总）

- pnjl -> relaxtime: 3
- relaxtime -> integration: 6
- relaxtime -> utils: 1

## 跨目录依赖明细

- pnjl/workflows/TransportWorkflow.jl -> relaxtime/OneLoopIntegrals.jl (pnjl -> relaxtime)
- pnjl/workflows/TransportWorkflow.jl -> relaxtime/RelaxationTime.jl (pnjl -> relaxtime)
- pnjl/workflows/TransportWorkflow.jl -> relaxtime/TransportCoefficients.jl (pnjl -> relaxtime)
- relaxtime/AverageScatteringRate.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegrals.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegrals.jl -> integration/IntervalQuadratureStrategies.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegralsAniso.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegralsAniso.jl -> integration/IntervalQuadratureStrategies.jl (relaxtime -> integration)
- relaxtime/TransportCoefficients.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/ScatteringAmplitude.jl -> utils/ParticleSymbols.jl (relaxtime -> utils)

## 违规点（基于依赖矩阵）

- 未发现违规

## 调整建议

- 当前依赖符合矩阵，建议继续保持。

