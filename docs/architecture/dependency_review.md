# 依赖审计报告
生成时间：2026-02-19T17:33:36.818

来源：docs/architecture/dependencies.mmd

## 跨目录依赖清单（汇总）

- relaxtime -> integration: 9
- relaxtime -> utils: 2

## 跨目录依赖明细

- relaxtime/AFieldBuilder.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/AverageScatteringRate.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/MesonMass.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegrals.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegrals.jl -> integration/IntervalQuadratureStrategies.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegralsAniso.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/OneLoopIntegralsAniso.jl -> integration/IntervalQuadratureStrategies.jl (relaxtime -> integration)
- relaxtime/TransportCoefficients.jl -> integration/GaussLegendre.jl (relaxtime -> integration)
- relaxtime/TransportCoefficients.jl -> integration/PhaseSpaceSampling.jl (relaxtime -> integration)
- relaxtime/MesonPropagator.jl -> utils/ParticleSymbols.jl (relaxtime -> utils)
- relaxtime/ScatteringAmplitude.jl -> utils/ParticleSymbols.jl (relaxtime -> utils)

## 违规点（基于依赖矩阵）

- 未发现违规

## 调整建议

- 当前依赖符合矩阵，建议继续保持。

