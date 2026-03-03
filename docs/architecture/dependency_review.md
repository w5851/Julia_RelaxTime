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

## 2026-03-02 阶段D补充（Phase 主流程去旧路径）

- 主流程入口（`scripts/pnjl/calculate_phase_structure.jl`、`Models.run_phase_pipeline`）已统一走 `src/models/phase/*`，不再依赖 `src/models/pnjl/analysis/PhaseTransition.jl`。
- `src/models/pnjl/analysis/PhaseTransition.jl` 保留为兼容层：用于历史调用方短期过渡，不作为新流程依赖。
- 迁移约束：新增/修改 phase 业务逻辑应仅落位 `src/models/phase/*`；兼容层仅允许最小转发和兼容提示，不承载新业务编排。

## 2026-03-03 阶段D收尾（删除兼容层）

- 已删除 `src/models/pnjl/analysis/PhaseTransition.jl`，`PNJL` 分析导出直接复用 `src/models/phase/PhaseCore.jl` 与 `src/models/phase/PhaseIO.jl`。
- 已移除仅服务兼容层的单测（`tests/unit/pnjl/test_phase_transition.jl`、`tests/unit/models/test_phase_transition_forwarding_smoke.jl`），主流程 smoke 继续保留。

