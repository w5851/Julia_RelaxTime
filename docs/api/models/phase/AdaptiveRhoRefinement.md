# 自适应 rho 加密

本文档描述相图主题中用于提升 CEP 与 Maxwell 稳定性的自适应 `rho` 加密辅助层。核心实现位于 [src/models/phase/AdaptiveRhoRefinement.jl](src/models/phase/AdaptiveRhoRefinement.jl#L3)。

## 适用场景

当已有 `μ(ρ)` 样本在某些区段出现低斜率、稀疏采样或 CEP 附近判据不稳定时，可以先生成补样建议，再让上层重新执行扫描。

当前核心接口：

- `AdaptiveRhoConfig`
- `suggest_refinement_points`
- `merge_rho_values`

这些接口是数据层辅助，而不是独立的求解器入口。

## `AdaptiveRhoConfig`

结构体定义位于 [src/models/phase/AdaptiveRhoRefinement.jl](src/models/phase/AdaptiveRhoRefinement.jl#L5)。构造器位于 [src/models/phase/AdaptiveRhoRefinement.jl](src/models/phase/AdaptiveRhoRefinement.jl#L12)。

```julia
AdaptiveRhoConfig(; slope_tol=5.0, min_gap=0.002, max_points=64, digits=6)
```

关键参数：

- `slope_tol`
- `min_gap`
- `max_points`
- `digits`

## `suggest_refinement_points`

实现位于 [src/models/phase/AdaptiveRhoRefinement.jl](src/models/phase/AdaptiveRhoRefinement.jl#L20)。

该函数根据现有 `(ρ, μ)` 样本估计需要补样的区段，并返回新的 `rho` 候选点。

它适合用于：

- CEP 直接策略下的局部补样
- 低斜率区段的二次扫描
- pipeline 内的自适应 `rho` 扩充

## `merge_rho_values`

实现位于 [src/models/phase/AdaptiveRhoRefinement.jl](src/models/phase/AdaptiveRhoRefinement.jl#L50)。

该函数负责把新增 `rho` 点与原始网格合并，并按给定精度去重与排序。

## 示例

```julia
curve = [(0.10, 325.0), (0.15, 323.0), (0.25, 323.2), (0.40, 335.0)]
rho_vals = [p[1] for p in curve]
mu_vals = [p[2] for p in curve]

config = AdaptiveRhoConfig(; slope_tol=3.0, min_gap=0.01)
extra_rho = suggest_refinement_points(rho_vals, mu_vals; config=config)
rho_grid = merge_rho_values(rho_vals, extra_rho; digits=config.digits)
```

## 与主流程的关系

相图 pipeline 在 CEP 直接策略下会调用该模块，集成点位于 [src/models/phase/PhasePipeline.jl](src/models/phase/PhasePipeline.jl#L181) 与 [src/models/phase/PhasePipeline.jl](src/models/phase/PhasePipeline.jl#L209)。

因此，这一层的职责是为上层重新计算提供更好的网格，而不是直接返回相变点。

production 模式还会在粗网格已经判为 `valid` 或 `invalid` 时执行至少一层嵌套细网格比较。
这与“根据低斜率区段提出补点”的辅助算法是两个互补层次：前者用 Maxwell/spinodal 几何量误差
决定该温度切片能否作为正式结果，后者负责在单层扫描内部改善采样。启用
`rho_geometry_convergence=true` 时，`cep_max_refine_level` 至少为 1；若达到最大层仍不能满足
`rho_position_tol_MeV`、`rho_density_tol` 与 `rho_maxwell_area_tol`，该切片记为 `unknown`。

粗细网格的实际误差、层数和原因写入 `phase_grid_convergence.csv`，不能通过放宽 CEP 分类阈值来掩盖网格未收敛。

如果你想先理解完整主流程，应优先阅读 [Overview.md](docs/api/models/phase/Overview.md)。
