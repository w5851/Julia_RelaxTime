# Audit boundary

输入文件的 calculation SHA、workflow head SHA、source manifest hash 和本地 producer head
写入 `manifest.json`。曲线只按已有 `(xi,T,rho)` 记录重建；没有插值、补点或 equilibrium
重算。`PhaseCore.maxwell_construction` 使用每个 sweep 的显式 `tol_area`，而
`PhaseGridConvergence._compare_phase_geometry` 使用独立的 outer area gate，并记录两者的
实际值。因此“area residual”与“geometry convergence”均可追溯，且不会因读取旧 slice
summary 而把其 gate 状态冒充为新证书。
