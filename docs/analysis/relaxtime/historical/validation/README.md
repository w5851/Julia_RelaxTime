# Historical Validation and Error-Diagnostic Figure Evidence

## 角色

本目录保存已从 `data/outputs/figures/relaxtime/validation` 迁入的历史诊断图。它归入 `analysis/legacy`，用于保留过去的 Julia/Fortran 对照证据，不是当前 `strict` 图，也不是已经建立的 external validation gate。

## 图像

- `mixed_identity_track_fortran_swap_compare_muB600.png`

图像包含四个面板：在文件名所示 `mu_B=600 MeV` 条件下，对比 Julia 与 Fortran（swap-validated）的 `eta`/`eta'` 质量和宽度随 `T [MeV]` 的变化。图中局部跳变和振荡属于需要保留上下文的历史诊断现象，不应单独解释为收敛、物理相变或当前结果的误差界。

## Provenance 边界

资产 registry 记录该 PNG 的 SHA-256 为 `56bcfa8ec14da8f56d0c46496b9f64e02c68652584bd9a0dc7e2616ca9039ba4`，但没有发现可关联的 `plot_manifest` 或 tracked generator。因而本说明只记录文件名和图面可见的比较对象，不推断未被证据支持的求解器设置、输入数据或数值结论。

相关的当前治理边界见 `docs/guides/sop/workflows/meson_thermodynamics.md`：仓库尚未把这类历史对照自动升级为 external literature validation。迁移不修改该 PNG 及任何结果 CSV。
