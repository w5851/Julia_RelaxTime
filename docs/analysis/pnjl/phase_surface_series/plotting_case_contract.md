# PNJL 三维相图绘图案例合同

公共尺寸、字体、图例及 PNG/PDF 交付流程见
[绘图 SOP](../../../guides/sop/workflows/figure_production.md)。
本页规定该图族的物理筛选和显示边界；具体输入与输出以各 case manifest 为准。

## 物理筛选

- 先筛选物理状态，再做三角化或排版。序参量偏导峰只是 response-peak 候选；
  只有处于同一 `xi` 切片的 CEP 化学势侧、且不落入 Maxwell 一阶区，才可绘为 crossover。
- 同一 `(xi, mu_q)` 不同时标为 crossover 和 Maxwell。`mu_q > mu_CEP` 时，
  响应峰只保留在诊断表，物理面只保留 Maxwell；不得混入物理 crossover 面。
- CEP 只有通过 strict gate 才可标为 confirmed endpoint。只有 bracket 时，
  显示 bracket 或明确的 `estimated_midpoint`，不把中点写成已确认单值。
- 相邻原生采样跨越 CEP 筛选边界时保留 gap，并在 manifest／派生表记录；
  不补线、补点或连接不同物理分支。support、收敛和几何资格分别保留。

## 诊断显示

- `visualization-only closed` 是 `audit` 的显示子模式：可统一 finite/converged
  Maxwell 行的颜色，但 geometry/interpolation 未闭合状态必须保留；任何仅用于
  显示的三角化上限都要声明，不产生新数据，不授予 strict 或 phase-reference 资格。
- 高于 CEP 的 response peak 若需展示，放入独立 diagnostic overlay；
  论文物理候选面排除这些点，不使用灰色叉号将其混入 crossover 面。
- 诊断空洞来源时使用 `diagnostic_no_triangulation`：仅连接原生有序 support
  的相邻线段，超过采样门限的 gap 断开。该视图不填面、不合成点，也不把
  unresolved 行升级为 Maxwell boundary。

诊断结果、历史图版和运行事实保留在各 `analysis/`、`figure_layer/` case；
不由本合同追溯改写已有数据、图像或审核状态。
