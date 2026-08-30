# C2 phase surfaces v7：crossover 派生补全层

本包以 v6 为唯一输入，生成第七版后处理结果。它不重跑 C2、不调用 equilibrium solver，也不写入 `data/reference`。v6 的 native crossover、Maxwell 和 spinodal 行均保留；本包新增的行明确标记为 `interpolated_noncertified` 或 `boundary_constrained_endpoint_interpolated_noncertified`。

固定 calculation SHA 为 `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`，v6 manifest SHA-256 为 `8f19f73f7654663504dfa27a9a7e338df13e090bfc7501b8c201742defe4cb0b`。同一 xi 切片内的原生支持间隔用 `4` 等分的分段线性插值填充；相邻 xi 只在共同 native crossover 支持区间内生成 midpoint xi。超过共同支持的尾部保留为 unresolved，不做普通外推。

CEP 采用 `estimated_midpoint` 作为派生层的边界估计：它只用于把 crossover 终止在物理 CEP 边界，不等同于 strict CEP 单值求解；完整 bracket 保存在 `tables/cep_boundary_estimates.csv`。任何 `mu_q > mu_CEP` 的响应峰都不会被重新标成 crossover。

Maxwell 面在本包只复制 v6 native rows，包含原有 `grid_unresolved`/geometry 状态；没有用插值伪造 `maxwell_area` 或 geometry certificate。Maxwell 近 CEP 补点是独立后续任务。

两张图分别用于完整派生显示和层级审计。完整显示图用同一橙色系表示 crossover，但派生层用半透明/虚线；审计图额外区分同 xi 与 xi 方向派生。两张图均不使用三角化。
