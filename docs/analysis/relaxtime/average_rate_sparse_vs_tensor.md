平均散射率：张量积 vs 稀疏网格（实验）

脚本：`tests/analysis/relaxtime/average_rate_sparse_vs_tensor.jl`

用途：在同一过程上对比 Gauss-Legendre 张量积积分与 Smolyak 稀疏网格，在“评估次数相近”条件下的精度表现。

## 如何运行

```bash
julia --project tests/analysis/relaxtime/average_rate_sparse_vs_tensor.jl
```

## 目前的设置

- 过程：`:ssbar_to_uubar`
- 基准：张量积节点 p=10, angle=6, phi=8（作为参考值）
- 对比：张量积节点 p=4, angle=3, phi=3（评估次数≈144，接近稀疏网格）
- 稀疏网格：Smolyak level=3（1D 采用 Gauss-Legendre，阶数=2*level-1，当前评估次数≈131）
- `n_sigma_points=16`；σ(s) 使用 AverageScatteringRate 的拟合缓存；ρ 仍用默认节点求。

## 输出解读

脚本打印三组 ω 以及相对误差（相对于基准）：

- `Matched tensor`：与稀疏网格评估次数相近的张量积网格。
- `Sparse grid`：稀疏网格的相对误差与评估次数。
- `Eval counts`：两者各自的积分点（函数评估）数量，用于对比“相同点数”下的精度。

## 最近一次运行结果（示例）

- Reference (tensor) p=10,6,8 → ω ≈ 2.981622e+00，evals=28800，time≈0.778 s，per-eval≈2.70e-05 s
- Matched tensor p=4,3,3 → ω ≈ 2.980953e+00，rel_err≈2.24e-04，evals=432，time≈0.001 s，per-eval≈2.59e-06 s
- Sparse grid level=3 → ω ≈ 3.025002e+00，rel_err≈1.46e-02，evals=131，time≈0.086 s，per-eval≈6.60e-04 s
- Eval counts: tensor matched 432 vs sparse 131

观察：在相近评估次数下（百级），张量积给出的相对误差更小且单点评估更快（稀疏网格未去重、系数组合带来额外开销）。

## 节点敏感性扫描（张量积，仅比较 vs 参考值）

参考同上（p=10,6,8）。以下保持其余维度固定，仅扫描单一维度：

- 扫 p（angle=4, phi=6）：
	- p=3 → rel_err≈3.10e-02；p=4 → 1.27e-04；p=5 → 4.74e-04；p=6 → 3.05e-05；p=8 → 6.05e-06；p=10 → 4.75e-06。
	- 结论：p 维最敏感，需 4–6 节点以上，8 节点已与参考非常接近。
- 扫 angle（p=6, phi=6）：
	- angle=2 → rel_err≈7.53e-05；angle=3 → 4.09e-05；angle=4 → 3.05e-05；angle=6 → 2.47e-05；angle=8 → 2.33e-05。
	- 结论：angle 维次敏感，3–4 节点即可达 1e-5 量级，相对 p 更快收敛。
- 扫 phi（p=6, angle=4）：
	- phi=2 → rel_err≈3.81e-03；phi=3 → 3.77e-04；phi=4 → 2.27e-05；phi=6 → 3.05e-05；phi=8 → 3.07e-05。
	- 结论：phi 维中等敏感，4 节点已足够；再加节点收益甚微。

综合：p 维需要最多节点，其次 angle，phi 最少；若在百级评估预算内，可优先提高 p 节点数，再适度增加 angle，phi 用 4 左右即可。

## 注意

- 稀疏网格实现为实验性质的 Smolyak 组合，使用非嵌套的 Gauss-Legendre 阶数序列；未做节点去重。
- 参考网格较密但仍有限，若需要更严格的基准可提高 p/angle/phi 节点并相应调整 level。
- 若要测试其他过程或节点/level，可直接修改脚本顶部的参数。