# 双端点连续性扫描：历史思路与当前对照

历史 `PNJL.DualBranchScan` 模块已随 PR #61 退役。其思路是从低参数端的强子样初值正向
扫描，再从高参数端的夸克样初值反向扫描；每一步使用前一点的解，在重叠区比较热力学势。
旧模块调用单点求解入口，内部可能继续展开多个种子，因此旧耗时不能代表严格单初值的双端点算法。

当前生产入口见 [Models 求解总览](../models/solver/Overview.md)。2026-11-11 的
[双端点对照](../../analysis/pnjl/algorithmic_feasibility/palc_matched_20261108_analysis/README.md)
使用相同的 Models 残差，每端只给一个初值，每个目标只求根一次，失败后停止该方向；
不调用种子目录。脚本为 `scripts/analysis/pnjl_bifurcation_spike/run_endpoint_comparison.jl`。

固定 T 扫 μ 时，双端点通常能覆盖两个外侧分支，适合比较压力交叉；中间不稳定支需要
另行追踪。固定 T 扫 ρ 时，μ(ρ) 的回弯未必造成状态关于 ρ 的折返，普通参数延拓可以
直接沿完整均匀态曲线推进。性能和分支覆盖的实测结论见上述对照。
