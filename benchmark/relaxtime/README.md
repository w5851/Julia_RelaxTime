# tests/perf/relaxtime

这一组脚本围绕 **弛豫时间(relaxtime)** 主流程的性能热点：

- 单圈积分与各向异性修正（策略：QUADGK / GL / HYBRID）
- 散射振幅与总截面
- 总传播子与极化函数等中间量

## 脚本速览

- `timing_analysis.jl`：对 `OneLoopIntegralsCorrection.tilde_B0_correction_k_positive` 做分项计时（策略/节点数/节点变换开销）。
- `performance_comparison.jl`：策略 × 节点数的误差/用时对比表。
- `comprehensive_benchmark.jl`：更大参数集的精度+用时综合对比。
- `compare_quadgk_vs_hybrid.jl`：QuadGK(近似“无限制”) vs HYBRID(n=32) 对比。
- `profile_hybrid_overhead.jl` / `profile_hybrid_v2.jl`：拆分 HYBRID 内部阶段（根查找/区间构建/节点生成/被积函数求值/完整调用）。

- `test_total_propagator_performance.jl`：总传播子计算相关性能测试（并对比极化函数）。
- `test_total_cross_section_performance.jl`：总截面完整链路性能测试（可选“完整计算”开关）。
- `test_average_scattering_rate_performance.jl`：平均散射率的轻量性能烟囱测试。

## 运行

从仓库根目录：

- `julia --project=. tests/perf/relaxtime/timing_analysis.jl`
- `julia --project=. tests/perf/relaxtime/performance_comparison.jl`
- `julia --project=. tests/perf/relaxtime/test_total_cross_section_performance.jl`

如果脚本会生成 Markdown/JSON，一般落在同目录或 `tests/perf/results/`（以各脚本 docstring 为准）。
