# benchmark/relaxtime

本目录用于 relaxtime、散射、传播子和介子热力学性能/数值 oracle 探针。

当前推荐入口：

- `benchmark_quadgk_oracle_smoke.jl`：用隔离 QuadGK 对照当前固定 hybrid B0 实部，并报告 oracle error estimate 与相对差；不设 correctness gate。
- `benchmark_total_cross_section_hotspot.jl`：总截面热点路径。
- `benchmark_average_scattering_rate_N60_w0cdf_pchip.jl`：平均散射率与截面缓存路径。
- `benchmark_bulk_viscosity_derivatives.jl`：体黏滞导数路径。
- `bench_bu_double_integral_phase_e3.jl`：BU 双积分点计算。
- `bench_meson_density_ad_diagnostic_chain.jl`：介子数密度 AD 诊断链。

历史 `benchmark_adaptive_nodes*.jl`、`benchmark_with_roots.jl` 和 `comprehensive_benchmark.jl` 保留为旧策略演化记录，不是稳定入口；它们的旧 strategy keyword 不代表当前 API。

环境初始化和显式叠加命令见上级 [benchmark README](../README.md)。
