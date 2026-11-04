# benchmark/relaxtime

本目录用于 relaxtime、散射、传播子和介子热力学性能/数值 oracle 探针。

当前推荐入口：

- `benchmark_quadgk_oracle_smoke.jl`：用隔离 QuadGK 对照当前固定 hybrid B0 实部，并报告 oracle error estimate 与相对差；不设 correctness gate。
- `benchmark_total_cross_section_hotspot.jl`：总截面热点路径。
- `benchmark_transport_cache_stages.jl`：当前 thermal-weighted 缓存的热测度、特征发现和自适应建表分项成本。
- `benchmark_transport_production_entrypoint.jl`：当前默认扫描入口的两点重复与完整计时；读取已审阅的生产基线，保留首次完整预热和各点全部通道诊断。
- `transport_runner_timing.jl`：原样调用 phase-guided CLI，记录加载/扫描的 wall、GC、JIT 与累计分配，以及实际主机和线程信息；`TRANSPORT_RUNNER_TIMING_PATH` 指定新的 TOML 输出。JIT 已包含在各阶段耗时中，完整进程时间与峰值 RSS 仍由外层 GNU time 记录。
- `benchmark_transport_phase_branches.jl`：局部分支续接与通用多种子池的同点对照。
- `benchmark_bulk_viscosity_derivatives.jl`：体黏滞导数路径。
- `bench_bu_double_integral_phase_e3.jl`：BU 双积分点计算。
- `bench_meson_density_ad_diagnostic_chain.jl`：介子数密度 AD 诊断链。

历史 `benchmark_adaptive_nodes*.jl`、`benchmark_with_roots.jl` 和 `comprehensive_benchmark.jl` 保留为旧策略演化记录，不是稳定入口；它们的旧 strategy keyword 不代表当前 API。

旧 w0cdf/PCHIP benchmark 已随旧实现退役。`benchmark_transport_threshold_complete.jl`
包含新旧完整路径对照，重放该历史实验需使用 `complete_timing_p48_v1/source` 快照，
不能在当前源码上调用已删除的旧 API。已验收的完整成本和重跑前检查见
[输运采样指南](../../docs/guides/scripts/transport_sampling_qualification.md)。

环境初始化和显式叠加命令见上级 [benchmark README](../README.md)。
