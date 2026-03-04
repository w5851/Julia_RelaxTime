# benchmark/

性能基准测试目录，使用 `BenchmarkTools.jl` 体系。

## 设计原则

- **关注点分离**：正确性测试在 `tests/`，性能基准在 `benchmark/`
- **可比较**：使用 `@benchmark` 产出 `Trial` 对象，支持 `judge(new, old)` 做 PR 间回归检测
- **依赖隔离**：独立 `Project.toml`，不污染主项目启动时间

## 目录结构

```
benchmark/
├── Project.toml      # BenchmarkTools 依赖
├── benchmarks.jl     # PkgBenchmark 统一入口
├── pnjl/             # PNJL 求解器基准
│   └── bench_*.jl
└── relaxtime/        # 输运/散射基准
    └── bench_*.jl
```

## 运行方式

```bash
# PkgBenchmark（推荐，支持 judge 比较）
julia --project=benchmark -e 'using PkgBenchmark; benchmarkpkg(".")'

# 手动运行单个基准
julia --project=benchmark benchmark/relaxtime/bench_total_cross_section.jl

# 比较两次运行
julia --project=benchmark -e '
    using PkgBenchmark
    old = benchmarkpkg(".", BenchmarkConfig(id="main"))
    new = benchmarkpkg(".", BenchmarkConfig(id="HEAD"))
    judge(new, old)
'
```

## 命名约定

- 基准文件：`bench_<feature>.jl`（前缀 `bench_` 以区分于 `test_`）
- 每个文件向全局 `SUITE["<subsystem>"]` 注册 `@benchmarkable`

## 与 perf 的区别

| | benchmark/ | scripts/perf/ |
|---|---|---|
| 目的 | 监测性能退化 | 诊断性能瓶颈 |
| 工具 | `BenchmarkTools.@benchmark` | `Profile.@profile` / flamegraph |
| 结果 | 可序列化 `Trial` 对象 | SVG / allocation 报告 |
| CI | 可选周期调度 | 手动执行 |
