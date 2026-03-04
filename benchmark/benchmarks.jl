"""
benchmark/benchmarks.jl — PkgBenchmark 入口

运行方式（推荐）：
    julia --project=benchmark -e 'using PkgBenchmark; benchmarkpkg(".")'

手动运行：
    julia --project=benchmark benchmark/benchmarks.jl

设计原则：
- 每个子系统一个 BenchmarkGroup
- 使用 @benchmarkable 注册，PkgBenchmark 自动收集
- 结果可序列化比较：judge(new_results, old_results)
"""

using BenchmarkTools

const SUITE = BenchmarkGroup()

# ── Subsystem groups ──────────────────────────────────────────────

SUITE["relaxtime"] = BenchmarkGroup()
SUITE["pnjl"] = BenchmarkGroup()

# ── Load subsystem benchmarks ────────────────────────────────────

const BENCH_DIR = @__DIR__

# relaxtime benchmarks
relaxtime_dir = joinpath(BENCH_DIR, "relaxtime")
if isdir(relaxtime_dir)
    for f in sort(readdir(relaxtime_dir; join=true))
        endswith(f, ".jl") && startswith(basename(f), "bench_") && include(f)
    end
end

# pnjl benchmarks
pnjl_dir = joinpath(BENCH_DIR, "pnjl")
if isdir(pnjl_dir)
    for f in sort(readdir(pnjl_dir; join=true))
        endswith(f, ".jl") && startswith(basename(f), "bench_") && include(f)
    end
end
