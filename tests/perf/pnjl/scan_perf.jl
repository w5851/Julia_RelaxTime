#!/usr/bin/env julia
"""
    scan_perf.jl

PNJL T-μ 和 T-ρ 扫描性能测试。
排除 JIT 编译影响，使用 BenchmarkTools 进行严谨测试。

输出：
- JSON 格式结果：output/perf/pnjl_scan.json
- Markdown 报告：output/perf/pnjl_scan.md
"""

using Printf
using Statistics
using BenchmarkTools
using JSON3
using Dates

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PERF_OUTPUT_DIR = joinpath(PROJECT_ROOT, "tests", "perf", "results", "pnjl")
const JSON_OUTPUT = joinpath(PERF_OUTPUT_DIR, "scan_benchmark.json")
const MARKDOWN_OUTPUT = joinpath(PERF_OUTPUT_DIR, "scan_benchmark.md")

# 加载 PNJL 模块
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .PNJL: solve, FixedMu, FixedRho, ContinuitySeed, update!
using .PNJL.TrhoScan: run_trho_scan

# ============================================================================
# 测试配置
# ============================================================================

# T-μ 扫描配置
const TMU_CONFIG = (
    T_range = 50.0:25.0:200.0,   # 7 个温度点
    μ_range = 0.0:50.0:300.0,    # 7 个化学势点
    xi = 0.0,
)

# T-ρ 扫描配置
const TRHO_CONFIG = (
    T_values = [80.0, 100.0, 120.0, 140.0, 160.0],  # 5 个温度点
    ρ_range = collect(0.0:0.2:3.0),                  # 16 个密度点
    xi = 0.0,
)

const DEFAULT_SAMPLES = 10  # BenchmarkTools 采样次数

# ============================================================================
# T-μ 扫描
# ============================================================================

function run_tmu_scan_benchmark()
    """执行一次完整的 T-μ 扫描"""
    xi = TMU_CONFIG.xi
    mode = FixedMu()
    seed = ContinuitySeed()
    
    n_success = 0
    
    for T_MeV in TMU_CONFIG.T_range
        for μ_MeV in TMU_CONFIG.μ_range
            T_fm = T_MeV / ħc_MeV_fm
            μ_fm = μ_MeV / ħc_MeV_fm
            result = solve(mode, T_fm, μ_fm; xi=xi, seed_strategy=seed)
            if result.converged
                n_success += 1
                update!(seed, result.solution)
            end
        end
    end
    
    return n_success
end

# ============================================================================
# T-ρ 扫描
# ============================================================================

function run_trho_scan_benchmark()
    """执行一次完整的 T-ρ 扫描"""
    output_path = tempname() * ".csv"
    
    result = run_trho_scan(
        T_values = TRHO_CONFIG.T_values,
        rho_values = TRHO_CONFIG.ρ_range,
        xi_values = [TRHO_CONFIG.xi],
        output_path = output_path,
        overwrite = true,
        reverse_rho = true
    )
    
    isfile(output_path) && rm(output_path)
    return result.success
end

# ============================================================================
# 结果结构
# ============================================================================

struct ScanBenchResult
    label::String
    n_points::Int
    stats::NamedTuple
    convergence_rate::Float64
end

function summarize_trial(trial::BenchmarkTools.Trial, n_points::Int)
    times_ms = trial.times ./ 1.0e6
    return (;
        samples = length(times_ms),
        total_min_ms = minimum(times_ms),
        total_median_ms = median(times_ms),
        total_mean_ms = mean(times_ms),
        total_max_ms = maximum(times_ms),
        per_point_min_ms = minimum(times_ms) / n_points,
        per_point_median_ms = median(times_ms) / n_points,
        per_point_mean_ms = mean(times_ms) / n_points,
        allocs = trial.allocs,
        memory_mb = trial.memory / 1024 / 1024,
    )
end

# ============================================================================
# 报告生成
# ============================================================================

function write_reports(results::Vector{ScanBenchResult})
    mkpath(PERF_OUTPUT_DIR)
    timestamp = Dates.format(Dates.now(), Dates.ISODateTimeFormat)
    
    # JSON 输出
    json_payload = (
        generated_at = timestamp,
        samples = DEFAULT_SAMPLES,
        tmu_config = TMU_CONFIG,
        trho_config = TRHO_CONFIG,
        benchmarks = [
            (
                label = r.label,
                n_points = r.n_points,
                convergence_rate = r.convergence_rate,
                stats = r.stats,
            ) for r in results
        ],
    )
    open(JSON_OUTPUT, "w") do io
        JSON3.pretty(io, json_payload)
    end
    
    # Markdown 输出
    open(MARKDOWN_OUTPUT, "w") do io
        println(io, "# PNJL Scan Benchmark")
        println(io, "Generated at: $timestamp")
        println(io, "")
        println(io, "## Configuration")
        println(io, "- T-μ scan: T=$(first(TMU_CONFIG.T_range))-$(last(TMU_CONFIG.T_range)) MeV, μ=$(first(TMU_CONFIG.μ_range))-$(last(TMU_CONFIG.μ_range)) MeV")
        println(io, "- T-ρ scan: T=$(TRHO_CONFIG.T_values) MeV, ρ=$(first(TRHO_CONFIG.ρ_range))-$(last(TRHO_CONFIG.ρ_range)) ρ₀")
        println(io, "- Samples per benchmark: $DEFAULT_SAMPLES")
        println(io, "")
        println(io, "## Results")
        println(io, "")
        println(io, "| Scan Type | Points | Convergence | Min (ms/pt) | Median (ms/pt) | Mean (ms/pt) |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: |")
        for r in results
            s = r.stats
            println(io, "| $(r.label) | $(r.n_points) | $(round(r.convergence_rate * 100, digits=1))% | $(round(s.per_point_min_ms, digits=2)) | $(round(s.per_point_median_ms, digits=2)) | $(round(s.per_point_mean_ms, digits=2)) |")
        end
        println(io, "")
        println(io, "## Memory Usage")
        println(io, "")
        for r in results
            println(io, "- $(r.label): $(round(r.stats.memory_mb, digits=2)) MB")
        end
    end
    
    println("Saved benchmark reports to:")
    println("  JSON -> $JSON_OUTPUT")
    println("  Markdown -> $MARKDOWN_OUTPUT")
end

# ============================================================================
# 主程序
# ============================================================================

function main()
    println()
    println("=" ^ 60)
    println("  PNJL Scan Performance Benchmark")
    println("=" ^ 60)
    println()
    
    results = ScanBenchResult[]
    
    # ========== T-μ 扫描 ==========
    n_tmu = length(TMU_CONFIG.T_range) * length(TMU_CONFIG.μ_range)
    println("T-μ Scan Benchmark")
    println("-" ^ 40)
    @printf("  Config: T=%.0f-%.0f MeV, μ=%.0f-%.0f MeV\n",
            first(TMU_CONFIG.T_range), last(TMU_CONFIG.T_range),
            first(TMU_CONFIG.μ_range), last(TMU_CONFIG.μ_range))
    @printf("  Points: %d\n", n_tmu)
    println()
    
    # 预热
    print("  Warming up... ")
    n_success_warmup = run_tmu_scan_benchmark()
    println("done")
    
    # 基准测试
    print("  Benchmarking ($DEFAULT_SAMPLES samples)... ")
    bench_tmu = @benchmarkable run_tmu_scan_benchmark() evals=1
    trial_tmu = run(bench_tmu; samples=DEFAULT_SAMPLES)
    println("done")
    
    stats_tmu = summarize_trial(trial_tmu, n_tmu)
    convergence_tmu = n_success_warmup / n_tmu
    
    @printf("  Results:\n")
    @printf("    Min:    %.2f ms/point\n", stats_tmu.per_point_min_ms)
    @printf("    Median: %.2f ms/point\n", stats_tmu.per_point_median_ms)
    @printf("    Mean:   %.2f ms/point\n", stats_tmu.per_point_mean_ms)
    @printf("    Convergence: %.1f%%\n", convergence_tmu * 100)
    println()
    
    push!(results, ScanBenchResult("T-μ scan", n_tmu, stats_tmu, convergence_tmu))
    
    # ========== T-ρ 扫描 ==========
    n_trho = length(TRHO_CONFIG.T_values) * length(TRHO_CONFIG.ρ_range)
    println("T-ρ Scan Benchmark")
    println("-" ^ 40)
    @printf("  Config: T=%s MeV, ρ=%.1f-%.1f ρ₀\n",
            join(TRHO_CONFIG.T_values, ","),
            first(TRHO_CONFIG.ρ_range), last(TRHO_CONFIG.ρ_range))
    @printf("  Points: %d\n", n_trho)
    println()
    
    # 预热
    print("  Warming up... ")
    n_success_warmup = run_trho_scan_benchmark()
    println("done")
    
    # 基准测试
    print("  Benchmarking ($DEFAULT_SAMPLES samples)... ")
    bench_trho = @benchmarkable run_trho_scan_benchmark() evals=1
    trial_trho = run(bench_trho; samples=DEFAULT_SAMPLES)
    println("done")
    
    stats_trho = summarize_trial(trial_trho, n_trho)
    convergence_trho = n_success_warmup / n_trho
    
    @printf("  Results:\n")
    @printf("    Min:    %.2f ms/point\n", stats_trho.per_point_min_ms)
    @printf("    Median: %.2f ms/point\n", stats_trho.per_point_median_ms)
    @printf("    Mean:   %.2f ms/point\n", stats_trho.per_point_mean_ms)
    @printf("    Convergence: %.1f%%\n", convergence_trho * 100)
    println()
    
    push!(results, ScanBenchResult("T-ρ scan", n_trho, stats_trho, convergence_trho))
    
    # ========== 输出报告 ==========
    write_reports(results)
    
    println()
    println("=" ^ 60)
    println("  Benchmark Complete")
    println("=" ^ 60)
end

main()
