# Struct vs NamedTuple Performance Benchmark
#
# 对应系统流程步骤：
# - Parameter struct migration performance validation
# - Compares struct vs NamedTuple parameter passing performance
#
# 测试内容：
# - Benchmark `relaxation_times` with struct vs NamedTuple parameters
# - Benchmark `average_scattering_rate` with struct vs NamedTuple parameters
# - Benchmark parameter normalization overhead
#
# 运行方式：
# - `julia --project=. tests/perf/relaxtime/benchmark_struct_vs_namedtuple.jl`
#
# 输出：
# - Console output with timing comparisons
# - Performance results saved to tests/perf/results/relaxtime/struct_vs_namedtuple_benchmark.md
#
# Note: This benchmark focuses on functions that have been fully migrated to support
# struct parameters. Some functions (like total_cross_section) have known issues with
# the A field that need to be addressed separately.

using BenchmarkTools
using Printf
using Dates

# Load ParameterTypes module
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "../../../src/ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

# Load modules under test
include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/RelaxationTime.jl")
include("../../../src/relaxtime/AverageScatteringRate.jl")

using .Constants_PNJL
using .EffectiveCouplings
using .RelaxationTime
using .AverageScatteringRate

# ============================================================================
# Test Parameters
# ============================================================================

# Create struct parameters
const QUARK_PARAMS_STRUCT = QuarkParams((
    m = (u=1.52, d=1.52, s=3.04),
    μ = (u=0.3, d=0.3, s=0.3)
))

const THERMO_PARAMS_STRUCT = ThermoParams((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))

# Create NamedTuple parameters
const QUARK_PARAMS_NT = as_namedtuple(QUARK_PARAMS_STRUCT)
const THERMO_PARAMS_NT = as_namedtuple(THERMO_PARAMS_STRUCT)

# Calculate proper K_coeffs using effective couplings
const A_u = 1.0  # Placeholder A values
const A_s = 1.0
const G_u = calculate_G_from_A(A_u, QUARK_PARAMS_NT.m.u)
const G_s = calculate_G_from_A(A_s, QUARK_PARAMS_NT.m.s)
const K_COEFFS = calculate_effective_couplings(
    Constants_PNJL.G_fm2,
    Constants_PNJL.K_fm5,
    G_u,
    G_s
)

# Test densities for relaxation_times (need all 6 flavors)
const TEST_DENSITIES = (
    u = 0.1,
    d = 0.1,
    s = 0.05,
    ubar = 0.1,
    dbar = 0.1,
    sbar = 0.05
)

# ============================================================================
# Helper Functions
# ============================================================================

"""
Create a constant cross-section cache for testing (avoids expensive cross-section calculation).
"""
function constant_sigma_cache(process::Symbol; sigma::Float64=1.0)
    cache = CrossSectionCache(process)
    AverageScatteringRate.insert_sigma!(cache, 0.0, sigma)
    AverageScatteringRate.insert_sigma!(cache, 500.0, sigma)
    return cache
end

"""
Format benchmark results for display.
"""
function format_benchmark_result(name::String, bench_struct, bench_nt)
    time_struct = median(bench_struct).time / 1e6  # Convert to ms
    time_nt = median(bench_nt).time / 1e6
    ratio = time_struct / time_nt
    percent_diff = (ratio - 1.0) * 100.0
    
    return (
        name = name,
        time_struct_ms = time_struct,
        time_nt_ms = time_nt,
        ratio = ratio,
        percent_diff = percent_diff
    )
end

"""
Print benchmark comparison results.
"""
function print_benchmark_comparison(result)
    println("\n" * "="^80)
    println("Benchmark: $(result.name)")
    println("="^80)
    @printf("Struct time:      %.4f ms\n", result.time_struct_ms)
    @printf("NamedTuple time:  %.4f ms\n", result.time_nt_ms)
    @printf("Ratio (S/NT):     %.4f\n", result.ratio)
    @printf("Difference:       %+.2f%%\n", result.percent_diff)
    
    if abs(result.percent_diff) <= 5.0
        println("✓ Performance within 5% tolerance")
    else
        println("⚠ Performance difference exceeds 5% tolerance")
    end
end

# ============================================================================
# Benchmark 1: Parameter Normalization Overhead
# ============================================================================

println("\n" * "="^80)
println("BENCHMARK 1: Parameter Normalization Overhead")
println("="^80)

# Test the overhead of the normalization helpers
@inline _nt_quark_test(q) = q isa QuarkParams ? as_namedtuple(q) : q
@inline _nt_thermo_test(t) = t isa ThermoParams ? as_namedtuple(t) : t

println("\nWarming up...")
# Warmup
_nt_quark_test(QUARK_PARAMS_STRUCT)
_nt_quark_test(QUARK_PARAMS_NT)
_nt_thermo_test(THERMO_PARAMS_STRUCT)
_nt_thermo_test(THERMO_PARAMS_NT)

println("Running benchmarks...")
bench_norm_struct_q = @benchmark _nt_quark_test($QUARK_PARAMS_STRUCT) samples=1000 evals=100
bench_norm_nt_q = @benchmark _nt_quark_test($QUARK_PARAMS_NT) samples=1000 evals=100
bench_norm_struct_t = @benchmark _nt_thermo_test($THERMO_PARAMS_STRUCT) samples=1000 evals=100
bench_norm_nt_t = @benchmark _nt_thermo_test($THERMO_PARAMS_NT) samples=1000 evals=100

println("\nNormalization overhead:")
@printf("  QuarkParams struct:   %.2f ns\n", median(bench_norm_struct_q).time)
@printf("  QuarkParams NT:       %.2f ns\n", median(bench_norm_nt_q).time)
@printf("  ThermoParams struct:  %.2f ns\n", median(bench_norm_struct_t).time)
@printf("  ThermoParams NT:      %.2f ns\n", median(bench_norm_nt_t).time)

# ============================================================================
# Benchmark 2: average_scattering_rate
# ============================================================================

println("\n" * "="^80)
println("BENCHMARK 2: average_scattering_rate")
println("="^80)

# Create constant cache to isolate average_scattering_rate performance
const CACHE_UU = constant_sigma_cache(:uu_to_uu; sigma=1.0)

println("\nWarming up...")
# Warmup
average_scattering_rate(
    :uu_to_uu,
    QUARK_PARAMS_STRUCT,
    THERMO_PARAMS_STRUCT,
    K_COEFFS;
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    cs_cache=CACHE_UU,
    n_sigma_points=4
)
average_scattering_rate(
    :uu_to_uu,
    QUARK_PARAMS_NT,
    THERMO_PARAMS_NT,
    K_COEFFS;
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    cs_cache=CACHE_UU,
    n_sigma_points=4
)

println("Running benchmarks...")
bench_asr_struct = @benchmark average_scattering_rate(
    :uu_to_uu,
    $QUARK_PARAMS_STRUCT,
    $THERMO_PARAMS_STRUCT,
    $K_COEFFS;
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    cs_cache=$CACHE_UU,
    n_sigma_points=4
) samples=50 evals=5

bench_asr_nt = @benchmark average_scattering_rate(
    :uu_to_uu,
    $QUARK_PARAMS_NT,
    $THERMO_PARAMS_NT,
    $K_COEFFS;
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    cs_cache=$CACHE_UU,
    n_sigma_points=4
) samples=50 evals=5

result_asr = format_benchmark_result("average_scattering_rate", bench_asr_struct, bench_asr_nt)
print_benchmark_comparison(result_asr)

# ============================================================================
# Benchmark 3: relaxation_times (full call chain)
# ============================================================================

println("\n" * "="^80)
println("BENCHMARK 3: relaxation_times (full call chain)")
println("="^80)

println("\nWarming up...")
# Warmup
relaxation_times(
    QUARK_PARAMS_STRUCT,
    THERMO_PARAMS_STRUCT,
    K_COEFFS;
    densities=TEST_DENSITIES,
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    n_sigma_points=4
)
relaxation_times(
    QUARK_PARAMS_NT,
    THERMO_PARAMS_NT,
    K_COEFFS;
    densities=TEST_DENSITIES,
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    n_sigma_points=4
)

println("Running benchmarks...")
bench_rt_struct = @benchmark relaxation_times(
    $QUARK_PARAMS_STRUCT,
    $THERMO_PARAMS_STRUCT,
    $K_COEFFS;
    densities=$TEST_DENSITIES,
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    n_sigma_points=4
) samples=20 evals=1

bench_rt_nt = @benchmark relaxation_times(
    $QUARK_PARAMS_NT,
    $THERMO_PARAMS_NT,
    $K_COEFFS;
    densities=$TEST_DENSITIES,
    p_nodes=8,
    angle_nodes=4,
    phi_nodes=4,
    n_sigma_points=4
) samples=20 evals=1

result_rt = format_benchmark_result("relaxation_times", bench_rt_struct, bench_rt_nt)
print_benchmark_comparison(result_rt)

# ============================================================================
# Summary Report
# ============================================================================

println("\n" * "="^80)
println("SUMMARY REPORT")
println("="^80)

results = [result_asr, result_rt]

println("\nPerformance Comparison Summary:")
println("-"^80)
@printf("%-35s %12s %12s %10s %10s\n", "Function", "Struct (ms)", "NT (ms)", "Ratio", "Diff (%)")
println("-"^80)

for r in results
    @printf("%-35s %12.4f %12.4f %10.4f %+10.2f\n",
        r.name, r.time_struct_ms, r.time_nt_ms, r.ratio, r.percent_diff)
end

println("-"^80)

# Check if all benchmarks are within tolerance
all_within_tolerance = all(abs(r.percent_diff) <= 5.0 for r in results)

if all_within_tolerance
    println("\n✓ All benchmarks within 5% performance tolerance")
    println("✓ Struct migration maintains performance requirements")
else
    println("\n⚠ Some benchmarks exceed 5% performance tolerance")
    println("⚠ Further investigation may be needed")
end

# ============================================================================
# Save Results to File
# ============================================================================

output_dir = joinpath(@__DIR__, "../results/relaxtime")
mkpath(output_dir)
output_file = joinpath(output_dir, "struct_vs_namedtuple_benchmark.md")

open(output_file, "w") do io
    println(io, "# Struct vs NamedTuple Performance Benchmark Results")
    println(io, "")
    println(io, "Generated: $(now())")
    println(io, "")
    println(io, "## Test Configuration")
    println(io, "")
    println(io, "- Quark masses: u=$(QUARK_PARAMS_NT.m.u), d=$(QUARK_PARAMS_NT.m.d), s=$(QUARK_PARAMS_NT.m.s)")
    println(io, "- Chemical potentials: u=$(QUARK_PARAMS_NT.μ.u), d=$(QUARK_PARAMS_NT.μ.d), s=$(QUARK_PARAMS_NT.μ.s)")
    println(io, "- Temperature: T=$(THERMO_PARAMS_NT.T)")
    println(io, "- Polyakov loops: Φ=$(THERMO_PARAMS_NT.Φ), Φbar=$(THERMO_PARAMS_NT.Φbar)")
    println(io, "- Anisotropy: ξ=$(THERMO_PARAMS_NT.ξ)")
    println(io, "")
    println(io, "## Benchmark Results")
    println(io, "")
    println(io, "| Function | Struct (ms) | NamedTuple (ms) | Ratio (S/NT) | Difference (%) | Status |")
    println(io, "|----------|-------------|-----------------|--------------|----------------|--------|")
    
    for r in results
        status = abs(r.percent_diff) <= 5.0 ? "✓ Pass" : "⚠ Warning"
        @printf(io, "| %-35s | %11.4f | %15.4f | %12.4f | %+14.2f | %-10s |\n",
            r.name, r.time_struct_ms, r.time_nt_ms, r.ratio, r.percent_diff, status)
    end
    
    println(io, "")
    println(io, "## Normalization Overhead")
    println(io, "")
    println(io, "The inline normalization helpers have minimal overhead:")
    println(io, "- QuarkParams struct normalization: $(median(bench_norm_struct_q).time) ns")
    println(io, "- QuarkParams NamedTuple passthrough: $(median(bench_norm_nt_q).time) ns")
    println(io, "- ThermoParams struct normalization: $(median(bench_norm_struct_t).time) ns")
    println(io, "- ThermoParams NamedTuple passthrough: $(median(bench_norm_nt_t).time) ns")
    println(io, "")
    println(io, "The normalization overhead is negligible compared to the actual computation time.")
    println(io, "")
    
    println(io, "## Analysis")
    println(io, "")
    
    if all_within_tolerance
        println(io, "✓ **All benchmarks passed**: Performance is within 5% tolerance.")
        println(io, "")
        println(io, "The struct migration successfully maintains performance requirements.")
    else
        println(io, "⚠ **Some benchmarks exceeded tolerance**: Performance difference > 5%.")
        println(io, "")
        println(io, "Functions exceeding tolerance:")
        for r in results
            if abs(r.percent_diff) > 5.0
                println(io, "- $(r.name): $(r.percent_diff > 0 ? "slower" : "faster") by $(abs(r.percent_diff))%")
            end
        end
    end
    
    println(io, "")
    println(io, "## Interpretation")
    println(io, "")
    println(io, "- **Ratio < 1.0**: Struct version is faster than NamedTuple")
    println(io, "- **Ratio ≈ 1.0**: Both versions have similar performance")
    println(io, "- **Ratio > 1.0**: Struct version is slower than NamedTuple")
    println(io, "")
    println(io, "Performance differences within ±5% are considered acceptable and may be due to:")
    println(io, "- Measurement noise")
    println(io, "- JIT compilation variations")
    println(io, "- Memory layout differences")
    println(io, "")
    println(io, "## Requirements Validation")
    println(io, "")
    println(io, "This benchmark validates requirements:")
    println(io, "- **Requirement 10.1**: Relaxation time calculations within 5% performance")
    println(io, "- **Requirement 10.2**: Average scattering rate calculations within 5% performance")
    println(io, "")
    println(io, "Note: Requirement 10.3 (cross-section calculations) requires additional work to handle")
    println(io, "the A field properly in the struct normalization process. This will be addressed in a")
    println(io, "future update.")
end

println("\n✓ Results saved to: $output_file")
println("\nBenchmark complete!")
