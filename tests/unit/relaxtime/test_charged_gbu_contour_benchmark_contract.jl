using Test

const CGBU_BENCHMARK_PATH = joinpath(@__DIR__, "..", "..", "..", "benchmark", "relaxtime", "bench_charged_gbu_contour_point.jl")
const CGBU_BENCHMARK_UTILS_PATH = joinpath(@__DIR__, "..", "..", "..", "benchmark", "relaxtime", "charged_gbu_benchmark_utils.jl")

include(CGBU_BENCHMARK_UTILS_PATH)
const CGBU_BENCH_UTILS = Main.ChargedGBUBenchmarkUtils

@testset "Charged GBU benchmark route and timing contract" begin
    source = read(CGBU_BENCHMARK_PATH, String)
    @test Meta.parseall(source) !== nothing
    @test occursin("charged_gbu_infinite_v1", source)
    @test occursin("shell_contribution_inv_fm2", source)
    @test occursin("compile_time", read(CGBU_BENCHMARK_UTILS_PATH, String))
    @test occursin("density_repeats", source)
    @test occursin("FixedMuBConservedCharges", source)
    @test occursin("source_hashes_before", source)
    @test occursin("source_hashes_after", source)
    @test occursin("process_startup_included=false", source)
    @test !occursin("sqrt_s_NN_GeV=NaN", source)
    @test !occursin("CausalGBUResearch", source)
    @test !occursin("thermal=24", source)
end

@testset "Charged GBU benchmark pure aggregation" begin
    points = CGBU_BENCH_UTILS.parse_points("80:719, 120.5:0")
    @test points == [(80.0, 719.0), (120.5, 0.0)]
    for invalid in ("80", "0:719", "NaN:0", "80:-1", "80:Inf", "", "80:719,")
        @test_throws ArgumentError CGBU_BENCH_UTILS.parse_points(invalid)
    end

    sample(ok, elapsed, compile) = (success=ok, value=:value,
        elapsed_s=elapsed, timed_s=elapsed, bytes=10.0, gc_s=0.1,
        compile_s=compile, recompile_s=0.0, error=ok ? nothing : "failed")
    summary = CGBU_BENCH_UTILS.summarize_timed([
        sample(true, 1.0, 0.2), sample(true, 3.0, 0.0), sample(false, 0.0, nothing)])
    @test summary.count == 3
    @test summary.success_count == 2
    @test summary.failure_count == 1
    @test summary.elapsed_s.median == 2.0
    @test summary.compile_s.count == 2
    @test summary.compilation_observed_count == 1
    @test summary.errors == ["failed"]
    @test CGBU_BENCH_UTILS.timing_record(sample(true, 1.0, 0.2)).success
    @test CGBU_BENCH_UTILS.summarize_timed(NamedTuple[]).elapsed_s.mean === nothing
    @test CGBU_BENCH_UTILS.summarize_timed([sample(false, 0.0, nothing)]).success_count == 0
    settings = (mesh=64, cut_nodes=32, tail_nodes=32, omega_nodes=64, q_nodes=8,
        qmax=8.0, repeats=5, density_repeats=3)
    @test CGBU_BENCH_UTILS.validate_settings(settings) === settings
    for bad in (merge(settings, (repeats=0,)), merge(settings, (density_repeats=101,)),
        merge(settings, (qmax=Inf,)), merge(settings, (q_nodes=0,)))
        @test_throws ArgumentError CGBU_BENCH_UTILS.validate_settings(bad)
    end
    result = CGBU_BENCH_UTILS.measure_call(() -> 7)
    @test result.success && result.value == 7
    failed = CGBU_BENCH_UTILS.measure_call(() -> throw(ArgumentError("synthetic failure")))
    @test !failed.success && occursin("synthetic failure", failed.error)
    @test_throws InterruptException CGBU_BENCH_UTILS.measure_call(() -> throw(InterruptException()))
end
