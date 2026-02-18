#!/usr/bin/env julia

using JSON3

const ROOT = pwd()
const SINGLE_POINT_JSON = joinpath(ROOT, "tests", "perf", "results", "pnjl", "single_point_benchmark.json")
const SCAN_JSON = joinpath(ROOT, "tests", "perf", "results", "pnjl", "scan_benchmark.json")

const DEFAULT_SINGLE_POINT_MEDIAN_MS = 2500.0
const SINGLE_POINT_MEDIAN_LIMITS = Dict(
    "Newton (default line search)" => 8000.0,
    "Newton + BackTracking line search" => 500.0,
    "Trust-region" => 300.0,
)
const MAX_SCAN_PER_POINT_MEDIAN_MS = 1200.0
const MIN_CONVERGENCE_RATE = 0.90

function require_file(path::String)
    isfile(path) || error("Benchmark result file not found: $(path)")
end

function as_float(x)
    return x isa Number ? Float64(x) : parse(Float64, string(x))
end

function check_single_point()
    require_file(SINGLE_POINT_JSON)
    payload = JSON3.read(read(SINGLE_POINT_JSON, String))
    benches = payload.benchmarks
    isempty(benches) && error("single_point_benchmark.json has no benchmarks")

    violations = String[]
    for b in benches
        label = String(b.label)
        median_ms = as_float(b.stats.median)
        limit = get(SINGLE_POINT_MEDIAN_LIMITS, label, DEFAULT_SINGLE_POINT_MEDIAN_MS)
        if median_ms > limit
            push!(violations, "single-point $(label): median=$(median_ms)ms > $(limit)ms")
        end
    end
    return violations
end

function check_scan()
    require_file(SCAN_JSON)
    payload = JSON3.read(read(SCAN_JSON, String))
    benches = payload.benchmarks
    isempty(benches) && error("scan_benchmark.json has no benchmarks")

    violations = String[]
    for b in benches
        label = String(b.label)
        median_ms = as_float(b.stats.per_point_median_ms)
        conv = as_float(b.convergence_rate)
        if median_ms > MAX_SCAN_PER_POINT_MEDIAN_MS
            push!(violations, "scan $(label): per_point_median=$(median_ms)ms > $(MAX_SCAN_PER_POINT_MEDIAN_MS)ms")
        end
        if conv < MIN_CONVERGENCE_RATE
            push!(violations, "scan $(label): convergence=$(conv) < $(MIN_CONVERGENCE_RATE)")
        end
    end
    return violations
end

function main()
    violations = String[]
    append!(violations, check_single_point())
    append!(violations, check_scan())

    if !isempty(violations)
        println("[benchmark-thresholds] FAILED: $(length(violations)) violation(s)")
        for v in violations
            println(" - $(v)")
        end
        exit(1)
    end

    println("[benchmark-thresholds] OK: all thresholds satisfied")
    println("  DEFAULT_SINGLE_POINT_MEDIAN_MS=$(DEFAULT_SINGLE_POINT_MEDIAN_MS)")
    println("  MAX_SCAN_PER_POINT_MEDIAN_MS=$(MAX_SCAN_PER_POINT_MEDIAN_MS)")
    println("  MIN_CONVERGENCE_RATE=$(MIN_CONVERGENCE_RATE)")
end

main()
