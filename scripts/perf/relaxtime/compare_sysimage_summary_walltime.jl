#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using JSON3

if length(ARGS) < 2
    error("usage: julia --project=. scripts/perf/relaxtime/compare_sysimage_summary_walltime.jl <summary-a.json> <summary-b.json>")
end

const SUMMARY_A = normpath(ARGS[1])
const SUMMARY_B = normpath(ARGS[2])

function _load_summary(path::String)
    isfile(path) || error("missing summary file: $(path)")
    return JSON3.read(read(path, String))
end

@inline _ms(x) = Float64(x)
@inline _pct(delta, base) = iszero(base) ? NaN : 100.0 * delta / base

function _print_metric(label::String, a::Float64, b::Float64)
    delta = b - a
    pct = _pct(delta, a)
    sign = delta > 0 ? "+" : ""
    println("  $(label): $(round(a; digits=1)) -> $(round(b; digits=1)) ($(sign)$(round(delta; digits=1)), $(sign)$(round(pct; digits=2))%)")
end

a = _load_summary(SUMMARY_A)
b = _load_summary(SUMMARY_B)

println("[sysimage-summary-walltime]")
println("  A=$(SUMMARY_A)")
println("  B=$(SUMMARY_B)")
println("  workload_a=$(get(a, :workload, missing))")
println("  workload_b=$(get(b, :workload, missing))")
println("  label_a=$(get(a, :label, missing))")
println("  label_b=$(get(b, :label, missing))")

println("\n[wall_ms]")
_print_metric("no-sys", _ms(a.no_sys.wall_ms), _ms(b.no_sys.wall_ms))
_print_metric("with-sys", _ms(a.with_sys.wall_ms), _ms(b.with_sys.wall_ms))
_print_metric("improvement", _ms(a.improvement.wall_ms), _ms(b.improvement.wall_ms))

println("\n[trace_lines]")
_print_metric("no-sys", _ms(a.no_sys.trace_lines), _ms(b.no_sys.trace_lines))
_print_metric("with-sys", _ms(a.with_sys.trace_lines), _ms(b.with_sys.trace_lines))
_print_metric("improvement", _ms(a.improvement.trace_lines), _ms(b.improvement.trace_lines))

println("\n[focus_trace_lines]")
_print_metric("no-sys", _ms(a.no_sys.focus_trace_lines), _ms(b.no_sys.focus_trace_lines))
_print_metric("with-sys", _ms(a.with_sys.focus_trace_lines), _ms(b.with_sys.focus_trace_lines))
_print_metric("improvement", _ms(a.improvement.focus_trace_lines), _ms(b.improvement.focus_trace_lines))
