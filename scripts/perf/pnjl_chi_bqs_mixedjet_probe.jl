using Printf
using Statistics

const SCRIPT_START = time()
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS_PATH = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")

if !isdefined(Main, :Models)
    Base.include(Main, MODELS_PATH)
end

using .Models

const PNJL = Models.pnjl_module()

function _arg_value(name::String, default::String)
    prefix = "--" * name * "="
    for arg in ARGS
        startswith(arg, prefix) && return arg[length(prefix) + 1:end]
    end
    return default
end

function _parse_orders(raw::String)
    cases = NTuple{3, Int}[]
    for part in split(raw, ';'; keepempty=false)
        cols = split(strip(part), ','; keepempty=false)
        length(cols) == 3 || throw(ArgumentError("orders entries must look like nB,nQ,nS; got $(part)"))
        orders = (parse(Int, strip(cols[1])), parse(Int, strip(cols[2])), parse(Int, strip(cols[3])))
        all(n -> n >= 0, orders) || throw(ArgumentError("orders must be non-negative, got $(orders)"))
        count(!iszero, orders) >= 2 || throw(ArgumentError("mixedjet probe expects mixed orders with at least two active axes, got $(orders)"))
        push!(cases, orders)
    end
    isempty(cases) && throw(ArgumentError("at least one orders case is required"))
    return cases
end

function _rss_bytes()
    try
        return Sys.maxrss()
    catch
        return -1
    end
end

@inline function _chi_case(T_fm::Float64, muB_fm::Float64, muQ_fm::Float64, muS_fm::Float64, orders::NTuple{3, Int}; kwargs...)
    return PNJL.chi_BQS(
        T_fm,
        muB_fm,
        muQ_fm,
        muS_fm;
        orders=orders,
        derivative_backend=:mixedjet,
        kwargs...,
    )
end

function _run_case(
    orders::NTuple{3, Int},
    T_fm::Float64,
    muB_fm::Float64,
    muQ_fm::Float64,
    muS_fm::Float64,
    next_delta_T_fm::Float64,
    repeats::Int;
    kwargs...,
)
    GC.gc()
    first = @timed _chi_case(T_fm, muB_fm, muQ_fm, muS_fm, orders; kwargs...)

    GC.gc()
    next_point = @timed _chi_case(T_fm + next_delta_T_fm, muB_fm, muQ_fm, muS_fm, orders; kwargs...)

    warm_times = Float64[]
    warm_allocs = Int[]
    for _ in 1:repeats
        GC.gc()
        t = @timed _chi_case(T_fm, muB_fm, muQ_fm, muS_fm, orders; kwargs...)
        push!(warm_times, t.time)
        push!(warm_allocs, t.bytes)
    end

    return (
        orders=orders,
        total_order=sum(orders),
        active_axes=count(!iszero, orders),
        first_call_s=first.time,
        first_alloc_bytes=first.bytes,
        same_process_next_point_s=next_point.time,
        same_process_next_point_alloc_bytes=next_point.bytes,
        warm_median_s=median(warm_times),
        warm_min_s=minimum(warm_times),
        warm_alloc_median_bytes=median(warm_allocs),
        value=first.value,
        peak_rss_bytes=_rss_bytes(),
    )
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "nB,nQ,nS,total_order,active_axes,first_call_s,first_alloc_bytes,same_process_next_point_s,same_process_next_point_alloc_bytes,warm_median_s,warm_min_s,warm_alloc_median_bytes,peak_rss_bytes,value")
        for row in rows
            @printf(
                io,
                "%d,%d,%d,%d,%d,%.9g,%d,%.9g,%d,%.9g,%.9g,%.0f,%d,%.16e\n",
                row.orders[1],
                row.orders[2],
                row.orders[3],
                row.total_order,
                row.active_axes,
                row.first_call_s,
                row.first_alloc_bytes,
                row.same_process_next_point_s,
                row.same_process_next_point_alloc_bytes,
                row.warm_median_s,
                row.warm_min_s,
                row.warm_alloc_median_bytes,
                row.peak_rss_bytes,
                row.value,
            )
        end
    end
end

orders_cases = _parse_orders(_arg_value("orders", "1,1,0;2,1,0;1,1,1;2,1,1"))
T_fm = parse(Float64, _arg_value("T_fm", "0.57"))
muB_fm = parse(Float64, _arg_value("muB_fm", "0.18"))
muQ_fm = parse(Float64, _arg_value("muQ_fm", "0.05"))
muS_fm = parse(Float64, _arg_value("muS_fm", "0.02"))
next_delta_T_fm = parse(Float64, _arg_value("next_delta_T_fm", "0.0025"))
p_num = parse(Int, _arg_value("p_num", "8"))
t_num = parse(Int, _arg_value("t_num", "4"))
repeats = parse(Int, _arg_value("repeats", "5"))
linear_solve = Symbol(_arg_value("linear_solve", "auto"))
series_residual_tol = parse(Float64, _arg_value("series_residual_tol", "1e-7"))
output = strip(_arg_value("output", ""))

repeats >= 1 || throw(ArgumentError("repeats must be >= 1, got $(repeats)"))

kwargs = (;
    xi=0.0,
    p_num=p_num,
    t_num=t_num,
    linear_solve=linear_solve,
    series_residual_tol=series_residual_tol,
)

rows = [
    _run_case(orders, T_fm, muB_fm, muQ_fm, muS_fm, next_delta_T_fm, repeats; kwargs...)
    for orders in orders_cases
]

println("PNJL chi_BQS mixed Taylor jet performance probe")
println("T_fm=$(T_fm)")
println("muB_fm=$(muB_fm)")
println("muQ_fm=$(muQ_fm)")
println("muS_fm=$(muS_fm)")
println("next_delta_T_fm=$(next_delta_T_fm)")
println("p_num=$(p_num)")
println("t_num=$(t_num)")
println("repeats=$(repeats)")
println("linear_solve=$(linear_solve)")
println("series_residual_tol=$(series_residual_tol)")

for row in rows
    println("case_orders=$(row.orders)")
    println("  total_order=$(row.total_order)")
    println("  active_axes=$(row.active_axes)")
    println("  first_call_s=$(row.first_call_s)")
    println("  first_alloc_bytes=$(row.first_alloc_bytes)")
    println("  same_process_next_point_s=$(row.same_process_next_point_s)")
    println("  same_process_next_point_alloc_bytes=$(row.same_process_next_point_alloc_bytes)")
    println("  warm_median_s=$(row.warm_median_s)")
    println("  warm_min_s=$(row.warm_min_s)")
    println("  warm_alloc_median_bytes=$(row.warm_alloc_median_bytes)")
    println("  peak_rss_bytes=$(row.peak_rss_bytes)")
    println("  value=$(row.value)")
end

println("cold_wall_s=$(time() - SCRIPT_START)")

if !isempty(output)
    _write_csv(output, rows)
    println("csv_output=$(abspath(output))")
end
