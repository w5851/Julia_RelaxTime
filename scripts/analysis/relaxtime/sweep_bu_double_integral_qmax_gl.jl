"""
Phase E3: GL + 硬截断下的 q_max 外推诊断。

目标：
- 固定 omega_max 与节点数
- 扫描 q_max
- 观察外层 q 壳贡献与双积分估计是否趋于稳定
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const SCRIPT_PATH = joinpath(@__DIR__, "diagnose_bu_double_integral_minimal.jl")

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "bu_phase_shift_minimal",
    "bu_double_integral_qmax_sweep.csv",
)

const DEFAULT_QMAX_VALUES = [4.0, 6.0, 8.0, 10.0]

function _parse_float_list(raw::AbstractString)
    vals = Float64[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("list cannot be empty"))
    return vals
end

function _selected_qmax_values()
    raw = strip(get(ENV, "PHASE_E_QMAX_SWEEP", ""))
    isempty(raw) && return DEFAULT_QMAX_VALUES
    return _parse_float_list(raw)
end

function _run_one(qmax::Float64, tmp_output::String)
    old_qmax = get(ENV, "PHASE_E_Q_MAX", nothing)
    old_output = get(ENV, "PHASE_E_OUTPUT", nothing)
    try
        ENV["PHASE_E_Q_MAX"] = string(qmax)
        ENV["PHASE_E_OUTPUT"] = tmp_output
        Base.include(Main, SCRIPT_PATH)
        Main.main()
    finally
        if old_qmax === nothing
            pop!(ENV, "PHASE_E_Q_MAX"; default=nothing)
        else
            ENV["PHASE_E_Q_MAX"] = old_qmax
        end
        if old_output === nothing
            pop!(ENV, "PHASE_E_OUTPUT"; default=nothing)
        else
            ENV["PHASE_E_OUTPUT"] = old_output
        end
    end
end

function _read_rows(path::String)
    rows = NamedTuple[]
    open(path, "r") do io
        header = split(strip(readline(io)), ',')
        idx = Dict(name => i for (i, name) in enumerate(header))
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            cols = split(s, ',')
            push!(rows, (
                T_MeV=parse(Float64, cols[idx["T_MeV"]]),
                meson=cols[idx["meson"]],
                q_max=parse(Float64, cols[idx["qmax_request"]]),
                omega_max=parse(Float64, cols[idx["omega_max_request"]]),
                q_nodes=parse(Int, cols[idx["q_nodes"]]),
                omega_nodes=parse(Int, cols[idx["omega_nodes"]]),
                q_integral_estimate=parse(Float64, cols[idx["q_integral_estimate"]]),
                omega_shell_at_qmax=parse(Float64, cols[idx["omega_shell_at_qmax"]]),
                double_integral_estimate=parse(Float64, cols[idx["double_integral_estimate"]]),
            ))
        end
    end
    return rows
end

@inline function _fmt(x)
    return string(x)
end

function main()
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
    mkpath(dirname(output))

    qmax_values = _selected_qmax_values()
    tmp_output = joinpath(dirname(output), "_tmp_bu_double_integral_qmax.csv")
    all_rows = NamedTuple[]

    for qmax in qmax_values
        _run_one(qmax, tmp_output)
        append!(all_rows, _read_rows(tmp_output))
    end

    rm(tmp_output; force=true)

    open(output, "w") do io
        println(io, "qmax_request,T_MeV,meson,q_max,omega_max,q_nodes,omega_nodes,q_integral_estimate,omega_shell_at_qmax,double_integral_estimate")
        for r in sort(all_rows; by=x -> (x.q_max, x.T_MeV, x.meson))
            println(io, join((
                _fmt(r.q_max),
                _fmt(r.T_MeV),
                r.meson,
                _fmt(r.q_max),
                _fmt(r.omega_max),
                _fmt(r.q_nodes),
                _fmt(r.omega_nodes),
                _fmt(r.q_integral_estimate),
                _fmt(r.omega_shell_at_qmax),
                _fmt(r.double_integral_estimate),
            ), ','))
        end
    end

    println("Wrote qmax sweep CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
