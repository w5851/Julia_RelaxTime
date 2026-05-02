"""
Phase E3 前的最小 BU omega 积分试探诊断。

目标：
- 复用 Phase E 相位诊断输出
- 先固定 q = 0
- 评估 omega 积分核在当前窗口下是否可积、是否趋于收敛

当前使用分部积分口径的被积函数诊断量：
    I(omega) = g(omega) * (1 + g(omega)) * delta_M(omega)

注意：
- 这里只是 omega 方向的试探诊断
- 不是完整 BU 数密度
- 不包含 q 积分
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const DEFAULT_DETAIL = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "bu_phase_shift_minimal",
    "bu_phase_shift_detail.csv",
)

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "bu_phase_shift_minimal",
    "bu_omega_integrand_minimal.csv",
)

function _parse_float(x::AbstractString)
    return parse(Float64, x)
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
                T_MeV=_parse_float(cols[idx["T_MeV"]]),
                meson=Symbol(cols[idx["meson"]]),
                q=_parse_float(cols[idx["q"]]),
                omega=_parse_float(cols[idx["omega"]]),
                bose_g=_parse_float(cols[idx["bose_g"]]),
                phase_unwrapped=_parse_float(cols[idx["phase_unwrapped"]]),
            ))
        end
    end
    return rows
end

@inline function _fmt(x)
    return string(x)
end

function _group_key(r)
    return (r.T_MeV, r.meson, r.q)
end

function _trapz(xs::Vector{Float64}, ys::Vector{Float64})
    length(xs) == length(ys) || throw(ArgumentError("xs and ys must have same length"))
    length(xs) >= 2 || return 0.0
    acc = 0.0
    for i in 1:length(xs)-1
        dx = xs[i + 1] - xs[i]
        acc += 0.5 * dx * (ys[i] + ys[i + 1])
    end
    return acc
end

function main()
    input = length(ARGS) >= 1 ? abspath(ARGS[1]) : DEFAULT_DETAIL
    output = length(ARGS) >= 2 ? abspath(ARGS[2]) : DEFAULT_OUTPUT

    rows = _read_rows(input)
    rows = filter(r -> iszero(r.q), rows)
    groups = Dict{Tuple{Float64,Symbol,Float64},Vector{NamedTuple}}()
    for r in rows
        push!(get!(groups, _group_key(r), NamedTuple[]), r)
    end

    mkpath(dirname(output))
    open(output, "w") do io
        println(io, "T_MeV,meson,q,omega_max,omega_integral_estimate,integrand_tail")
        for (key, group_rows) in sort(collect(groups); by=x -> x[1])
            sorted_rows = sort(group_rows; by=r -> r.omega)
            xs = Float64[r.omega for r in sorted_rows]
            ys = Float64[r.bose_g * (1.0 + r.bose_g) * r.phase_unwrapped for r in sorted_rows]
            total = _trapz(xs, ys)
            tail = ys[end]
            println(io, join((
                _fmt(key[1]), String(key[2]), _fmt(key[3]), _fmt(xs[end]), _fmt(total), _fmt(tail)
            ), ','))
        end
    end

    println("Wrote omega-integrand diagnostic CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
