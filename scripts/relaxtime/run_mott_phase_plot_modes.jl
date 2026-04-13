#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "plot_scan_csv.py")

const OBSERVABLES = ["M_pi", "M_K", "Gamma_pi", "Gamma_K", "M_u_plus_M_d", "M_u_plus_M_s"]
const MODE_AB_XIS = (-0.3, 0.0, 0.3)
const MODE_AB_GROUPS = [
    ("M_K", ["M_K", "Gamma_K", "M_u_plus_M_s"]),
    ("M_pi", ["M_pi", "Gamma_pi", "M_u_plus_M_d"]),
]

function _print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_mott_phase_plot_modes.jl --in <derived_csv> --out-dir <fig_dir>")
end

function _parse_args(args::Vector{String})
    input_csv = nothing
    out_dir = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end

        if arg == "--in"
            input_csv = require_value()
        elseif arg == "--out-dir"
            out_dir = require_value()
        elseif arg in ("-h", "--help")
            _print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    input_csv === nothing && throw(ArgumentError("--in is required"))
    out_dir === nothing && throw(ArgumentError("--out-dir is required"))
    return String(input_csv), String(out_dir)
end

@inline function _number_tag(v::Real)
    x = Float64(v)
    if isapprox(x, round(x); atol=1e-12, rtol=0.0)
        return string(Int(round(x)))
    end
    return replace(string(x), "." => "p")
end

function _collect_xi_values(csv_path::String)
    vals = Set{Float64}()
    header = String[]
    open(csv_path, "r") do io
        header_seen = false
        xi_idx = 0
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !header_seen
                header = [strip(x) for x in split(s, ',')]
                xi_idx = findfirst(==("xi"), header)
                xi_idx === nothing && throw(ArgumentError("input csv missing xi column"))
                header_seen = true
                continue
            end
            parts = split(s, ',')
            xi_idx > length(parts) && continue
            x = tryparse(Float64, strip(parts[xi_idx]))
            x === nothing && continue
            push!(vals, x)
        end
    end
    return sort(collect(vals))
end

function _run_python(args::Vector{String})
    try
        run(Cmd(vcat("python", args)))
        return
    catch
    end
    run(Cmd(vcat("python3", args)))
end

@inline function _is_mode_ab_xi(x::Real)
    xf = Float64(x)
    return any(t -> isapprox(xf, t; atol=1e-12, rtol=0.0), MODE_AB_XIS)
end

function _write_mode_ab_long_csv(input_csv::String, ys_list::Vector{String})
    out_csv = joinpath(mktempdir(), "mode_ab_long.csv")
    open(input_csv, "r") do src
        open(out_csv, "w") do dst
            println(dst, "T_MeV,series,value")
            header_seen = false
            t_idx = 0
            xi_idx = 0
            y_idxs = Int[]

            for line in eachline(src)
                s = strip(line)
                isempty(s) && continue
                startswith(s, "#") && continue

                if !header_seen
                    cols = [strip(x) for x in split(s, ',')]
                    t_idx = findfirst(==("T_MeV"), cols)
                    xi_idx = findfirst(==("xi"), cols)
                    (t_idx === nothing || xi_idx === nothing) && throw(ArgumentError("input csv missing T_MeV or xi column"))
                    for y in ys_list
                        yi = findfirst(==(y), cols)
                        yi === nothing && throw(ArgumentError("input csv missing $y column"))
                        push!(y_idxs, yi)
                    end
                    header_seen = true
                    continue
                end

                parts = split(s, ',')
                max_idx = maximum(vcat([t_idx, xi_idx], y_idxs))
                length(parts) < max_idx && continue

                t = tryparse(Float64, strip(parts[t_idx]))
                xi = tryparse(Float64, strip(parts[xi_idx]))
                (t === nothing || xi === nothing || !_is_mode_ab_xi(xi)) && continue

                xi_tag = _number_tag(xi)
                for (y, yi) in zip(ys_list, y_idxs)
                    yv = tryparse(Float64, strip(parts[yi]))
                    yv === nothing && continue
                    println(dst, string(t, ',', y, "__xi", xi_tag, ',', yv))
                end
            end
        end
    end
    return out_csv
end

function _mode_a(input_csv::String, out_root::String, xis::Vector{Float64})
    mode_a_dir = joinpath(out_root, "mode_a")
    mkpath(mode_a_dir)

    ys = join(OBSERVABLES, ',')
    tmp_out = mktempdir()
    args = String[
        PLOT_SCRIPT,
        "--mode", "lines",
        "--csv", input_csv,
        "--x", "T_MeV",
        "--ys", ys,
        "--split", "xi",
        "--multi-y",
        "--out-dir", tmp_out,
    ]
    _run_python(args)

    for xi in xis
        split_dir = joinpath(tmp_out, "xi=$(xi)")
        if !isdir(split_dir)
            continue
        end
        src = joinpath(split_dir, "multi_y_M_pi_M_K_Gamma_pi_Gamma_K_M_u_plus_M_d_M_u_plus_M_s_vs_T_MeV.png")
        if !isfile(src)
            files = readdir(split_dir; join=true)
            pngs = filter(f -> endswith(lowercase(f), ".png"), files)
            isempty(pngs) && continue
            src = pngs[1]
        end
        dst = joinpath(mode_a_dir, "mott_mode_a__xi$(_number_tag(xi)).png")
        cp(src, dst; force=true)
    end
end

function _mode_b(input_csv::String, out_root::String)
    mode_b_dir = joinpath(out_root, "mode_b")
    mkpath(mode_b_dir)

    ys = join(OBSERVABLES, ',')
    tmp_out = mktempdir()
    args = String[
        PLOT_SCRIPT,
        "--mode", "lines",
        "--csv", input_csv,
        "--x", "T_MeV",
        "--ys", ys,
        "--group", "xi",
        "--out-dir", tmp_out,
    ]
    _run_python(args)

    for obs in OBSERVABLES
        src = joinpath(tmp_out, "$(obs)_vs_T_MeV.png")
        isfile(src) || continue
        dst = joinpath(mode_b_dir, "mott_mode_b__$(obs).png")
        cp(src, dst; force=true)
    end
end

function _mode_ab(input_csv::String, out_root::String)
    mode_ab_dir = joinpath(out_root, "mode_ab")
    mkpath(mode_ab_dir)

    for (label, ys_list) in MODE_AB_GROUPS
        long_csv = _write_mode_ab_long_csv(input_csv, ys_list)
        tmp_out = mktempdir()
        args = String[
            PLOT_SCRIPT,
            "--mode", "lines",
            "--csv", long_csv,
            "--x", "T_MeV",
            "--ys", "value",
            "--group", "series",
            "--out-dir", tmp_out,
        ]
        _run_python(args)

        src = joinpath(tmp_out, "value_vs_T_MeV.png")
        if !isfile(src)
            pngs = filter(f -> endswith(lowercase(f), ".png"), readdir(tmp_out; join=true))
            isempty(pngs) && continue
            src = pngs[1]
        end

        dst = joinpath(mode_ab_dir, "mott_mode_ab__$(label)__xi3.png")
        cp(src, dst; force=true)
    end
end

function main()
    input_csv, out_dir = _parse_args(ARGS)
    isfile(PLOT_SCRIPT) || throw(ArgumentError("plot script not found: $PLOT_SCRIPT"))
    isfile(input_csv) || throw(ArgumentError("input csv not found: $input_csv"))

    mkpath(out_dir)
    xis = _collect_xi_values(input_csv)
    _mode_a(input_csv, out_dir, xis)
    _mode_b(input_csv, out_dir)
    _mode_ab(input_csv, out_dir)
    println("Wrote plot modes under: ", out_dir)
end

main()
