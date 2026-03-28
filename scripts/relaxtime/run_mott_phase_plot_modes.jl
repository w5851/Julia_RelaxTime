#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "plot_scan_csv.py")

const OBSERVABLES = ["M_pi", "M_K", "Gamma_pi", "Gamma_K", "M_u_plus_M_d", "M_u_plus_M_s"]

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

function main()
    input_csv, out_dir = _parse_args(ARGS)
    isfile(PLOT_SCRIPT) || throw(ArgumentError("plot script not found: $PLOT_SCRIPT"))
    isfile(input_csv) || throw(ArgumentError("input csv not found: $input_csv"))

    mkpath(out_dir)
    xis = _collect_xi_values(input_csv)
    _mode_a(input_csv, out_dir, xis)
    _mode_b(input_csv, out_dir)
    println("Wrote plot modes under: ", out_dir)
end

main()
