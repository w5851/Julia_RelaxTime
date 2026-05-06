"""
显式外部离散路径点列上的介子数密度 workflow 扫描。

特点：
- 统一通过 `Models.run_external_path_meson_density_scan` 进入；
- 路径输入来自 CSV，而不是内部 crossover / freezeout 生成器；
- 适合 Friesen 2019 phase-line 一类外部 `(T, mu_B)` 点列复现。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_external_path_meson_density_scan

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_external_path_meson_density_scan.jl [options]\n")
    println("Options:")
    println("  --input-csv <path>                    external path CSV")
    println("  --output <path>                       output CSV")
    println("  --source-fig <name>                   filter source_fig")
    println("  --case-id <name>                      filter case_id")
    println("  --line-style <solid|dashed>           filter line_style")
    println("  --flavor-chemical-profile <name>      flavor chemical profile")
    println("  --meson-chemical-profile <name>       meson chemical profile")
    println("  --regime <stable|phase_shift_current|phase_shift_gbu|strict_bw_stage1|strict_bw_stage2>")
    println("  --xi <float>                          anisotropy xi")
    println("  --overwrite                           overwrite output path")
    println("  --p-num <int>                         equilibrium momentum nodes")
    println("  --t-num <int>                         equilibrium theta nodes")
    println("  --max-iter <int>                      solver / mass iterations")
    println("  --qmax <float>                        q upper bound")
    println("  --q-nodes <int>                       q nodes")
    println("  --omega-min <float>                   phase-shift omega lower bound")
    println("  --omega-max <float>                   omega upper bound")
    println("  --omega-nodes <int>                   omega nodes")
    println("  --eta <float>                         phase-shift eta")
    println("  -h, --help                            show help")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :input_csv => raw"D:\Desktop\paper\dev\outputs\formalized\friesen2019_phase_lines\friesen2019_phase_lines.csv",
        :output => Models.ExternalPathMesonDensityScan.DEFAULT_EXTERNAL_PATH_MESON_DENSITY_OUTPUT_PATH,
        :source_fig => nothing,
        :case_id => nothing,
        :line_style => nothing,
        :flavor_chemical_profile => "friesen2019_mu_s_0p55",
        :meson_chemical_profile => "friesen2019_kplus_over_piplus_neutral",
        :regime => :stable,
        :xi => 0.0,
        :overwrite => false,
        :p_num => 24,
        :t_num => 8,
        :max_iter => 40,
        :qmax => 12.0,
        :q_nodes => 48,
        :omega_min => 0.05,
        :omega_max => 10.0,
        :omega_nodes => 48,
        :eta => 1e-6,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            val = args[i + 1]
            i += 1
            return val
        end
        if arg == "--input-csv"
            opts[:input_csv] = require_value()
        elseif arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--source-fig"
            opts[:source_fig] = require_value()
        elseif arg == "--case-id"
            opts[:case_id] = require_value()
        elseif arg == "--line-style"
            opts[:line_style] = require_value()
        elseif arg == "--flavor-chemical-profile"
            opts[:flavor_chemical_profile] = require_value()
        elseif arg == "--meson-chemical-profile"
            opts[:meson_chemical_profile] = require_value()
        elseif arg == "--regime"
            opts[:regime] = Symbol(require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--qmax"
            opts[:qmax] = parse(Float64, require_value())
        elseif arg == "--q-nodes"
            opts[:q_nodes] = parse(Int, require_value())
        elseif arg == "--omega-min"
            opts[:omega_min] = parse(Float64, require_value())
        elseif arg == "--omega-max"
            opts[:omega_max] = parse(Float64, require_value())
        elseif arg == "--omega-nodes"
            opts[:omega_nodes] = parse(Int, require_value())
        elseif arg == "--eta"
            opts[:eta] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return opts
end

function read_points(path::String; source_fig=nothing, case_id=nothing, line_style=nothing)
    lines = readlines(path)
    isempty(lines) && error("input CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    idx = Dict(name => i for (i, name) in enumerate(header))
    required = ("source_fig", "case_id", "line_style", "point_index", "muB_MeV", "T_MeV")
    for key in required
        haskey(idx, key) || error("input CSV missing column: $key")
    end
    points = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        row_source = strip(cols[idx["source_fig"]])
        row_case = strip(cols[idx["case_id"]])
        row_style = strip(cols[idx["line_style"]])
        source_fig !== nothing && row_source != source_fig && continue
        case_id !== nothing && row_case != case_id && continue
        line_style !== nothing && row_style != line_style && continue
        push!(points, (
            source_fig=row_source,
            case_id=row_case,
            line_style=row_style,
            point_index=parse(Int, cols[idx["point_index"]]),
            muB_MeV=parse(Float64, cols[idx["muB_MeV"]]),
            T_MeV=parse(Float64, cols[idx["T_MeV"]]),
        ))
    end
    isempty(points) && error("no rows selected from input CSV")
    return points
end

function main(args::Vector{String})
    opts = parse_args(args)
    points = read_points(
        String(opts[:input_csv]);
        source_fig=opts[:source_fig],
        case_id=opts[:case_id],
        line_style=opts[:line_style],
    )
    run_external_path_meson_density_scan(
        points=points,
        xi=Float64(opts[:xi]),
        flavor_chemical_profile_name=String(opts[:flavor_chemical_profile]),
        meson_chemical_profile_name=String(opts[:meson_chemical_profile]),
        regime=Symbol(opts[:regime]),
        output_path=String(opts[:output]),
        overwrite=Bool(opts[:overwrite]),
        p_num=Int(opts[:p_num]),
        t_num=Int(opts[:t_num]),
        max_iter=Int(opts[:max_iter]),
        strict_bw_qmax=Float64(opts[:qmax]),
        strict_bw_q_nodes=Int(opts[:q_nodes]),
        strict_bw_omega_max=Float64(opts[:omega_max]),
        strict_bw_omega_nodes=Int(opts[:omega_nodes]),
        phase_shift_qmax=Float64(opts[:qmax]),
        phase_shift_q_nodes=Int(opts[:q_nodes]),
        phase_shift_omega_min=Float64(opts[:omega_min]),
        phase_shift_omega_max=Float64(opts[:omega_max]),
        phase_shift_omega_nodes=Int(opts[:omega_nodes]),
        phase_shift_eta=Float64(opts[:eta]),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
