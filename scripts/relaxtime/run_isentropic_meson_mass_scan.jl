"""
沿 fixed-σ 等熵线运行介子质量 workflow。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_isentropic_meson_mass_scan

struct ScanOptions
    output::String
    T_values::Vector{Float64}
    sigma_target::Float64
    xi::Float64
    path_profile::String
    traversal::Symbol
    overwrite::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_isentropic_meson_mass_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>                     输出 CSV")
    println("  --sigma-target <value>              目标 sigma = s/n_B")
    println("  --T <MeV>                           追加一个温度点（可多次）")
    println("  --T-list v1,v2,...                  用逗号分隔的温度点列表替换")
    println("  --xi <value>                        各向异性 xi (default 0.0)")
    println("  --path-profile <name>               path profile (default default)")
    println("  --traversal <as_given|T_ascending|T_descending>")
    println("  --overwrite                         覆盖输出文件")
    println("  --p-num <int>                       平衡态动量节点数 (default 24)")
    println("  --t-num <int>                       平衡态角度节点数 (default 8)")
    println("  --max-iter <int>                    solver/mass iterations 上限 (default 40)")
    println("  -h, --help                          显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :output => joinpath("data", "outputs", "results", "relaxtime", "meson_mass", "path_scan", "isentropic", "isentropic_meson_mass_path_scan.csv"),
        :T_values => Float64[140.0, 160.0, 180.0],
        :sigma_target => nothing,
        :xi => 0.0,
        :path_profile => "default",
        :traversal => :T_ascending,
        :overwrite => false,
        :p_num => 24,
        :t_num => 8,
        :max_iter => 40,
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

        if arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--sigma-target"
            opts[:sigma_target] = parse(Float64, require_value())
        elseif arg == "--T"
            if opts[:T_values] == Float64[140.0, 160.0, 180.0]
                opts[:T_values] = Float64[]
            end
            push!(opts[:T_values], parse(Float64, require_value()))
        elseif arg == "--T-list"
            raw = split(require_value(), ',')
            opts[:T_values] = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--path-profile"
            opts[:path_profile] = require_value()
        elseif arg == "--traversal"
            opts[:traversal] = Symbol(require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    opts[:sigma_target] === nothing && error("missing required option: --sigma-target")
    T_values = unique(Float64.(opts[:T_values]))
    isempty(T_values) && error("T grid must not be empty")

    return ScanOptions(
        String(opts[:output]),
        T_values,
        Float64(opts[:sigma_target]),
        Float64(opts[:xi]),
        String(opts[:path_profile]),
        Symbol(opts[:traversal]),
        Bool(opts[:overwrite]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
    )
end

function main(args::Vector{String})
    opts = parse_args(args)
    mkpath(dirname(opts.output))

    result = run_isentropic_meson_mass_scan(
        T_MeV_values=opts.T_values,
        sigma_target=opts.sigma_target,
        xi=opts.xi,
        path_profile_name=opts.path_profile,
        output_path=opts.output,
        traversal=opts.traversal,
        overwrite=opts.overwrite,
        p_num=opts.p_num,
        t_num=opts.t_num,
        max_iter=opts.max_iter,
    )

    println("isentropic meson-mass path scan completed")
    println("output_path=$(result.output_path)")
    println("points=$(result.points)")
    println("workflow_entry=$(result.workflow_entry)")
    return result
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
