"""
沿 freeze-out 路径运行介子质量 workflow。

默认提供一组正式产出口径：

- xi = 0.0
- sqrt(s_NN) in [1, 200] GeV
- 在 log10(sqrt(s_NN)) 上等间距
- n = 30
- case 名：default_baseline_freezeout_xi0_loggrid_1to200_n30
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_freezeout_meson_mass_scan

struct ScanOptions
    output::String
    sqrt_s_values::Vector{Float64}
    xi::Float64
    freezeout_profile::String
    path_profile::String
    traversal::Symbol
    overwrite::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
    case_name::String
end

@inline function default_case_name()
    return "default_baseline_freezeout_xi0_loggrid_1to200_n30"
end

@inline function default_output_path(case_name::AbstractString)
    return joinpath(
        "data", "outputs", "results", "relaxtime", "meson_mass", "path_scan", "freezeout",
        string(case_name, ".csv"),
    )
end

function build_log10_sqrts_grid(min_GeV::Real, max_GeV::Real, n::Integer)
    lo = Float64(min_GeV)
    hi = Float64(max_GeV)
    count = Int(n)
    lo > 0.0 || throw(ArgumentError("min_GeV must be positive, got $(min_GeV)"))
    hi > lo || throw(ArgumentError("max_GeV must be greater than min_GeV"))
    count >= 2 || throw(ArgumentError("n must be >= 2, got $(n)"))
    xs = range(log10(lo), log10(hi); length=count)
    return [10.0^x for x in xs]
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_freezeout_meson_mass_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>                     输出 CSV")
    println("  --case-name <name>                  case 名；默认 default_baseline_freezeout_xi0_loggrid_1to200_n30")
    println("  --sqrts <value>                     追加一个 sqrt(s_NN) 点（GeV，可多次）")
    println("  --sqrts-list v1,v2,...              用逗号分隔的 sqrt(s_NN) 列表替换")
    println("  --loggrid-min <GeV>                 生成 log10 等间距网格下界")
    println("  --loggrid-max <GeV>                 生成 log10 等间距网格上界")
    println("  --loggrid-n <int>                   生成 log10 等间距网格点数")
    println("  --xi <value>                        各向异性 xi (default 0.0)")
    println("  --freezeout-profile <name>          freeze-out profile (default default)")
    println("  --path-profile <name>               path profile (default baseline_freezeout)")
    println("  --traversal <as_given|sqrts_ascending|sqrts_descending|muB_descending|muB_ascending>")
    println("  --overwrite                         覆盖输出文件")
    println("  --p-num <int>                       平衡态动量节点数 (default 24)")
    println("  --t-num <int>                       平衡态角度节点数 (default 8)")
    println("  --max-iter <int>                    solver/mass iterations 上限 (default 40)")
    println("  -h, --help                          显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :case_name => default_case_name(),
        :output => nothing,
        :sqrts_values => build_log10_sqrts_grid(1.0, 200.0, 30),
        :xi => 0.0,
        :freezeout_profile => "default",
        :path_profile => "baseline_freezeout",
        :traversal => :sqrts_ascending,
        :overwrite => false,
        :p_num => 24,
        :t_num => 8,
        :max_iter => 40,
        :loggrid_min => nothing,
        :loggrid_max => nothing,
        :loggrid_n => nothing,
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
        elseif arg == "--case-name"
            opts[:case_name] = require_value()
        elseif arg == "--sqrts"
            if opts[:sqrts_values] == build_log10_sqrts_grid(1.0, 200.0, 30)
                opts[:sqrts_values] = Float64[]
            end
            push!(opts[:sqrts_values], parse(Float64, require_value()))
        elseif arg == "--sqrts-list"
            raw = split(require_value(), ',')
            opts[:sqrts_values] = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
        elseif arg == "--loggrid-min"
            opts[:loggrid_min] = parse(Float64, require_value())
        elseif arg == "--loggrid-max"
            opts[:loggrid_max] = parse(Float64, require_value())
        elseif arg == "--loggrid-n"
            opts[:loggrid_n] = parse(Int, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--freezeout-profile"
            opts[:freezeout_profile] = require_value()
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

    if opts[:loggrid_min] !== nothing || opts[:loggrid_max] !== nothing || opts[:loggrid_n] !== nothing
        lo = something(opts[:loggrid_min], 1.0)
        hi = something(opts[:loggrid_max], 200.0)
        n = something(opts[:loggrid_n], 30)
        opts[:sqrts_values] = build_log10_sqrts_grid(lo, hi, n)
    end

    sqrt_s_vals = unique(Float64.(opts[:sqrts_values]))
    isempty(sqrt_s_vals) && error("sqrt(s_NN) grid must not be empty")

    output = opts[:output] === nothing ? default_output_path(String(opts[:case_name])) : String(opts[:output])
    return ScanOptions(
        output,
        sqrt_s_vals,
        Float64(opts[:xi]),
        String(opts[:freezeout_profile]),
        String(opts[:path_profile]),
        Symbol(opts[:traversal]),
        Bool(opts[:overwrite]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        String(opts[:case_name]),
    )
end

function main(args::Vector{String})
    opts = parse_args(args)
    mkpath(dirname(opts.output))

    result = run_freezeout_meson_mass_scan(
        sqrt_s_NN_values=opts.sqrt_s_values,
        xi=opts.xi,
        freezeout_profile_name=opts.freezeout_profile,
        path_profile_name=opts.path_profile,
        output_path=opts.output,
        traversal=opts.traversal,
        overwrite=opts.overwrite,
        p_num=opts.p_num,
        t_num=opts.t_num,
        max_iter=opts.max_iter,
    )

    println("freezeout meson-mass path scan completed")
    println("case_name=$(opts.case_name)")
    println("output_path=$(result.output_path)")
    println("points=$(result.points)")
    println("workflow_entry=$(result.workflow_entry)")
    return result
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
