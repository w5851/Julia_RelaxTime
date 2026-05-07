module PhaseGuidedTransportScanCLI

struct PhaseGuidedScanOptions
    mode::Symbol
    outdir::String
    case_name::String
    xi_values::Vector{Float64}
    muB_values::Vector{Float64}
    alpha_T_values::Vector{Float64}
    T_values::Vector{Float64}
    dry_run::Bool
    overwrite::Bool
    resume::Bool
end

function print_usage(io::IO=stdout)
    println(io, "Usage: julia --project=. scripts/relaxtime/run_phase_guided_transport_scan.jl [options]\n")
    println(io, "Options:")
    println(io, "  --mode <a|b>                 扫描模式：a=fixed muB + T/T_phase, b=fixed T + sparse muB")
    println(io, "  --outdir <dir>               输出目录（可选；默认写入 canonical results 路径）")
    println(io, "  --case-name <name>           canonical case 名称（default manual_case）")
    println(io, "  --xi-list v1,v2,...          连续 xi 采样列表（required）")
    println(io, "  --muB-list v1,v2,...         mode a/b 使用的离散 muB 列表（required）")
    println(io, "  --alphaT-list v1,v2,...      mode a 的 T/T_phase 倍率列表（default 1.0,1.1,1.2)")
    println(io, "  --T-list v1,v2,...           mode b 的固定温度列表（required for mode b）")
    println(io, "  --dry-run                    仅生成 sampling plan/README/provenance，不执行计算")
    println(io, "  --resume                     若 result CSV 已存在则跳过已有点")
    println(io, "  --overwrite                  覆盖已有输出")
    println(io, "  -h, --help                   显示帮助")
end

function _parse_float_list(raw::AbstractString)
    vals = Float64[]
    for item in split(String(raw), ',')
        s = strip(item)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    return unique(sort(vals))
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :mode => nothing,
        :outdir => nothing,
        :case_name => "manual_case",
        :xi_values => nothing,
        :muB_values => nothing,
        :alpha_T_values => Float64[1.0, 1.1, 1.2],
        :T_values => Float64[],
        :dry_run => false,
        :overwrite => false,
        :resume => false,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $(arg)")
            i += 1
            return args[i]
        end

        if arg == "--mode"
            mode_raw = lowercase(strip(require_value()))
            if mode_raw == "a"
                opts[:mode] = :mode_a_fixed_muB_phase_scaled
            elseif mode_raw == "b"
                opts[:mode] = :mode_b_fixed_T_sparse_muB
            else
                error("unsupported mode: $(mode_raw)")
            end
        elseif arg == "--outdir"
            opts[:outdir] = require_value()
        elseif arg == "--case-name"
            opts[:case_name] = strip(require_value())
        elseif arg == "--xi-list"
            opts[:xi_values] = _parse_float_list(require_value())
        elseif arg == "--muB-list"
            opts[:muB_values] = _parse_float_list(require_value())
        elseif arg == "--alphaT-list"
            opts[:alpha_T_values] = _parse_float_list(require_value())
        elseif arg == "--T-list"
            opts[:T_values] = _parse_float_list(require_value())
        elseif arg == "--dry-run"
            opts[:dry_run] = true
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--resume"
            opts[:resume] = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $(arg)")
        end
        i += 1
    end

    isnothing(opts[:mode]) && error("--mode is required")
    isnothing(opts[:xi_values]) && error("--xi-list is required")
    isnothing(opts[:muB_values]) && error("--muB-list is required")

    mode = opts[:mode]
    if mode == :mode_a_fixed_muB_phase_scaled
        isempty(opts[:alpha_T_values]) && error("mode a requires non-empty --alphaT-list")
    elseif mode == :mode_b_fixed_T_sparse_muB
        isempty(opts[:T_values]) && error("mode b requires non-empty --T-list")
    end

    mode_dir = mode == :mode_a_fixed_muB_phase_scaled ? "mode_a_fixed_muB_phase_scaled" : "mode_b_fixed_T_sparse_muB"
    case_name = isempty(strip(String(opts[:case_name]))) ? "manual_case" : String(opts[:case_name])
    outdir = isnothing(opts[:outdir]) ?
        joinpath("data", "outputs", "results", "relaxtime", "transport", "phase_guided", mode_dir, case_name) :
        String(opts[:outdir])

    return PhaseGuidedScanOptions(
        mode,
        outdir,
        case_name,
        Float64.(opts[:xi_values]),
        Float64.(opts[:muB_values]),
        Float64.(opts[:alpha_T_values]),
        Float64.(opts[:T_values]),
        Bool(opts[:dry_run]),
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
    )
end

export PhaseGuidedScanOptions
export parse_args
export print_usage

end # module PhaseGuidedTransportScanCLI
