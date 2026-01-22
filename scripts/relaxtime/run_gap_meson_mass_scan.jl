"""
批量扫描 (T, μ_B, ξ) 网格，串联：PNJL 平衡求解 → 介子质量/宽度（MesonMass）→ Mott 阈值与 gap。

输出：CSV（默认 data/outputs/results/relaxtime/gap_meson_mass_scan.csv）

单位约定：
- CLI 参数的温度/化学势使用 MeV（更符合扫描习惯）；脚本内部换算到 fm⁻¹。
- PNJL.solve 使用“夸克化学势” μ_q（fm⁻¹），因此 μ_q = μ_B/3。

示例：
  julia --project=. scripts/relaxtime/run_gap_meson_mass_scan.jl --tmin 120 --tmax 220 --tstep 10 --mubmin 0 --mubmax 800 --mubstep 200 --xi-list 0.0 --overwrite

说明：
- 本脚本输出每个点的阈值与 gap；Mott 点可通过 gap 过零在后处理中定位。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "MesonMassWorkflow.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .MesonMassWorkflow: solve_gap_and_meson_point, DEFAULT_MESONS
using .ScanCSV: ScanCSV

struct ScanOptions
    output::String
    xi_values::Vector{Float64}
    tmin_mev::Float64
    tmax_mev::Float64
    tstep_mev::Float64
    mubmin_mev::Float64
    mubmax_mev::Float64
    mubstep_mev::Float64
    overwrite::Bool
    resume::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_gap_meson_mass_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>             输出 CSV (default data/outputs/results/relaxtime/gap_meson_mass_scan.csv)")
    println("  --xi <value>                追加一个 ξ 值（可多次传入）")
    println("  --xi-list v1,v2,...         用逗号分隔的 ξ 列表替换")
    println("  --tmin/--tmax/--tstep <MeV> 温度范围与步长")
    println("  --mubmin/--mubmax/--mubstep <MeV> 重子化学势 μ_B 范围与步长")
    println("  --overwrite                 覆盖输出文件")
    println("  --no-resume                 禁用跳过逻辑，强制重算")
    println("  --p-num <int>               能隙求解动量节点数 (default 12)")
    println("  --t-num <int>               能隙求解角度节点数 (default 6)")
    println("  --max-iter <int>            NLsolve iterations 上限 (default 40)")
    println("  -h, --help                  显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :output => joinpath("data", "outputs", "results", "relaxtime", "gap_meson_mass_scan.csv"),
        :xi_values => Float64[0.0],
        :tmin => 120.0,
        :tmax => 220.0,
        :tstep => 10.0,
        :mubmin => 0.0,
        :mubmax => 800.0,
        :mubstep => 200.0,
        :overwrite => false,
        :resume => true,
        :p_num => 12,
        :t_num => 6,
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
        elseif arg == "--xi"
            val = parse(Float64, require_value())
            if opts[:xi_values] == Float64[0.0]
                opts[:xi_values] = Float64[]
            end
            push!(opts[:xi_values], val)
        elseif arg == "--xi-list"
            raw = split(require_value(), ',')
            vals = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
            opts[:xi_values] = vals
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--tstep"
            opts[:tstep] = parse(Float64, require_value())
        elseif arg == "--mubmin"
            opts[:mubmin] = parse(Float64, require_value())
        elseif arg == "--mubmax"
            opts[:mubmax] = parse(Float64, require_value())
        elseif arg == "--mubstep"
            opts[:mubstep] = parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-resume"
            opts[:resume] = false
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    xi_vals = unique(sort(Float64.(opts[:xi_values])))
    isempty(xi_vals) && (xi_vals = Float64[0.0])

    tstep = Float64(opts[:tstep]); tstep > 0 || error("tstep must be positive")
    mubstep = Float64(opts[:mubstep]); mubstep > 0 || error("mubstep must be positive")

    return ScanOptions(
        String(opts[:output]),
        xi_vals,
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        tstep,
        Float64(opts[:mubmin]),
        Float64(opts[:mubmax]),
        mubstep,
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
    )
end

@inline function _format(x)
    x isa Bool && return (x ? "true" : "false")
    return string(x)
end

function _header_cols(mesons::Tuple{Vararg{Symbol}})
    cols = String[
        "T_MeV", "muB_MeV", "T_fm", "muB_fm", "mu_fm", "xi", "Phi", "Phibar", "m_u", "m_d", "m_s",
    ]

    for m in mesons
        push!(cols, "M_" * String(m))
        push!(cols, "Gamma_" * String(m))
        push!(cols, "converged_" * String(m))
        push!(cols, "residual_" * String(m))

        if m in (:eta, :eta_prime, :sigma, :sigma_prime)
            push!(cols, "thr_uu_" * String(m))
            push!(cols, "thr_ss_" * String(m))
            push!(cols, "thr_min_" * String(m))
            push!(cols, "gap_uu_" * String(m))
            push!(cols, "gap_ss_" * String(m))
            push!(cols, "gap_min_" * String(m))
        else
            push!(cols, "threshold_" * String(m))
            push!(cols, "gap_" * String(m))
        end
    end

    return cols
end

function _row_to_values(cols::Vector{String}, row::Dict{String,Any})
    return [_format(get(row, c, "")) for c in cols]
end

function main()
    opts = parse_args(ARGS)

    outpath = opts.output
    mkpath(dirname(outpath))

    mesons = DEFAULT_MESONS
    cols = _header_cols(mesons)

    key_cols = ["T_MeV", "muB_MeV", "xi"]
    existing = (isfile(outpath) && opts.resume && !opts.overwrite) ? ScanCSV.read_existing_keys(outpath, key_cols) : Set{Tuple{Vararg{Float64}}}()

    if opts.overwrite && isfile(outpath)
        rm(outpath)
        empty!(existing)
    end

    is_new = !isfile(outpath)

    open(outpath, is_new ? "w" : "a") do io
        if is_new
            ScanCSV.write_metadata(io, Dict(
                "format" => "scan_csv_v1",
                "script" => "scripts/relaxtime/run_gap_meson_mass_scan.jl",
                "mesons" => join(String.(mesons), ","),
                "note" => "mu_fm denotes quark chemical potential (muB/3)"
            ))
            ScanCSV.write_header(io, cols)
        end

        for xi in opts.xi_values
            T = opts.tmin_mev
            while T <= opts.tmax_mev + 1e-9
                muB = opts.mubmin_mev
                while muB <= opts.mubmax_mev + 1e-9
                    key = (Float64(T), Float64(muB), Float64(xi))
                    if key in existing
                        muB += opts.mubstep_mev
                        continue
                    end

                    T_fm = T / ħc_MeV_fm
                    muB_fm = muB / ħc_MeV_fm
                    mu_fm = muB_fm / 3.0

                    # 平衡求解 + 介子质量
                    res = solve_gap_and_meson_point(T_fm, mu_fm;
                        xi=xi,
                        mesons=mesons,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        solver_kwargs=(; iterations=opts.max_iter),
                    )

                    qp = res.quark_params
                    tp = res.thermo_params

                    row = Dict{String,Any}(
                        "T_MeV" => T,
                        "muB_MeV" => muB,
                        "T_fm" => T_fm,
                        "muB_fm" => muB_fm,
                        "mu_fm" => mu_fm,
                        "xi" => xi,
                        "Phi" => tp.Φ,
                        "Phibar" => tp.Φbar,
                        "m_u" => qp.m.u,
                        "m_d" => qp.m.d,
                        "m_s" => qp.m.s,
                    )

                    for m in mesons
                        mr = res.meson_results[m]
                        row["M_" * String(m)] = mr.mass
                        row["Gamma_" * String(m)] = mr.gamma
                        row["converged_" * String(m)] = mr.converged
                        row["residual_" * String(m)] = mr.residual

                        if m in (:eta, :eta_prime, :sigma, :sigma_prime)
                            thr = mr.threshold
                            gaps = mr.gaps
                            row["thr_uu_" * String(m)] = thr.uu
                            row["thr_ss_" * String(m)] = thr.ss
                            row["thr_min_" * String(m)] = thr.min
                            row["gap_uu_" * String(m)] = gaps.uu
                            row["gap_ss_" * String(m)] = gaps.ss
                            row["gap_min_" * String(m)] = gaps.min
                        else
                            row["threshold_" * String(m)] = mr.threshold
                            row["gap_" * String(m)] = mr.gap
                        end
                    end

                    vals = _row_to_values(cols, row)
                    println(io, join(vals, ','))
                    push!(existing, key)

                    muB += opts.mubstep_mev
                end
                T += opts.tstep_mev
            end
        end
    end

    println("Wrote scan CSV: ", outpath)
end

main()
