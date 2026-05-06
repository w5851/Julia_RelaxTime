"""
Friesen 2019 Fig.5 top panel charged ratio reproduction.

特点：
- 路径来自外部提取的 Friesen 2019 phase-line CSV；
- 每个路径点只求一次 meson workflow；
- 在同一个 meson_point 上同时后处理 K+/pi+ 与 K-/pi-；
- 全程仍通过 Models workflow 入口，不在脚本层重写物理求解。
"""

function _bootstrap_profile_env!(args::Vector{String})
    pnjl_profile = "friesen2019_poly"
    physics_profile = "friesen2019_base"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--pnjl-profile" && i < length(args)
            pnjl_profile = args[i + 1]
            i += 1
        elseif arg == "--physics-profile" && i < length(args)
            physics_profile = args[i + 1]
            i += 1
        end
        i += 1
    end
    ENV["PNJL_PARAM_PROFILE"] = pnjl_profile
    ENV["PHYSICS_PARAM_PROFILE"] = physics_profile
    return nothing
end

_bootstrap_profile_env!(ARGS)

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models
using Main.Constants_PNJL: ħc_MeV_fm
using Dates

const DEFAULT_TARGET_CSV = raw"D:\Desktop\paper\dev\outputs\formalized\friesen2019_kpi_ratio_curves\friesen2019_fig5_top_charged_ratio_curves_mu055.csv"
const DEFAULT_PHASE_LINE_CSV = raw"D:\Desktop\paper\dev\outputs\formalized\friesen2019_phase_lines\friesen2019_phase_lines.csv"
const DEFAULT_OUTPUT_DIR = joinpath(
    "data", "outputs", "results", "relaxtime", "meson_density", "crossover_validation",
    "friesen2019_fig5_top_mu055_stable_batched",
)

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :target_csv => DEFAULT_TARGET_CSV,
        :phase_line_csv => DEFAULT_PHASE_LINE_CSV,
        :output_dir => DEFAULT_OUTPUT_DIR,
        :source_fig => "fig3_left",
        :case_id => "case1_pnjl",
        :pnjl_profile => get(ENV, "PNJL_PARAM_PROFILE", "friesen2019_poly"),
        :physics_profile => get(ENV, "PHYSICS_PARAM_PROFILE", "friesen2019_base"),
        :xi => 0.0,
        :flavor_chemical_profile => "friesen2019_mu_s_0p55",
        :plus_meson_chemical_profile => "friesen2019_kplus_over_piplus_mu_pi_135_rule",
        :minus_meson_chemical_profile => "friesen2019_kminus_over_piminus_mu_pi_135_rule",
        :overwrite => false,
        :p_num => 12,
        :t_num => 4,
        :max_iter => 20,
        :stable_num_q_nodes => 96,
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
        if arg == "--target-csv"
            opts[:target_csv] = require_value()
        elseif arg == "--phase-line-csv"
            opts[:phase_line_csv] = require_value()
        elseif arg == "--output-dir"
            opts[:output_dir] = require_value()
        elseif arg == "--source-fig"
            opts[:source_fig] = require_value()
        elseif arg == "--case-id"
            opts[:case_id] = require_value()
        elseif arg == "--pnjl-profile"
            opts[:pnjl_profile] = require_value()
        elseif arg == "--physics-profile"
            opts[:physics_profile] = require_value()
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--flavor-chemical-profile"
            opts[:flavor_chemical_profile] = require_value()
        elseif arg == "--plus-meson-chemical-profile"
            opts[:plus_meson_chemical_profile] = require_value()
        elseif arg == "--minus-meson-chemical-profile"
            opts[:minus_meson_chemical_profile] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--stable-q-nodes"
            opts[:stable_num_q_nodes] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/relaxtime/run_friesen2019_fig5_top_charged_scan.jl [options]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return opts
end

function read_phase_line_points(path::String, source_fig::String, case_id::String)
    lines = readlines(path)
    isempty(lines) && error("phase-line CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    idx = Dict(name => i for (i, name) in enumerate(header))
    required = ("source_fig", "case_id", "line_style", "point_index", "muB_MeV", "T_MeV")
    for key in required
        haskey(idx, key) || error("phase-line CSV missing column: $key")
    end

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        strip(cols[idx["source_fig"]]) == source_fig || continue
        strip(cols[idx["case_id"]]) == case_id || continue
        push!(rows, (
            source_fig=strip(cols[idx["source_fig"]]),
            case_id=strip(cols[idx["case_id"]]),
            line_style=strip(cols[idx["line_style"]]),
            point_index=parse(Int, cols[idx["point_index"]]),
            muB_MeV=parse(Float64, cols[idx["muB_MeV"]]),
            T_MeV=parse(Float64, cols[idx["T_MeV"]]),
        ))
    end
    isempty(rows) && error("no phase-line rows found for source_fig=$(source_fig), case_id=$(case_id)")
    sort!(rows; by=r -> (r.muB_MeV, r.line_style == "dashed" ? 0 : 1, r.point_index))
    return rows
end

function read_target_csv(path::String)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    idx = Dict(name => i for (i, name) in enumerate(header))
    xs_plus = Float64[]
    ys_plus = Float64[]
    xs_minus = Float64[]
    ys_minus = Float64[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        curve_id = strip(cols[idx["curve_id"]])
        x = parse(Float64, cols[idx["T_over_muB"]])
        if curve_id == "kplus_over_piplus"
            y = strip(cols[idx["Kplus_over_piplus"]])
            isempty(y) || (push!(xs_plus, x); push!(ys_plus, parse(Float64, y)))
        elseif curve_id == "kminus_over_piminus"
            y = strip(cols[idx["Kminus_over_piminus"]])
            isempty(y) || (push!(xs_minus, x); push!(ys_minus, parse(Float64, y)))
        end
    end
    return (plus=(xs=xs_plus, ys=ys_plus), minus=(xs=xs_minus, ys=ys_minus))
end

function interpolate_rows(rows, x::Float64, field::Symbol)
    n = length(rows)
    n == 0 && return nothing
    x < rows[1].T_over_muB && return nothing
    x > rows[end].T_over_muB && return nothing
    for i in 1:n
        rows[i].T_over_muB == x && return Float64(getproperty(rows[i], field))
    end
    for i in 1:(n - 1)
        left = rows[i]
        right = rows[i + 1]
        if left.T_over_muB <= x <= right.T_over_muB
            dx = right.T_over_muB - left.T_over_muB
            dx > 0.0 || return Float64(getproperty(left, field))
            w = (x - left.T_over_muB) / dx
            return Float64(getproperty(left, field)) + w * (Float64(getproperty(right, field)) - Float64(getproperty(left, field)))
        end
    end
    return nothing
end

function write_comparison_csv(path::String, rows, xs::Vector{Float64}, ys::Vector{Float64}, field::Symbol, label::String)
    open(path, "w") do io
        println(io, "curve_label,x_value,target_y,model_y,abs_diff,rel_diff")
        for (x, y) in zip(xs, ys)
            my = interpolate_rows(rows, x, field)
            if my === nothing
                println(io, join((label, string(x), string(y), "NaN", "NaN", "NaN"), ','))
            else
                ad = abs(my - y)
                rd = iszero(y) ? NaN : ad / abs(y)
                println(io, join((label, string(x), string(y), string(my), string(ad), string(rd)), ','))
            end
        end
    end
end

function summarize(xs::Vector{Float64}, ys::Vector{Float64}, rows, field::Symbol)
    model = Float64[]
    adiffs = Float64[]
    rdiffs = Float64[]
    in_range = 0
    for (x, y) in zip(xs, ys)
        my = interpolate_rows(rows, x, field)
        my === nothing && continue
        in_range += 1
        push!(model, my)
        ad = abs(my - y)
        push!(adiffs, ad)
        iszero(y) || push!(rdiffs, ad / abs(y))
    end
    return (
        target_points=length(xs),
        in_range_points=in_range,
        target_min=minimum(ys),
        target_max=maximum(ys),
        model_min=isempty(model) ? NaN : minimum(model),
        model_max=isempty(model) ? NaN : maximum(model),
        mean_abs_diff=isempty(adiffs) ? NaN : sum(adiffs) / length(adiffs),
        max_abs_diff=isempty(adiffs) ? NaN : maximum(adiffs),
        mean_rel_diff=isempty(rdiffs) ? NaN : sum(rdiffs) / length(rdiffs),
        max_rel_diff=isempty(rdiffs) ? NaN : maximum(rdiffs),
    )
end

function main(args::Vector{String})
    opts = parse_args(args)
    outdir = String(opts[:output_dir])
    if !Bool(opts[:overwrite]) && isdir(outdir)
        error("output_dir already exists: $outdir")
    end
    mkpath(outdir)

    target = read_target_csv(String(opts[:target_csv]))
    flavor_profile = Models.FlavorChemicalProfiles.load_flavor_chemical_profile(profile=String(opts[:flavor_chemical_profile]))
    plus_profile = Models.MesonChemicalProfiles.load_meson_chemical_profile(profile=String(opts[:plus_meson_chemical_profile]))
    minus_profile = Models.MesonChemicalProfiles.load_meson_chemical_profile(profile=String(opts[:minus_meson_chemical_profile]))
    selected = read_phase_line_points(
        String(opts[:phase_line_csv]),
        String(opts[:source_fig]),
        String(opts[:case_id]),
    )

    workflow_csv = joinpath(outdir, "workflow_scan.csv")
    open(workflow_csv, "w") do io
        println(io, "path_source,path_case_id,path_line_style,path_point_index,muq_MeV,muB_MeV,T_MeV,T_over_muB,m_pi_plus_MeV,m_K_plus_MeV,m_pi_minus_MeV,m_K_minus_MeV,n_pi_plus,n_K_plus,kplus_over_piplus,n_pi_minus,n_K_minus,kminus_over_piminus,equilibrium_converged,equilibrium_iterations,equilibrium_residual_norm,message")

        continuation_state = nothing
        for row in selected
            T_MeV = Float64(row.T_MeV)
            muB_MeV = Float64(row.muB_MeV)
            muq_MeV = muB_MeV / 3.0
            T_fm = T_MeV / ħc_MeV_fm
            muq_fm = muq_MeV / ħc_MeV_fm
            flavor_mev = Models.FlavorChemicalProfiles.flavor_mu_profile_MeV(flavor_profile, muq_MeV)
            flavor_fm = Models.FlavorChemicalProfiles.flavor_mu_profile_fm(flavor_profile, muq_MeV)
            plus_chem = Models.MesonChemicalProfiles.meson_chemical_profile_fm(plus_profile; flavor_mev=flavor_mev)
            minus_chem = Models.MesonChemicalProfiles.meson_chemical_profile_fm(minus_profile; flavor_mev=flavor_mev)
            x = muB_MeV > 0.0 ? T_MeV / muB_MeV : NaN
            try
                meson_point = Models.solve_gap_and_meson_point(
                    T_fm,
                    muq_fm;
                    xi=Float64(opts[:xi]),
                    mesons=(:pi_plus, :K_plus, :pi_minus, :K_minus),
                    continuation_state=continuation_state,
                    mixed_branch_align=:strict_sign_binding,
                    flavor_mu_override=(flavor_fm.mu_u_fm, flavor_fm.mu_d_fm, flavor_fm.mu_s_fm),
                    p_num=Int(opts[:p_num]),
                    t_num=Int(opts[:t_num]),
                    solver_kwargs=(iterations=Int(opts[:max_iter]),),
                    mass_kwargs=(iterations=Int(opts[:max_iter]),),
                )
                continuation_state = meson_point.continuation_state

                plus = Models.solve_meson_density_from_meson_point(
                    meson_point;
                    pi_channel=:pi_plus,
                    k_channel=:K_plus,
                    μ_pi=plus_chem.mu_pi_fm,
                    μ_K=plus_chem.mu_K_fm,
                    d_pi=1,
                    d_K=1,
                    num_q_nodes=Int(opts[:stable_num_q_nodes]),
                )
                minus = Models.solve_meson_density_from_meson_point(
                    meson_point;
                    pi_channel=:pi_minus,
                    k_channel=:K_minus,
                    μ_pi=minus_chem.mu_pi_fm,
                    μ_K=minus_chem.mu_K_fm,
                    d_pi=1,
                    d_K=1,
                    num_q_nodes=Int(opts[:stable_num_q_nodes]),
                )

                eq = meson_point.equilibrium
                eq_converged = hasproperty(eq, :converged) ? Bool(getproperty(eq, :converged)) : false
                eq_iterations = hasproperty(eq, :iterations) ? getproperty(eq, :iterations) : -1
                eq_residual = hasproperty(eq, :residual_norm) ? getproperty(eq, :residual_norm) : NaN
                message = ""

                println(io, join((
                    row.source_fig,
                    row.case_id,
                    row.line_style,
                    string(row.point_index),
                    string(muq_MeV),
                    string(muB_MeV),
                    string(T_MeV),
                    string(x),
                    string(plus.m_pi * ħc_MeV_fm),
                    string(plus.m_K * ħc_MeV_fm),
                    string(minus.m_pi * ħc_MeV_fm),
                    string(minus.m_K * ħc_MeV_fm),
                    string(plus.n_pi),
                    string(plus.n_K),
                    string(plus.kpi_ratio),
                    string(minus.n_pi),
                    string(minus.n_K),
                    string(minus.kpi_ratio),
                    string(eq_converged),
                    string(eq_iterations),
                    string(eq_residual),
                    message,
                ), ','))
            catch err
                println(io, join((
                    row.source_fig,
                    row.case_id,
                    row.line_style,
                    string(row.point_index),
                    string(muq_MeV),
                    string(muB_MeV),
                    string(T_MeV),
                    string(x),
                    "NaN","NaN","NaN","NaN",
                    "NaN","NaN","NaN","NaN","NaN","NaN",
                    "false","-1","NaN",
                    replace(sprint(showerror, err), ',' => ';'),
                ), ','))
            end
        end
    end

    lines = readlines(workflow_csv)
    header = split(strip(lines[1]), ',')
    idx = Dict(name => i for (i, name) in enumerate(header))
    rows = NamedTuple[]
    for line in lines[2:end]
        cols = split(strip(line), ',')
        isempty(cols[1]) && continue
        x = try parse(Float64, cols[idx["T_over_muB"]]) catch; NaN end
        yp = try parse(Float64, cols[idx["kplus_over_piplus"]]) catch; NaN end
        ym = try parse(Float64, cols[idx["kminus_over_piminus"]]) catch; NaN end
        isfinite(x) && isfinite(yp) && isfinite(ym) || continue
        push!(rows, (T_over_muB=x, kplus_over_piplus=yp, kminus_over_piminus=ym))
    end
    sort!(rows; by=r -> r.T_over_muB)

    plus_cmp = joinpath(outdir, "comparison_kplus_over_piplus.csv")
    minus_cmp = joinpath(outdir, "comparison_kminus_over_piminus.csv")
    write_comparison_csv(plus_cmp, rows, target.plus.xs, target.plus.ys, :kplus_over_piplus, "K+/pi+")
    write_comparison_csv(minus_cmp, rows, target.minus.xs, target.minus.ys, :kminus_over_piminus, "K-/pi-")

    plus_summary = summarize(target.plus.xs, target.plus.ys, rows, :kplus_over_piplus)
    minus_summary = summarize(target.minus.xs, target.minus.ys, rows, :kminus_over_piminus)
    readme = joinpath(outdir, "README.md")
    open(readme, "w") do io
        println(io, "# Friesen 2019 Fig.5 top charged scan")
        println(io)
        println(io, "更新日期：$(Dates.today())")
        println(io)
        println(io, "## 运行口径")
        println(io)
        println(io, "- phase-line csv: `$(opts[:phase_line_csv])`")
        println(io, "- source_fig / case_id: `$(opts[:source_fig]) / $(opts[:case_id])`")
        println(io, "- pnjl / physics profile: `$(opts[:pnjl_profile]) / $(opts[:physics_profile])`")
        println(io, "- flavor profile: `$(opts[:flavor_chemical_profile])`")
        println(io, "- plus meson profile: `$(opts[:plus_meson_chemical_profile])`")
        println(io, "- minus meson profile: `$(opts[:minus_meson_chemical_profile])`")
        println(io, "- output: `workflow_scan.csv`, `comparison_kplus_over_piplus.csv`, `comparison_kminus_over_piminus.csv`")
        println(io)
        println(io, "## K+/pi+")
        println(io, "- target points: `$(plus_summary.target_points)`")
        println(io, "- in-range points: `$(plus_summary.in_range_points)`")
        println(io, "- target range: `$(plus_summary.target_min)` -> `$(plus_summary.target_max)`")
        println(io, "- model range: `$(plus_summary.model_min)` -> `$(plus_summary.model_max)`")
        println(io, "- mean rel diff: `$(plus_summary.mean_rel_diff)`")
        println(io, "- max rel diff: `$(plus_summary.max_rel_diff)`")
        println(io)
        println(io, "## K-/pi-")
        println(io, "- target points: `$(minus_summary.target_points)`")
        println(io, "- in-range points: `$(minus_summary.in_range_points)`")
        println(io, "- target range: `$(minus_summary.target_min)` -> `$(minus_summary.target_max)`")
        println(io, "- model range: `$(minus_summary.model_min)` -> `$(minus_summary.model_max)`")
        println(io, "- mean rel diff: `$(minus_summary.mean_rel_diff)`")
        println(io, "- max rel diff: `$(minus_summary.max_rel_diff)`")
    end

    println("friesen fig5 top charged scan completed")
    println("output_dir=$(outdir)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
