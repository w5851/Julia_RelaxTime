"""
Phase E5 / meson-density point diagnostic workflow probe.

目标：
- 直接走当前 `Models.solve_gap_and_meson_point -> Models.solve_phase_shift_point_diagnostic_from_meson_point`
  workflow 合同；
- 在固定 `(T, xi, q, ω)` 点上导出 `AD / formula / FD` 对照；
- 判断当前 AD 结果更像“完整导数”还是“冻结积分拓扑后的局部导数”。

输出：
- CSV: `data/outputs/results/relaxtime/analysis/meson_density_phase_point_ad_vs_fd.csv`
- Markdown summary: `data/outputs/results/relaxtime/analysis/meson_density_phase_point_ad_vs_fd_summary.md`
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models: solve_gap_and_meson_point, solve_phase_shift_point_diagnostic_from_meson_point

const DEFAULT_CSV = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "meson_density_phase_point_ad_vs_fd.csv",
)

const DEFAULT_MD = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "meson_density_phase_point_ad_vs_fd_summary.md",
)

const DEFAULT_T_MEV = [208.0, 210.0, 212.0]
const DEFAULT_XI = [0.0, 0.1]
const DEFAULT_Q = [0.0, 0.5, 1.0]
const DEFAULT_OMEGA = [0.3, 0.6, 1.0]
const DEFAULT_MESONS = (:pi, :K)
const DEFAULT_SCHEME = :current
const DEFAULT_FD_STEP = 1e-5
const DEFAULT_ETA = 1e-6

@inline _fmt_bool(x::Bool) = x ? "true" : "false"
@inline _fmt(x) = x isa Bool ? _fmt_bool(x) : string(x)

function _parse_float_list(raw::AbstractString)
    vals = Float64[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("parsed float list is empty"))
    return vals
end

function _parse_symbol_list(raw::AbstractString)
    vals = Symbol[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, Symbol(s))
    end
    isempty(vals) && throw(ArgumentError("parsed symbol list is empty"))
    return Tuple(vals)
end

_selected_t_mev() = isempty(strip(get(ENV, "MESON_PHASE_AD_T_MEV", ""))) ? DEFAULT_T_MEV : _parse_float_list(ENV["MESON_PHASE_AD_T_MEV"])
_selected_xi() = isempty(strip(get(ENV, "MESON_PHASE_AD_XI", ""))) ? DEFAULT_XI : _parse_float_list(ENV["MESON_PHASE_AD_XI"])
_selected_q() = isempty(strip(get(ENV, "MESON_PHASE_AD_Q", ""))) ? DEFAULT_Q : _parse_float_list(ENV["MESON_PHASE_AD_Q"])
_selected_omega() = isempty(strip(get(ENV, "MESON_PHASE_AD_OMEGA", ""))) ? DEFAULT_OMEGA : _parse_float_list(ENV["MESON_PHASE_AD_OMEGA"])
_selected_mesons() = isempty(strip(get(ENV, "MESON_PHASE_AD_MESONS", ""))) ? DEFAULT_MESONS : _parse_symbol_list(ENV["MESON_PHASE_AD_MESONS"])

function _write_csv(path::String, rows::Vector{NamedTuple})
    isempty(rows) && error("no rows to write")
    headers = collect(keys(rows[1]))
    open(path, "w") do io
        println(io, join(string.(headers), ","))
        for row in rows
            vals = [_fmt(getproperty(row, h)) for h in headers]
            println(io, join(vals, ","))
        end
    end
end

function _failed_row(T_MeV, T_fm, xi, meson, q, ω, err)
    return (
        T_MeV=T_MeV,
        T_fm=T_fm,
        converged=true,
        status="error",
        error=replace(sprint(showerror, err), ',' => ';'),
        meson=meson,
        omega=Float64(ω),
        q=Float64(q),
        xi=Float64(xi),
        temperature=Float64(T_fm),
        eta=DEFAULT_ETA,
        fd_step=DEFAULT_FD_STEP,
        scheme=DEFAULT_SCHEME,
        reD=NaN,
        imD=NaN,
        phase=NaN,
        weighted_phase=NaN,
        dReD_domega=NaN,
        dImD_domega=NaN,
        dReD_fd=NaN,
        dImD_fd=NaN,
        dphase_ad=NaN,
        dphase_formula=NaN,
        dphase_fd=NaN,
        dphase_raw_fd=NaN,
        dphase_abs_diff=NaN,
        dphase_formula_abs_diff=NaN,
        dphase_raw_fd_abs_diff=NaN,
        dweighted_ad=NaN,
        dweighted_fd=NaN,
        dweighted_abs_diff=NaN,
    )
end

function _collect_rows(; scheme::Symbol, fd_step::Float64, eta::Float64)
    all_rows = NamedTuple[]
    for xi in _selected_xi(), T_MeV in _selected_t_mev()
        T_fm = T_MeV / ħc_MeV_fm
        meson_point = solve_gap_and_meson_point(
            T_fm,
            0.0;
            xi=xi,
            p_num=12,
            t_num=6,
        )
        Bool(meson_point.equilibrium.converged) || error("workflow point did not converge at T=$(T_MeV) MeV, xi=$(xi)")

        for meson in _selected_mesons(), q in _selected_q(), ω in _selected_omega()
            try
                diag = solve_phase_shift_point_diagnostic_from_meson_point(
                    meson_point;
                    mesons=(meson,),
                    q_values=[q],
                    omega_values=[ω],
                    scheme=scheme,
                    eta=eta,
                    fd_step=fd_step,
                )
                row = only(diag.rows)
                push!(all_rows, merge((
                    T_MeV=T_MeV,
                    T_fm=T_fm,
                    converged=Bool(meson_point.equilibrium.converged),
                    status="ok",
                    error="",
                ), row))
            catch err
                push!(all_rows, _failed_row(T_MeV, T_fm, xi, meson, q, ω, err))
            end
        end
    end
    return all_rows
end

function _max_by(rows::Vector{NamedTuple}, pred)
    best_idx = 1
    best_val = pred(rows[1])
    for i in 2:length(rows)
        val = pred(rows[i])
        if val > best_val
            best_idx = i
            best_val = val
        end
    end
    return rows[best_idx], best_val
end

function _summarize(rows::Vector{NamedTuple}, md_path::String; scheme::Symbol, fd_step::Float64, eta::Float64)
    ok_rows = filter(r -> r.status == "ok", rows)
    err_rows = filter(r -> r.status != "ok", rows)
    isempty(ok_rows) && error("no successful rows to summarize")
    max_formula_row, max_formula = _max_by(ok_rows, r -> r.dphase_formula_abs_diff)
    max_phase_row, max_phase = _max_by(ok_rows, r -> r.dphase_abs_diff)
    max_raw_row, max_raw = _max_by(ok_rows, r -> r.dphase_raw_fd_abs_diff)
    max_re_row, max_re = _max_by(ok_rows, r -> abs(r.dReD_domega - r.dReD_fd))
    max_im_row, max_im = _max_by(ok_rows, r -> abs(r.dImD_domega - r.dImD_fd))

    open(md_path, "w") do io
        println(io, "# Meson Density Phase Point AD vs FD Workflow Probe")
        println(io)
        println(io, "## Configuration")
        println(io)
        println(io, "- scheme: `", scheme, "`")
        println(io, "- eta: `", eta, "`")
        println(io, "- fd_step: `", fd_step, "`")
        println(io, "- T_MeV: `", join(_selected_t_mev(), ", "), "`")
        println(io, "- xi: `", join(_selected_xi(), ", "), "`")
        println(io, "- q: `", join(_selected_q(), ", "), "`")
        println(io, "- omega: `", join(_selected_omega(), ", "), "`")
        println(io, "- mesons: `", join(string.(collect(_selected_mesons())), ", "), "`")
        println(io)
        println(io, "## Top Findings")
        println(io)
        println(io, "- successful rows: `", length(ok_rows), "`")
        println(io, "- failed rows: `", length(err_rows), "`")
        println(io, "- `dphase_formula_abs_diff` max = `", max_formula, "` at `(T_MeV=", max_formula_row.T_MeV, ", xi=", max_formula_row.xi, ", meson=", max_formula_row.meson, ", q=", max_formula_row.q, ", omega=", max_formula_row.omega, ")`")
        println(io, "- `|dphase_ad - dphase_fd|` max = `", max_phase, "` at `(T_MeV=", max_phase_row.T_MeV, ", xi=", max_phase_row.xi, ", meson=", max_phase_row.meson, ", q=", max_phase_row.q, ", omega=", max_phase_row.omega, ")`")
        println(io, "- `|dphase_ad - dphase_raw_fd|` max = `", max_raw, "` at `(T_MeV=", max_raw_row.T_MeV, ", xi=", max_raw_row.xi, ", meson=", max_raw_row.meson, ", q=", max_raw_row.q, ", omega=", max_raw_row.omega, ")`")
        println(io, "- `|dReD_domega - dReD_fd|` max = `", max_re, "` at `(T_MeV=", max_re_row.T_MeV, ", xi=", max_re_row.xi, ", meson=", max_re_row.meson, ", q=", max_re_row.q, ", omega=", max_re_row.omega, ")`")
        println(io, "- `|dImD_domega - dImD_fd|` max = `", max_im, "` at `(T_MeV=", max_im_row.T_MeV, ", xi=", max_im_row.xi, ", meson=", max_im_row.meson, ", q=", max_im_row.q, ", omega=", max_im_row.omega, ")`")
        println(io)
        println(io, "## Reading")
        println(io)
        println(io, "1. `dphase_formula_abs_diff` should stay near machine precision if the AD path and the `Re/Im D` analytic derivative formula are internally consistent.")
        println(io, "2. If `dphase_abs_diff` remains visibly nonzero while (1) is tiny, the main gap is not `atan` branch handling alone; it already exists at the `ReD/ImD` derivative level.")
        println(io, "3. In the current implementation, that usually means AD is differentiating a frozen-topology integrand, while FD also samples how singularity intervals / segmentation move with `omega`.")
        if !isempty(err_rows)
            println(io)
            println(io, "## Unsupported Points")
            println(io)
            for row in err_rows
                println(io, "- `(T_MeV=", row.T_MeV, ", xi=", row.xi, ", meson=", row.meson, ", q=", row.q, ", omega=", row.omega, ")`: `", row.error, "`")
            end
        end
    end
end

function main()
    csv_path = isempty(ARGS) ? DEFAULT_CSV : abspath(ARGS[1])
    md_path = length(ARGS) >= 2 ? abspath(ARGS[2]) : DEFAULT_MD
    mkpath(dirname(csv_path))
    mkpath(dirname(md_path))

    rows = _collect_rows(
        scheme=DEFAULT_SCHEME,
        fd_step=DEFAULT_FD_STEP,
        eta=DEFAULT_ETA,
    )
    _write_csv(csv_path, rows)
    _summarize(rows, md_path; scheme=DEFAULT_SCHEME, fd_step=DEFAULT_FD_STEP, eta=DEFAULT_ETA)

    println("Wrote workflow AD/FD CSV: ", csv_path)
    println("Wrote workflow AD/FD summary: ", md_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
