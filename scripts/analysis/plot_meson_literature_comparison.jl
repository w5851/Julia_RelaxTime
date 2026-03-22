"""
Plot Julia vs digitized literature meson masses (pi/K).

Usage:
  julia --project=. scripts/analysis/plot_meson_literature_comparison.jl
  julia --project=. scripts/analysis/plot_meson_literature_comparison.jl --input <targets.csv> --out-csv <cmp.csv> --out-fig <cmp.png>
"""

using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const _MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

struct Options
    input_csv::String
    output_csv::String
    output_fig::String
end

function parse_args(args::Vector{String})
    input_csv = joinpath(
        PROJECT_ROOT,
        "tests", "validation", "data", "targets", "relaxtime", "literature", "meson",
        "relaxtime_meson_mass_literature_targets_v1.csv",
    )
    output_csv = joinpath(
        PROJECT_ROOT,
        "data", "outputs", "results", "relaxtime", "literature",
        "meson_mass_julia_vs_literature_comparison.csv",
    )
    output_fig = joinpath(
        PROJECT_ROOT,
        "data", "outputs", "figures", "relaxtime",
        "meson_mass_julia_vs_literature_comparison.png",
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input"
            i == length(args) && error("missing value for --input")
            i += 1
            input_csv = args[i]
        elseif arg == "--out-csv"
            i == length(args) && error("missing value for --out-csv")
            i += 1
            output_csv = args[i]
        elseif arg == "--out-fig"
            i == length(args) && error("missing value for --out-fig")
            i += 1
            output_fig = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/analysis/plot_meson_literature_comparison.jl [--input <targets.csv>] [--out-csv <cmp.csv>] [--out-fig <cmp.png>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return Options(input_csv, output_csv, output_fig)
end

function _load_targets(path::String)
    isfile(path) || error("target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            series=strip(cols[2]),
            meson=Symbol(strip(cols[3])),
            muB_MeV=parse(Float64, strip(cols[4])),
            xi=parse(Float64, strip(cols[5])),
            T_MeV=parse(Float64, strip(cols[6])),
            expected_mass_MeV=parse(Float64, strip(cols[7])),
            rtol=parse(Float64, strip(cols[8])),
            source=strip(cols[9]),
        ))
    end
    return rows
end

function _compute_mass_mev(T_MeV::Float64, muB_MeV::Float64, xi::Float64, meson::Symbol)
    T_fm = T_MeV / _HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / _HBARC_MEV_FM
    res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        muq_fm;
        xi=xi,
        mesons=(meson,),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=25,),
    )
    return Float64(res.meson_results[meson].mass) * _HBARC_MEV_FM
end

function _write_comparison_csv(path::String, rows)
    mkpath(dirname(path))
    header = [
        "target_id", "series", "meson", "muB_MeV", "xi", "T_MeV",
        "literature_mass_MeV", "julia_mass_MeV", "abs_diff_MeV", "rel_diff", "rtol",
    ]
    open(path, "w") do io
        println(io, join(header, ','))
        for r in rows
            println(io, join([
                r.target_id,
                r.series,
                String(r.meson),
                string(r.muB_MeV),
                string(r.xi),
                string(r.T_MeV),
                string(r.literature_mass_MeV),
                string(r.julia_mass_MeV),
                string(r.abs_diff_MeV),
                string(r.rel_diff),
                string(r.rtol),
            ], ','))
        end
    end
end

function _plot_and_save(path::String, rows)
    mkpath(dirname(path))
    function _series(meson::Symbol, muB::Float64)
        s = sort(filter(r -> r.meson == meson && isapprox(r.muB_MeV, muB; atol=1e-12), rows); by=r -> r.T_MeV)
        Ts = [r.T_MeV for r in s]
        y_lit = [r.literature_mass_MeV for r in s]
        y_julia = [r.julia_mass_MeV for r in s]
        y_rel_pct = [100.0 * r.rel_diff for r in s]
        return Ts, y_lit, y_julia, y_rel_pct
    end

    Ts_pi0, pi_l0, pi_j0, pi_r0 = _series(:pi, 0.0)
    Ts_pi6, pi_l6, pi_j6, pi_r6 = _series(:pi, 600.0)
    Ts_k0, k_l0, k_j0, k_r0 = _series(:K, 0.0)
    Ts_k6, k_l6, k_j6, k_r6 = _series(:K, 600.0)

    p_pi_mass = plot(Ts_pi0, pi_l0;
        label="lit muB=0",
        lw=2.5,
        marker=:circle,
        color=:black,
        xlabel="T [MeV]",
        ylabel="mass [MeV]",
        title="pi mass",
        legend=:topleft,
    )
    plot!(p_pi_mass, Ts_pi0, pi_j0; label="julia muB=0", lw=2.5, marker=:diamond, color=:dodgerblue)
    plot!(p_pi_mass, Ts_pi6, pi_l6; label="lit muB=600", lw=2.0, ls=:dash, marker=:circle, color=:gray40)
    plot!(p_pi_mass, Ts_pi6, pi_j6; label="julia muB=600", lw=2.0, ls=:dash, marker=:diamond, color=:firebrick)

    p_pi_rel = plot(Ts_pi0, pi_r0;
        label="muB=0",
        lw=2.5,
        marker=:utriangle,
        color=:dodgerblue,
        xlabel="T [MeV]",
        ylabel="rel diff [%]",
        title="pi relative difference",
        legend=:topright,
    )
    plot!(p_pi_rel, Ts_pi6, pi_r6; label="muB=600", lw=2.0, ls=:dash, marker=:utriangle, color=:firebrick)

    p_k_mass = plot(Ts_k0, k_l0;
        label="lit muB=0",
        lw=2.5,
        marker=:circle,
        color=:black,
        xlabel="T [MeV]",
        ylabel="mass [MeV]",
        title="K mass",
        legend=:topleft,
    )
    plot!(p_k_mass, Ts_k0, k_j0; label="julia muB=0", lw=2.5, marker=:diamond, color=:dodgerblue)
    plot!(p_k_mass, Ts_k6, k_l6; label="lit muB=600", lw=2.0, ls=:dash, marker=:circle, color=:gray40)
    plot!(p_k_mass, Ts_k6, k_j6; label="julia muB=600", lw=2.0, ls=:dash, marker=:diamond, color=:firebrick)

    p_k_rel = plot(Ts_k0, k_r0;
        label="muB=0",
        lw=2.5,
        marker=:utriangle,
        color=:dodgerblue,
        xlabel="T [MeV]",
        ylabel="rel diff [%]",
        title="K relative difference",
        legend=:topright,
    )
    plot!(p_k_rel, Ts_k6, k_r6; label="muB=600", lw=2.0, ls=:dash, marker=:utriangle, color=:firebrick)

    fig = plot(p_pi_mass, p_pi_rel, p_k_mass, p_k_rel; layout=(2, 2), size=(1200, 820))
    savefig(fig, path)
end

function main(args::Vector{String})
    opts = parse_args(args)
    targets = _load_targets(opts.input_csv)
    isempty(targets) && error("no targets loaded from $(opts.input_csv)")

    rows = NamedTuple[]
    for t in targets
        julia_mass = _compute_mass_mev(t.T_MeV, t.muB_MeV, t.xi, t.meson)
        lit_mass = t.expected_mass_MeV
        abs_diff = abs(julia_mass - lit_mass)
        rel_diff = abs_diff / max(abs(lit_mass), 1e-12)
        push!(rows, (
            target_id=t.target_id,
            series=t.series,
            meson=t.meson,
            muB_MeV=t.muB_MeV,
            xi=t.xi,
            T_MeV=t.T_MeV,
            literature_mass_MeV=lit_mass,
            julia_mass_MeV=julia_mass,
            abs_diff_MeV=abs_diff,
            rel_diff=rel_diff,
            rtol=t.rtol,
        ))
    end

    _write_comparison_csv(opts.output_csv, rows)
    _plot_and_save(opts.output_fig, rows)

    println("Wrote comparison CSV: " * opts.output_csv)
    println("Wrote comparison figure: " * opts.output_fig)
end

main(ARGS)
