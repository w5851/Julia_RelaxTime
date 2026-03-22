"""
Plot K-mass comparison for literature vs Julia vs fortran_mott (if available).

Usage:
  julia --project=. scripts/analysis/plot_k_mass_literature_julia_fortranmott.jl
  julia --project=. scripts/analysis/plot_k_mass_literature_julia_fortranmott.jl --fortran-csv <path>
"""

using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const _MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

const DEFAULT_TARGET_CSV = joinpath(
    PROJECT_ROOT,
    "tests", "validation", "data", "targets", "relaxtime", "literature", "meson",
    "relaxtime_meson_mass_literature_targets_v1.csv",
)

const DEFAULT_FORTRAN_CSVS = [
    joinpath(
        PROJECT_ROOT,
        "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson",
        "legacy_meson_scan_fortran_muB0_v1.csv",
    ),
    joinpath(
        PROJECT_ROOT,
        "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson",
        "legacy_meson_scan_fortran_muB600_v1.csv",
    ),
]

const DEFAULT_OUT_CSV = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "results", "relaxtime", "literature",
    "k_mass_literature_julia_fortranmott_comparison.csv",
)

const DEFAULT_OUT_FIG = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "figures", "relaxtime", "literature",
    "k_mass_literature_julia_fortranmott_comparison.png",
)

function _parse_args(args::Vector{String})
    target_csv = DEFAULT_TARGET_CSV
    fortran_csv = join(DEFAULT_FORTRAN_CSVS, ";")
    out_csv = DEFAULT_OUT_CSV
    out_fig = DEFAULT_OUT_FIG

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--target-csv"
            i == length(args) && error("missing value for --target-csv")
            i += 1
            target_csv = args[i]
        elseif arg == "--fortran-csv"
            i == length(args) && error("missing value for --fortran-csv")
            i += 1
            fortran_csv = args[i]
        elseif arg == "--out-csv"
            i == length(args) && error("missing value for --out-csv")
            i += 1
            out_csv = args[i]
        elseif arg == "--out-fig"
            i == length(args) && error("missing value for --out-fig")
            i += 1
            out_fig = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/analysis/plot_k_mass_literature_julia_fortranmott.jl [--target-csv <path>] [--fortran-csv <path>] [--out-csv <path>] [--out-fig <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return target_csv, fortran_csv, out_csv, out_fig
end

function _load_literature_k_targets(path::String)
    isfile(path) || error("missing target CSV: $path")
    lines = readlines(path)
    isempty(lines) && error("empty target CSV: $path")

    rows = NamedTuple[]
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        startswith(s, "#") && continue
        c = split(s, ',')
        length(c) == 9 || continue
        meson = Symbol(strip(c[3]))
        meson == :K || continue
        push!(rows, (
            muB_MeV=parse(Float64, strip(c[4])),
            xi=parse(Float64, strip(c[5])),
            T_MeV=parse(Float64, strip(c[6])),
            literature_mass_MeV=parse(Float64, strip(c[7])),
        ))
    end
    return rows
end

function _load_fortran_k(paths_spec::String)
    out = Dict{Tuple{Float64,Float64},Float64}()
    paths = [strip(p) for p in split(paths_spec, ';') if !isempty(strip(p))]
    for path in paths
        isfile(path) || continue
        lines = readlines(path)
        isempty(lines) && continue
        header = split(strip(lines[1]), ',')
        idx = Dict(name => i for (i, name) in enumerate(header))
        needed = ("meson", "muB_MeV", "xi", "T_MeV", "mass_MeV")
        all(k -> haskey(idx, k), needed) || continue

        for ln in lines[2:end]
            s = strip(ln)
            isempty(s) && continue
            c = split(s, ',')
            length(c) == length(header) || continue
            strip(c[idx["meson"]]) == "K" || continue
            muB = parse(Float64, strip(c[idx["muB_MeV"]]))
            xi = parse(Float64, strip(c[idx["xi"]]))
            isapprox(xi, 0.0; atol=1e-12) || continue
            T = parse(Float64, strip(c[idx["T_MeV"]]))
            mass = parse(Float64, strip(c[idx["mass_MeV"]]))
            key = (muB, T)
            haskey(out, key) || (out[key] = mass)
        end
    end
    return out
end

function _julia_k_mass(T_MeV::Float64, muB_MeV::Float64)
    T_fm = T_MeV / _HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / _HBARC_MEV_FM
    res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        muq_fm;
        xi=0.0,
        mesons=(:K,),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=25,),
    )
    return Float64(res.meson_results[:K].mass) * _HBARC_MEV_FM
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "muB_MeV,T_MeV,literature_K_MeV,julia_K_MeV,fortran_mott_K_MeV,rel_diff_julia_lit,rel_diff_fortran_lit")
        for r in rows
            println(io, join([
                string(r.muB_MeV),
                string(r.T_MeV),
                string(r.lit),
                string(r.julia),
                isnan(r.fortran) ? "NaN" : string(r.fortran),
                string(r.rel_julia),
                isnan(r.rel_fortran) ? "NaN" : string(r.rel_fortran),
            ], ','))
        end
    end
end

function main(args::Vector{String})
    target_csv, fortran_csv, out_csv, out_fig = _parse_args(args)

    targets = _load_literature_k_targets(target_csv)
    isempty(targets) && error("no K targets in: $target_csv")
    fortran_map = _load_fortran_k(fortran_csv)

    rows = NamedTuple[]
    for t in sort(targets; by=x -> (x.muB_MeV, x.T_MeV))
        julia_mass = _julia_k_mass(t.T_MeV, t.muB_MeV)
        key = (t.muB_MeV, t.T_MeV)
        fortran_mass = get(fortran_map, key, NaN)
        rel_j = abs(julia_mass - t.literature_mass_MeV) / max(abs(t.literature_mass_MeV), 1e-12)
        rel_f = isfinite(fortran_mass) ? abs(fortran_mass - t.literature_mass_MeV) / max(abs(t.literature_mass_MeV), 1e-12) : NaN
        push!(rows, (
            muB_MeV=t.muB_MeV,
            T_MeV=t.T_MeV,
            lit=t.literature_mass_MeV,
            julia=julia_mass,
            fortran=fortran_mass,
            rel_julia=rel_j,
            rel_fortran=rel_f,
        ))
    end

    _write_csv(out_csv, rows)

    function _series(muB::Float64)
        s = filter(r -> isapprox(r.muB_MeV, muB; atol=1e-12), rows)
        Ts = [r.T_MeV for r in s]
        lit = [r.lit for r in s]
        julia = [r.julia for r in s]
        fortran = [r.fortran for r in s]
        rel_j = [100.0 * r.rel_julia for r in s]
        rel_f = [100.0 * r.rel_fortran for r in s]
        return Ts, lit, julia, fortran, rel_j, rel_f
    end

    T0, L0, J0, F0, RJ0, RF0 = _series(0.0)
    T6, L6, J6, F6, RJ6, RF6 = _series(600.0)

    p1 = plot(T0, L0; label="lit muB=0", lw=2.5, marker=:circle, color=:black,
        xlabel="T [MeV]", ylabel="K mass [MeV]", title="K mass: muB=0", legend=:topleft)
    plot!(p1, T0, J0; label="julia muB=0", lw=2.5, marker=:diamond, color=:dodgerblue)
    if any(isfinite, F0)
        plot!(p1, T0, F0; label="fortran_mott muB=0", lw=2.0, ls=:dash, marker=:star5, color=:darkgreen)
    end

    p2 = plot(T6, L6; label="lit muB=600", lw=2.5, marker=:circle, color=:black,
        xlabel="T [MeV]", ylabel="K mass [MeV]", title="K mass: muB=600", legend=:topleft)
    plot!(p2, T6, J6; label="julia muB=600", lw=2.5, marker=:diamond, color=:firebrick)
    if any(isfinite, F6)
        plot!(p2, T6, F6; label="fortran_mott muB=600", lw=2.0, ls=:dash, marker=:star5, color=:darkgreen)
    else
        annotate!(p2, (T6[1], maximum(L6), text("fortran_mott muB=600 unavailable", 9, :darkgreen)))
    end

    p3 = plot(T0, RJ0; label="julia-lit (muB=0)", lw=2.5, marker=:utriangle, color=:dodgerblue,
        xlabel="T [MeV]", ylabel="rel diff [%]", title="Relative difference vs literature", legend=:topright)
    if any(isfinite, RF0)
        plot!(p3, T0, RF0; label="fortran-lit (muB=0)", lw=2.0, ls=:dash, marker=:utriangle, color=:darkgreen)
    end
    plot!(p3, T6, RJ6; label="julia-lit (muB=600)", lw=2.0, ls=:dash, marker=:utriangle, color=:firebrick)
    if any(isfinite, RF6)
        plot!(p3, T6, RF6; label="fortran-lit (muB=600)", lw=2.0, ls=:dot, marker=:utriangle, color=:olive)
    end

    fig = plot(p1, p2, p3; layout=(3, 1), size=(980, 1180))
    mkpath(dirname(out_fig))
    savefig(fig, out_fig)

    println("Wrote comparison CSV: " * out_csv)
    println("Wrote comparison figure: " * out_fig)
    println("fortran_mott rows loaded: " * string(length(fortran_map)))
end

main(ARGS)
