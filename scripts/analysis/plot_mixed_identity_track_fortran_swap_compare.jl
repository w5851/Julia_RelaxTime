#!/usr/bin/env julia

using DelimitedFiles
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _read_csv_rows(path::String)
    lines = readlines(path)
    length(lines) >= 2 || error("csv has no data rows: $path")
    header = split(lines[1], ',')
    rows = Vector{NamedTuple}()
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) == length(header) || continue
        push!(rows, (
            T_MeV=parse(Float64, cols[1]),
            slot=Symbol(cols[2]),
            mapping=cols[3],
            j_mass_MeV=parse(Float64, cols[4]),
            j_gamma_MeV=parse(Float64, cols[5]),
            f_mass_MeV=parse(Float64, cols[6]),
            f_gamma_MeV=parse(Float64, cols[7]),
            mass_rel=parse(Float64, cols[8]),
            gamma_rel_floor1=parse(Float64, cols[9]),
            gamma_rel_raw=parse(Float64, cols[10]),
        ))
    end
    return rows
end

function _series(rows, slot::Symbol)
    s = sort(filter(r -> r.slot === slot, rows); by=r -> r.T_MeV)
    T = [r.T_MeV for r in s]
    jm = [r.j_mass_MeV for r in s]
    fm = [r.f_mass_MeV for r in s]
    jg = [r.j_gamma_MeV for r in s]
    fg = [r.f_gamma_MeV for r in s]
    return T, jm, fm, jg, fg
end

function main()
    in_csv = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics", "fortran_julia_identity_track_swap_compare_muB600.csv")
    isfile(in_csv) || error("input csv missing: $in_csv")

    rows = _read_csv_rows(in_csv)
    isempty(rows) && error("no rows in $in_csv")

    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "validation")
    mkpath(out_dir)

    T_eta, jm_eta, fm_eta, jg_eta, fg_eta = _series(rows, :eta)
    T_etap, jm_etap, fm_etap, jg_etap, fg_etap = _series(rows, :eta_prime)

    p1 = plot(T_eta, fm_eta; lw=2.2, color=:black, label="Fortran eta", xlabel="T [MeV]", ylabel="mass [MeV]", title="eta mass: Julia vs Fortran(swap-validated)")
    plot!(p1, T_eta, jm_eta; lw=2.0, color=:dodgerblue, ls=:dash, label="Julia eta")

    p2 = plot(T_etap, fm_etap; lw=2.2, color=:black, label="Fortran eta_prime", xlabel="T [MeV]", ylabel="mass [MeV]", title="eta' mass: Julia vs Fortran(swap-validated)")
    plot!(p2, T_etap, jm_etap; lw=2.0, color=:crimson, ls=:dash, label="Julia eta_prime")

    p3 = plot(T_eta, fg_eta; lw=2.2, color=:black, label="Fortran eta", xlabel="T [MeV]", ylabel="gamma [MeV]", title="eta gamma: Julia vs Fortran(swap-validated)")
    plot!(p3, T_eta, jg_eta; lw=2.0, color=:dodgerblue, ls=:dash, label="Julia eta")

    p4 = plot(T_etap, fg_etap; lw=2.2, color=:black, label="Fortran eta_prime", xlabel="T [MeV]", ylabel="gamma [MeV]", title="eta' gamma: Julia vs Fortran(swap-validated)")
    plot!(p4, T_etap, jg_etap; lw=2.0, color=:crimson, ls=:dash, label="Julia eta_prime")

    out_png = joinpath(out_dir, "mixed_identity_track_fortran_swap_compare_muB600.png")
    savefig(plot(p1, p2, p3, p4; layout=(2, 2), size=(1300, 900)), out_png)
    println("Wrote comparison figure: $out_png")
end

main()
