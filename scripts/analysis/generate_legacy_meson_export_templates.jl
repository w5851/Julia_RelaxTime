"""
Generate legacy meson export CSV template for Fortran.

This script creates schema-complete CSV files under:
  tests/validation/data/targets/relaxtime/legacy/meson/

Values are placeholders (NaN/0/fail) and should be replaced with
actual legacy outputs.

Usage:
  julia --project=. scripts/analysis/generate_legacy_meson_export_templates.jl
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const OUTDIR = joinpath(
    PROJECT_ROOT,
    "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson",
)

const MANDATORY_COLS = (
    "record_id",
    "source_impl",
    "T_MeV",
    "muB_MeV",
    "xi",
    "meson",
    "mass_MeV",
    "threshold_MeV",
    "gap_MeV",
    "mott_flag",
    "solver_status",
)

const MESONS = (
    "pi",
    "K",
    "eta",
    "eta_prime",
    "sigma_pi",
    "sigma_K",
    "sigma",
    "sigma_prime",
)

const T_ISO = (120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0)
const T_ANISO = (160.0, 180.0, 200.0, 220.0)

function xi_grid()
    xs = Float64[]
    x = -0.4
    while x <= 0.4 + 1e-12
        push!(xs, round(x; digits=10))
        x += 0.05
    end
    return xs
end

function placeholder_row(source_impl::String, T_MeV::Float64, xi::Float64, meson::String, idx::Int)
    rid = string(source_impl, "_", meson, "_", replace(string(T_MeV), "." => "p"), "_xi_", replace(string(xi), "." => "p"), "_", idx)
    return (
        rid,
        source_impl,
        string(T_MeV),
        "0.0",
        string(xi),
        meson,
        "NaN",
        "NaN",
        "NaN",
        "0",
        "fail",
    )
end

function write_template(path::String, source_impl::String)
    mkpath(dirname(path))
    rows = Vector{NTuple{11,String}}()
    idx = 0

    for T in T_ISO
        for meson in MESONS
            idx += 1
            push!(rows, placeholder_row(source_impl, T, 0.0, meson, idx))
        end
    end

    for T in T_ANISO
        for xi in xi_grid()
            for meson in MESONS
                idx += 1
                push!(rows, placeholder_row(source_impl, T, xi, meson, idx))
            end
        end
    end

    open(path, "w") do io
        println(io, join(MANDATORY_COLS, ','))
        for row in rows
            println(io, join(row, ','))
        end
    end

    println("Wrote template: " * path * " (" * string(length(rows)) * " rows)")
end

function main()
    fortran_path = joinpath(OUTDIR, "legacy_meson_scan_fortran_muB0_v1.csv")
    write_template(fortran_path, "fortran")
end

main()
