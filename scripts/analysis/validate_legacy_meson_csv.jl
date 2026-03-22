"""
Validate legacy meson CSV files against schema and numeric consistency.

Checks:
- required columns
- finite numeric columns policy by solver_status
- gap consistency: gap == mass - threshold
- mott flag consistency: abs(gap) <= tol

Usage:
  julia --project=. scripts/analysis/validate_legacy_meson_csv.jl
  julia --project=. scripts/analysis/validate_legacy_meson_csv.jl <csv1> [<csv2> ...]
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const DEFAULT_FILES = [
    joinpath(PROJECT_ROOT, "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB0_v1.csv"),
]

const REQUIRED = [
    "record_id", "source_impl", "T_MeV", "muB_MeV", "xi", "meson",
    "mass_MeV", "threshold_MeV", "gap_MeV", "mott_flag", "solver_status",
]

const MOTT_TOL = 1.0

_isnan_str(s::AbstractString) = lowercase(strip(s)) == "nan"

function _parse_float(s::AbstractString)
    t = strip(s)
    _isnan_str(t) && return NaN
    v = tryparse(Float64, t)
    v === nothing && error("invalid float value: $s")
    return v
end

function validate_file(path::String)
    isfile(path) || error("missing CSV: $path")
    lines = readlines(path)
    isempty(lines) && error("empty CSV: $path")

    header = split(strip(lines[1]), ',')
    length(header) >= length(REQUIRED) || error("header too short in $path")
    for c in REQUIRED
        c in header || error("missing required column '$c' in $path")
    end
    colidx = Dict(name => i for (i, name) in enumerate(header))

    total = 0
    warnings = 0
    for (k, line) in enumerate(lines[2:end])
        ln = k + 1
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) == length(header) || error("column count mismatch at $path:$ln")
        total += 1

        status = strip(cols[colidx["solver_status"]])
        mass = _parse_float(cols[colidx["mass_MeV"]])
        thr = _parse_float(cols[colidx["threshold_MeV"]])
        gap = _parse_float(cols[colidx["gap_MeV"]])
        mott = parse(Int, strip(cols[colidx["mott_flag"]]))

        if status in ("fail", "partial_threshold_only", "partial_threshold_nearestT")
            # allow NaN in mass/gap for partial/fail
            continue
        end

        if !(isfinite(mass) && isfinite(thr) && isfinite(gap))
            error("non-finite mass/threshold/gap for status=$status at $path:$ln")
        end

        if !isapprox(gap, mass - thr; atol=1e-8, rtol=1e-8)
            error("gap inconsistency at $path:$ln (gap=$(gap), mass-thr=$(mass-thr))")
        end

        expected = abs(gap) <= MOTT_TOL ? 1 : 0
        if mott != expected
            warnings += 1
        end
    end

    println("[validate] ", path)
    println("  rows: ", total)
    println("  mott-flag warnings: ", warnings)
    return nothing
end

function main(args::Vector{String})
    files = isempty(args) ? DEFAULT_FILES : args
    for f in files
        validate_file(f)
    end
    println("[validate] done")
end

main(ARGS)
