using Test

if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MOTT_MAP_SCRIPT = joinpath(_PROJECT_ROOT, "scripts", "analysis", "mott_reference_mapping.jl")

isfile(_MOTT_MAP_SCRIPT) || error("mott reference mapping script missing: $(_MOTT_MAP_SCRIPT)")

if !isdefined(Main, :MottReferenceMapping)
    Base.include(Main, _MOTT_MAP_SCRIPT)
end

using .MottReferenceMapping: load_reference_table, validate_reference_schema

const _FORTRAN_MUB600_FILE = validation_targets_path(
    "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB600_v1.csv",
)
const _LITERATURE_FILE = validation_targets_path(
    "relaxtime", "literature", "meson", "relaxtime_meson_mass_literature_targets_v1.csv",
)

function _fortran_k_targets_muB600(rows)
    out = Dict{Float64,Float64}()
    for r in rows
        String(r[:source_impl]) == "fortran" || continue
        String(r[:meson]) == "K" || continue
        isapprox(r[:muB_MeV], 600.0; atol=1e-12) || continue
        isapprox(r[:xi], 0.0; atol=1e-12) || continue
        T_MeV = Float64(r[:T_MeV])
        out[T_MeV] = Float64(r[:mass_MeV])
    end
    return out
end

function _literature_k_targets_muB600(rows)
    out = Dict{Float64,Float64}()
    for r in rows
        String(r[:meson]) == "K" || continue
        isapprox(r[:muB_MeV], 600.0; atol=1e-12) || continue
        isapprox(r[:xi], 0.0; atol=1e-12) || continue
        T_MeV = Float64(r[:T_MeV])
        out[T_MeV] = Float64(r[:target_mass_MeV])
    end
    return out
end

function _load_literature_meson_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid meson mass validation target row: $line")
        push!(rows, (
            meson=strip(cols[3]),
            muB_MeV=parse(Float64, strip(cols[4])),
            xi=parse(Float64, strip(cols[5])),
            T_MeV=parse(Float64, strip(cols[6])),
            target_mass_MeV=parse(Float64, strip(cols[7])),
        ))
    end
    return rows
end

@testset "Legacy Fortran K targets available at muB=600, xi=0" begin
    fortran_rows = validate_reference_schema(load_reference_table(_FORTRAN_MUB600_FILE))
    literature_rows = _load_literature_meson_targets(_LITERATURE_FILE)

    fortran_k = _fortran_k_targets_muB600(fortran_rows)
    literature_k = _literature_k_targets_muB600(literature_rows)

    expected_Ts = [160.0, 180.0, 200.0, 220.0, 240.0, 260.0]
    @test sort(collect(keys(literature_k))) == expected_Ts
    @test all(T -> haskey(fortran_k, T), expected_Ts)

    for T_MeV in expected_Ts
        @test haskey(fortran_k, T_MeV)
        @test isfinite(fortran_k[T_MeV])
        @test fortran_k[T_MeV] > 0.0
    end
end
