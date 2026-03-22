using Test

if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_SCRIPT = joinpath(_PROJECT_ROOT, "src", "models", "Models.jl")
const _CONSTANTS_SCRIPT = joinpath(_PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl")
const _MOTT_MAP_SCRIPT = joinpath(_PROJECT_ROOT, "scripts", "analysis", "mott_reference_mapping.jl")

isfile(_MODELS_SCRIPT) || error("Models script missing: $(_MODELS_SCRIPT)")
isfile(_CONSTANTS_SCRIPT) || error("Constants script missing: $(_CONSTANTS_SCRIPT)")
isfile(_MOTT_MAP_SCRIPT) || error("mott reference mapping script missing: $(_MOTT_MAP_SCRIPT)")

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_SCRIPT)
end

if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_SCRIPT)
end

if !isdefined(Main, :MottReferenceMapping)
    Base.include(Main, _MOTT_MAP_SCRIPT)
end

using .MottReferenceMapping: load_reference_table, validate_reference_schema

const _FORTRAN_FILES = (
    validation_targets_path("relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB0_v1.csv"),
    validation_targets_path("relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB600_v1.csv"),
)

const _MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

const _MESONS = (:pi, :K, :eta, :eta_prime, :sigma_pi, :sigma_K, :sigma, :sigma_prime)

const _RTOL_BY_MESON = Dict(
    :pi => 1e-2,
    :K => 1e-2,
    :eta => 4e-2,
    :eta_prime => 3e-1,
    :sigma_pi => 1e-2,
    :sigma_K => 1e-2,
    :sigma => 1e-2,
    :sigma_prime => 1e-2,
)

function _julia_symbol_for_legacy(meson::Symbol)
    if meson == :eta
        return :eta_prime
    elseif meson == :eta_prime
        return :eta
    elseif meson == :sigma
        return :sigma_prime
    elseif meson == :sigma_prime
        return :sigma
    end
    return meson
end

function _load_fortran_targets(files)
    out = Dict{Tuple{Float64,Float64,Symbol},Float64}()
    for path in files
        rows = validate_reference_schema(load_reference_table(path))
        for r in rows
            String(r[:source_impl]) == "fortran" || continue
            isapprox(r[:xi], 0.0; atol=1e-12) || continue
            meson = Symbol(r[:meson])
            meson in _MESONS || continue
            key = (Float64(r[:muB_MeV]), Float64(r[:T_MeV]), meson)
            out[key] = Float64(r[:mass_MeV])
        end
    end
    return out
end

function _julia_masses_at(T_MeV::Float64, muB_MeV::Float64)
    T_fm = T_MeV / _HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / _HBARC_MEV_FM
    res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        muq_fm;
        xi=0.0,
        mesons=_MESONS,
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=25,),
    )
    out = Dict{Symbol,Float64}()
    for m in _MESONS
        out[m] = Float64(res.meson_results[m].mass) * _HBARC_MEV_FM
    end
    return out
end

@testset "Julia vs Fortran all-8 meson numeric targets (muB=0 and 600, xi=0)" begin
    targets = _load_fortran_targets(_FORTRAN_FILES)

    expected_mus = (0.0, 600.0)
    expected_Ts = (120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0)

    for muB in expected_mus
        for T in expected_Ts
            julia_m = _julia_masses_at(T, muB)
            for meson in _MESONS
                key = (muB, T, meson)
                @test haskey(targets, key)
                legacy_mass = targets[key]
                actual = julia_m[_julia_symbol_for_legacy(meson)]
                @test isfinite(actual)
                @test isapprox(actual, legacy_mass; rtol=_RTOL_BY_MESON[meson], atol=0.0)
            end
        end
    end
end
