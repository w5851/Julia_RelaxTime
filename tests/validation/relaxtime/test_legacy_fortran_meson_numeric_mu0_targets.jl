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

const _FORTRAN_FILE = validation_targets_path(
    "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB0_v1.csv",
)
const _MESON_WORKFLOW_MOD = Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

function _fortran_mass_targets(rows, meson::String)
    out = Dict{Float64,Float64}()
    for r in rows
        String(r[:source_impl]) == "fortran" || continue
        String(r[:meson]) == meson || continue
        isapprox(r[:muB_MeV], 0.0; atol=1e-12) || continue
        isapprox(r[:xi], 0.0; atol=1e-12) || continue
        T_MeV = Float64(r[:T_MeV])
        mass_MeV = Float64(r[:mass_MeV])
        if !haskey(out, T_MeV)
            out[T_MeV] = mass_MeV
        end
    end
    return out
end

function _julia_meson_mass_mev(T_MeV::Float64, meson::Symbol)
    T_fm = T_MeV / _HBARC_MEV_FM
    res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        0.0;
        xi=0.0,
        mesons=(meson,),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=25,),
    )
    mr = res.meson_results[meson]
    return Float64(mr.mass) * _HBARC_MEV_FM
end

@testset "Julia vs legacy Fortran pi/K numeric targets (muB=0, xi=0)" begin
    rows = validate_reference_schema(load_reference_table(_FORTRAN_FILE))

    mesons = ("pi", "K")
    expected_Ts = [120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0]

    for meson_str in mesons
        targets = _fortran_mass_targets(rows, meson_str)
        @test sort(collect(keys(targets))) == expected_Ts

        meson = Symbol(meson_str)
        for T_MeV in expected_Ts
            legacy_mass = targets[T_MeV]
            julia_mass = _julia_meson_mass_mev(T_MeV, meson)

            @test isfinite(julia_mass)
            @test isapprox(julia_mass, legacy_mass; rtol=1e-2, atol=0.0)
        end
    end
end
