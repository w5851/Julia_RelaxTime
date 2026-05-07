using Test

if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

if !isdefined(Main, :_solve_relaxtime_literature_validation_equilibrium)
    include(joinpath(@__DIR__, "literature_validation_helpers.jl"))
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

using .MottReferenceMapping: estimate_mott_temperature

const RELAXTIME_MOTT_VALIDATION_DATA_PATH = validation_targets_path(
    "relaxtime",
    "literature",
    "mott",
    "relaxtime_mott_temperature_literature_targets_v1.csv",
)

const _MESON_WORKFLOW_MOD = Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

function _load_mott_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 8 || error("invalid mott validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            series=strip(cols[2]),
            meson=Symbol(strip(cols[3])),
            muB_MeV=parse(Float64, strip(cols[4])),
            xi=parse(Float64, strip(cols[5])),
            expected_T_mott_MeV=parse(Float64, strip(cols[6])),
            rtol=parse(Float64, strip(cols[7])),
            source=strip(cols[8]),
        ))
    end
    return rows
end

function _compute_relaxtime_mott_temperature_mev(muB_MeV::Float64, xi::Float64, meson::Symbol)
    Ts = collect(140.0:5.0:240.0)
    T_valid = Float64[]
    gaps = Float64[]
    for T_MeV in Ts
        T_fm = T_MeV / _HBARC_MEV_FM
        muq_fm = (muB_MeV / 3.0) / _HBARC_MEV_FM
        equilibrium = _solve_relaxtime_literature_validation_equilibrium(T_fm, muq_fm, xi)
        res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=xi,
            mesons=(meson,),
            p_num=8,
            t_num=4,
            seed_state=Vector(equilibrium.x_state),
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
        )
        gap_fm = Float64(res.meson_results[meson].gap)
        if isfinite(gap_fm)
            push!(T_valid, T_MeV)
            push!(gaps, gap_fm * _HBARC_MEV_FM)
        end
    end
    length(T_valid) >= 2 || error("insufficient points to estimate Mott temperature")
    est = estimate_mott_temperature(T_valid, gaps)
    return Float64(est.T_mott_MeV)
end

@testset "RelaxTime literature digitized Mott temperature targets" begin
    targets = _load_mott_targets(RELAXTIME_MOTT_VALIDATION_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        actual = _compute_relaxtime_mott_temperature_mev(row.muB_MeV, row.xi, row.meson)
        @test isfinite(actual)
        @test isapprox(actual, row.expected_T_mott_MeV; rtol=row.rtol, atol=0.0)
    end
end
