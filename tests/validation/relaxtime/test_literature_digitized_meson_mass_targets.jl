using Test

if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_SCRIPT = joinpath(_PROJECT_ROOT, "src", "models", "Models.jl")
const _CONSTANTS_SCRIPT = joinpath(_PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl")

isfile(_MODELS_SCRIPT) || error("Models script missing: $(_MODELS_SCRIPT)")
isfile(_CONSTANTS_SCRIPT) || error("Constants script missing: $(_CONSTANTS_SCRIPT)")

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_SCRIPT)
end

if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_SCRIPT)
end

const RELAXTIME_MESON_MASS_VALIDATION_DATA_PATH = validation_targets_path(
    "relaxtime",
    "literature",
    "meson",
    "relaxtime_meson_mass_literature_targets_v1.csv",
)

const _MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

function _load_meson_mass_targets(path::String)
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

function _compute_relaxtime_meson_mass_mev(T_MeV::Float64, muB_MeV::Float64, xi::Float64, meson::Symbol)
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

@testset "RelaxTime literature digitized meson mass targets" begin
    targets = _load_meson_mass_targets(RELAXTIME_MESON_MASS_VALIDATION_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        actual = _compute_relaxtime_meson_mass_mev(row.T_MeV, row.muB_MeV, row.xi, row.meson)
        @test isfinite(actual)
        @test isapprox(actual, row.expected_mass_MeV; rtol=row.rtol, atol=0.0)
    end
end
