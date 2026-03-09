using Test

if !isdefined(Main, :validation_targets_path)
    Base.include(Main, joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

const PROJECT_ROOT = VALIDATION_PROJECT_ROOT
const DATA_PATH = validation_targets_path(
    "pnjl",
    "literature",
    "pnjl_literature_textual_targets_v1.csv",
)
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
using .Models

function _load_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 8 || error("invalid validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            observable=strip(cols[2]),
            variable=strip(cols[3]),
            mu_MeV=parse(Float64, strip(cols[4])),
            expected_MeV=parse(Float64, strip(cols[5])),
            lower_bound_MeV=parse(Float64, strip(cols[6])),
            upper_bound_MeV=parse(Float64, strip(cols[7])),
            source=strip(cols[8]),
        ))
    end
    return rows
end

function _target_by_id(rows, target_id::String)
    for row in rows
        row.target_id == target_id && return row
    end
    error("validation target not found: $target_id")
end

function _detect_mu0_crossover(variable::Symbol)
    result = Models.detect_crossover(
        0.0,
        (150.0 / HBARC_MEV_FM, 220.0 / HBARC_MEV_FM);
        method=:inflection,
        variable=variable,
        xi=0.0,
        n_scan=16,
        tol=1e-4,
        max_iter=20,
        p_num=8,
        t_num=4,
        model_kind=:PNJL,
        solver_backend=:models,
    )

    result.found || error("crossover detection failed for variable=$(variable)")
    result.T_crossover === nothing && error("missing crossover temperature for variable=$(variable)")
    return Float64(result.T_crossover) * HBARC_MEV_FM
end

@testset "PNJL literature textual targets" begin
    targets = _load_targets(DATA_PATH)
    chiral_target = _target_by_id(targets, "gao_thesis_mu0_chiral")
    deconf_target = _target_by_id(targets, "gao_thesis_mu0_deconf")

    chiral_T = _detect_mu0_crossover(:phi_u)
    deconf_T = _detect_mu0_crossover(:Phi)

    @test chiral_T >= chiral_target.lower_bound_MeV
    @test chiral_T <= chiral_target.upper_bound_MeV
    @test deconf_T >= deconf_target.lower_bound_MeV
    @test deconf_T <= deconf_target.upper_bound_MeV

    @test chiral_T > deconf_T
end