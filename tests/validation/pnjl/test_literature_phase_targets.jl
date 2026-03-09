using Test

const PHASE_VALIDATION_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PHASE_VALIDATION_DATA_PATH = joinpath(PHASE_VALIDATION_PROJECT_ROOT, "tests", "validation", "data", "pnjl_literature_phase_targets_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PHASE_VALIDATION_PROJECT_ROOT, "src", "models", "Models.jl"))
end
using .Models

function _load_phase_targets(path::String)
    isfile(path) || error("validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid validation target row: $line")
        push!(rows, (
            target_id=strip(cols[1]),
            observable=strip(cols[2]),
            expected_T_MeV=parse(Float64, strip(cols[3])),
            expected_mu_MeV=parse(Float64, strip(cols[4])),
            lower_T_MeV=parse(Float64, strip(cols[5])),
            upper_T_MeV=parse(Float64, strip(cols[6])),
            lower_mu_MeV=parse(Float64, strip(cols[7])),
            upper_mu_MeV=parse(Float64, strip(cols[8])),
            source=strip(cols[9]),
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

function _compute_reference_cep()
    tmp = mktempdir()
    result = Models.run_phase_pipeline(
        :PNJL;
        T_grid=[120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0],
        rho_grid=collect(0.1:0.1:3.0),
        xi=0.0,
        output_dir=tmp,
        profile=:regression,
        solver_backend=:legacy,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=40,
        compute_crossover=false,
        cep_strategy=:interpolate,
        promote_reference=false,
    )

    result.cep.found || error("CEP detection failed for validation")
    return result.cep
end

@testset "PNJL literature phase targets" begin
    targets = _load_phase_targets(PHASE_VALIDATION_DATA_PATH)
    cep_target = _target_by_id(targets, "gao_thesis_pnjl_cep")
    cep = _compute_reference_cep()

    @test cep.T_cep_MeV >= cep_target.lower_T_MeV
    @test cep.T_cep_MeV <= cep_target.upper_T_MeV
    @test cep.mu_cep_MeV >= cep_target.lower_mu_MeV
    @test cep.mu_cep_MeV <= cep_target.upper_mu_MeV
end