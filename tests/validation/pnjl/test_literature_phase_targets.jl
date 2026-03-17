using Test

if !isdefined(Main, :_load_phase_targets)
    Base.include(Main, joinpath(@__DIR__, "pnjl_validation_helpers.jl"))
end

const PHASE_VALIDATION_PROJECT_ROOT = VALIDATION_PROJECT_ROOT
const PHASE_VALIDATION_DATA_PATH = validation_targets_path(
    "pnjl",
    "literature",
    "pnjl_literature_phase_targets_v1.csv",
)

using .Models

function _compute_reference_cep()
    tmp = mktempdir()
    result = Models.run_phase_pipeline(
        :PNJL;
        mode=:research,
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
