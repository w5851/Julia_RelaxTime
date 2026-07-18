using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

module DensePhaseReferenceContractTests

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DENSE_REFERENCE_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "build_dense_phase_reference.jl")

include(DENSE_REFERENCE_SCRIPT)

function _dense_synthetic_result(xi::Float64; nonlinear::Float64=0.0)
    offset = 2.0 * xi + nonlinear
    return Models.PhasePipelineResult(
        xi=xi,
        cep=Models.CEPResult(found=true, T_cep_MeV=130.0 + offset, mu_cep_MeV=295.0 + offset),
        first_order_boundary=[(
            T_MeV=100.0,
            mu_transition_MeV=300.0 + offset,
            rho_hadron=1.0,
            rho_quark=2.0,
            area_residual=5e-5,
            converged=true,
        )],
        spinodal=[(
            T_MeV=100.0,
            mu_spinodal_hadron_MeV=310.0 + offset,
            mu_spinodal_quark_MeV=290.0 + offset,
            rho_spinodal_hadron=0.8,
            rho_spinodal_quark=2.2,
        )],
    )
end

@testset "Dense phase reference contract" begin
    cfg = parse_args([
        "--xi-list", "-0.1,0.1",
        "--T-min", "10",
        "--T-max", "240",
        "--crossover-T-max", "235",
        "--adaptive-xi",
        "--xi-refine-levels", "1",
        "--T-refine-levels", "1",
    ])
    @test cfg.adaptive_xi
    @test cfg.xi_max_refine_level == 1
    @test cfg.temperature_max_refine_level == 1
    @test resolved_crossover_T_max_MeV(cfg) == 235.0

    inherited = DensePhaseReferenceConfig(T_min=10.0, T_max=240.0)
    @test resolved_crossover_T_max_MeV(inherited) == 240.0

    cache = Dict(
        -0.1 => _dense_synthetic_result(-0.1),
        0.1 => _dense_synthetic_result(0.1),
    )
    actual_evaluations = Float64[]
    evaluator = function (xi::Float64)
        key = round(xi; digits=12)
        if !haskey(cache, key)
            push!(actual_evaluations, key)
            cache[key] = _dense_synthetic_result(key; nonlinear=(key == 0.0 ? 0.2 : 0.0))
        end
        return cache[key]
    end
    resolved, audit = _adaptive_xi_refinement!([-0.1, 0.1], evaluator, cfg)
    @test resolved == [-0.1, 0.0, 0.1]
    @test actual_evaluations == [0.0]
    @test length(audit) == 1
    @test audit[1].axis == "xi"
    @test !audit[1].converged

    @test_throws ErrorException parse_args([
        "--xi-list", "-0.1,0.1",
        "--adaptive-xi",
        "--crossover-only",
    ])
end

end # module DensePhaseReferenceContractTests
