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
        "--crossover-mu0-only",
        "--adaptive-xi",
        "--xi-refine-levels", "1",
        "--T-refine-levels", "1",
        "--rho-hybrid-endpoint-policy", "three_crossing_endpoint_local_v2",
    ])
    @test cfg.adaptive_xi
    @test cfg.xi_max_refine_level == 1
    @test cfg.temperature_max_refine_level == 1
    @test resolved_crossover_T_max_MeV(cfg) == 235.0
    @test cfg.crossover_mu_only_zero
    @test cfg.rho_hybrid_endpoint_policy == :three_crossing_endpoint_local_v2
    @test manifest_config_payload(cfg)["rho_hybrid_endpoint_policy"] == "three_crossing_endpoint_local_v2"

    inherited = DensePhaseReferenceConfig(T_min=10.0, T_max=240.0)
    @test resolved_crossover_T_max_MeV(inherited) == 240.0

    low_temperature_grid = inclusive_step_grid(1.0, 240.0, 5.0; axis="temperature")
    @test first(low_temperature_grid) == 1.0
    @test low_temperature_grid[end - 1] == 236.0
    @test last(low_temperature_grid) == 240.0
    @test length(low_temperature_grid) == 49
    @test inclusive_step_grid(60.0, 240.0, 5.0; axis="temperature") == collect(60.0:5.0:240.0)
    @test inclusive_step_grid(0.0, 1.0, 0.3; axis="rho") == [0.0, 0.3, 0.6, 0.9, 1.0]
    @test_throws ErrorException parse_args(["--T-min", "0", "--T-max", "240"])
    @test_throws ErrorException parse_args(["--T-min", "10", "--T-max", "5"])

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

    mktempdir() do root
        boundary_path = joinpath(root, "boundary.csv")
        write_boundary_csv(boundary_path, [(
            xi=0.0,
            T_MeV=100.0,
            mu_transition_MeV=300.0,
            rho_hadron=1.0,
            rho_quark=2.0,
            area_residual=5e-5,
            converged=true,
        )])
        header = first(readlines(boundary_path))
        @test occursin("area_residual,converged", header)
    end

    qualified = _dense_record_with_xi((axis="temperature", xi=-99.0, level=1), 0.25)
    @test qualified.xi == 0.25

    mktempdir() do root
        path = joinpath(root, "grid.csv")
        write_grid_convergence_csv(path, [(
            axis="temperature",
            xi=0.3,
            T_MeV=10.0,
            level=1,
            left=5.0,
            right=15.0,
            midpoint=10.0,
            position_error_MeV=0.1,
            density_error=0.01,
            maxwell_area=1e-4,
            response_rtol=0.05,
            converged=false,
            reason="valid,unknown,\"valid\"\nreview",
        )])
        @test occursin("\"valid,unknown,\"\"valid\"\"\nreview\"", read(path, String))
    end
end

end # module DensePhaseReferenceContractTests
