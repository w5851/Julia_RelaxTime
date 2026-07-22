using Test
using JSON3

const XI_PLAN_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(XI_PLAN_PROJECT_ROOT, "src", "models", "Models.jl"))
end

module DenseReferenceXiRefinementPlanTests

using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "pnjl", "plan_dense_reference_xi_refinement.jl"))

function _write_reference(root::String, tag::String)
    mkpath(root)
    xis = [-0.1, 0.0, 0.1]
    open(joinpath(root, "boundary_$(tag).csv"), "w") do io
        println(io, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark,area_residual,converged,curve_parameter,plot_order_key")
        for xi in xis
            nonlinear = xi == 0.0 ? 0.2 : 0.0
            println(io, "$(xi),100,$(300 + 2xi + nonlinear),1,2,0.00005,true,100,100")
        end
    end
    open(joinpath(root, "cep_$(tag).csv"), "w") do io
        println(io, "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV,uncertainty_T_MeV,T_bracket_low_MeV,T_bracket_high_MeV,bracket_width_T_MeV")
        for xi in xis
            println(io, "$(xi),$(130 + 2xi),$(295 + 2xi),$(3 * (295 + 2xi)),0.05,129.95,130.05,0.1")
        end
    end
    open(joinpath(root, "spinodals_$(tag).csv"), "w") do io
        println(io, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark,curve_parameter,plot_order_key")
        for xi in xis
            println(io, "$(xi),100,$(310 + 2xi),$(290 + 2xi),0.8,2.2,100,100")
        end
    end
    open(joinpath(root, "crossover_$(tag).csv"), "w") do io
        println(io, "xi,mu_MeV,T_crossover_MeV,rho,method,converged,derivative,variable,curve_parameter,plot_order_key")
        for xi in xis
            println(io, "$(xi),0,$(180 + 2xi),0.3,peak,true,4,phi_u,0,0")
        end
    end
    open(joinpath(root, "phase_grid_convergence_$(tag).csv"), "w") do io
        println(io, "axis,xi,T_MeV,level,left,right,midpoint,position_error_MeV,density_error,maxwell_area,response_rtol,converged,reason")
    end
    manifest = Dict("config" => Dict("xi_values" => xis))
    write(joinpath(root, "phase_reference_$(tag)_manifest.json"), JSON3.write(manifest))
end

@testset "staged dense-reference xi refinement plan" begin
    mktempdir() do root
        tag = "test"
        _write_reference(root, tag)
        intervals_path = joinpath(root, "intervals.json")
        write(intervals_path, JSON3.write([["-0.1", "0.1"]]))
        cfg = DenseXiRefinementPlanConfig(
            reference_root=root,
            tag=tag,
            intervals_path=intervals_path,
            level=1,
            max_refine_level=2,
            position_tol_MeV=0.1,
            density_tol=0.01,
            maxwell_area_tol=1e-4,
            response_rtol=0.05,
            matrix_output=joinpath(root, "matrix.json"),
            intervals_output=joinpath(root, "next_intervals.json"),
            records_output=joinpath(root, "records.csv"),
            summary_output=joinpath(root, "summary.json"),
        )

        summary = build_refinement_plan(cfg)
        @test summary["evaluated_interval_count"] == 1
        @test summary["unconverged_interval_count"] == 1
        @test summary["next_interval_count"] == 2

        matrix = JSON3.read(read(cfg.matrix_output, String))
        @test length(matrix.include) == 2
        @test sort(Float64[parse(Float64, entry.xi) for entry in matrix.include]) == [-0.05, 0.05]
        records = read(cfg.records_output, String)
        @test occursin("xi,0.0,,1,-0.1,0.1,0.0", records)
        @test occursin("interpolation_tolerance_exceeded", records)

        _write_records(cfg.records_output, NamedTuple[(
            axis="xi",
            xi=0.0,
            T_MeV=nothing,
            level=1,
            left=-0.1,
            right=0.1,
            midpoint=0.0,
            position_error_MeV=0.2,
            density_error=0.01,
            maxwell_area=1e-4,
            response_rtol=0.05,
            converged=false,
            reason="valid,unknown,\"valid\"\nreview",
        )])
        @test occursin("\"valid,unknown,\"\"valid\"\"\nreview\"", read(cfg.records_output, String))
    end
end

end # module DenseReferenceXiRefinementPlanTests
