using Test
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "pnjl_bifurcation_spike", "run_production_comparison.jl"))

@testset "Production comparison serializes failed points" begin
    P = PALCProductionComparison
    payload = (; rows=[(; value=NaN, residual=Inf, converged=false)],
        nested=Dict("missing"=>-Inf), count=3, route=:rho)
    safe = P._json_finite(payload)
    parsed = P.JSON.parse(P.JSON.json(safe))
    @test parsed["rows"][1]["value"] === nothing
    @test parsed["rows"][1]["residual"] === nothing
    @test parsed["nested"]["missing"] === nothing
    @test parsed["count"] == 3
    @test parsed["route"] == "rho"
end

@testset "Phase CLI propagates continuation from configuration and arguments" begin
    if !isdefined(Main, :PhaseCliSupport)
        include(joinpath(@__DIR__, "..", "..", "..", "scripts", "pnjl", "phase_cli_support.jl"))
    end
    P = PhaseCliSupport
    @test P.parse_args(["--continuation=palc"]).continuation === :palc
    @test P.parse_args(["--continuation=legacy"]).continuation === :legacy
    @test_throws ArgumentError P.parse_args(["--continuation=unknown"])
    mktempdir() do directory
        config = joinpath(directory, "phase.toml")
        write(config, "[phase_pipeline]\ncontinuation = \"single_pass\"\n")
        @test P.parse_args(["--config=$config"]).continuation === :single_pass
        @test P.parse_args(["--config=$config", "--continuation=palc"]).continuation === :palc
        received = Ref{Any}()
        fake_models = (; RhoSupportConfig=Models.RhoSupportConfig,
            RhoHybridVerificationConfig=Models.RhoHybridVerificationConfig,
            run_phase_pipeline=(args...; kwargs...) -> begin
                received[] = (; kwargs...)
                throw(InterruptException())
            end)
        @test_throws InterruptException P.main(fake_models, directory, ["--continuation=palc"])
        @test received[].continuation === :palc
    end
end

@testset "Complete scan resume preserves CSV and skips callbacks" begin
    P = PALCProductionComparison
    for route in (:asym, :rho, :tmu), strategy in (:legacy, :single_pass, :palc)
        mktempdir() do directory
            file = joinpath(directory, "scan.csv")
            open(file, "w") do io
                println(io, "T,coordinate,xi")
                for coordinate in 1.0:5.0
                    println(io, "120.0,$coordinate,0.0")
                end
            end
            result = P.scan_contract_check(route, strategy, [120.0], [0.0], collect(1.0:5.0), file)
            @test result.skipped == 5
            @test result.byte_preserving
        end
    end
end

@testset "Density paths reject mismatched model parameters before solving" begin
    original = Models.create_model(:PNJL)
    modified = Models.PNJLModel(merge(original.consts, (; thermal_p_max_inv_fm=11.0)))
    @test_throws ArgumentError Models._path_options(modified, Models.FixedRhoPath(1.0, [0.5]))
    @test Models._path_options(original, Models.FixedRhoPath(1.0, [0.5])).residual_tol == 1e-6
    @test Models._path_options(modified, Models.FixedMuPath(1.0, [0.5])).residual_tol == 1e-6
end
