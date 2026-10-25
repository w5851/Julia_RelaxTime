using Test

const _SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_gap_transport_scan.jl")
if !isdefined(Main, :ScanOptions)
    include(_SCAN_SCRIPT)
end

@testset "temperature scan aliases preserve a nonuniform muB list" begin
    include(joinpath(dirname(_SCAN_SCRIPT),"scan_relaxation_times_vs_T.jl"))
    parsed=Main.RelaxationTemperatureScan.parse_args(["--out","tau.csv","--mub-list","0,450,900",
        "--tmin","150","--tmax","160","--tstep","10"])
    opts=parsed.scan_options
    @test opts.output=="tau.csv"
    @test !opts.compute_bulk
    plan=Main.build_scan_plan(opts;muB_values=parsed.muB_values)
    @test plan.muB_values==[0.0,450.0,900.0]
    @test plan.total==6
    @test_throws ArgumentError Main.build_scan_plan(opts;muB_values=[NaN])
    @test_throws ArgumentError Main.RelaxationTemperatureScan.parse_args(["--sigma-grid-n","60"])
end

@testset "run_gap_transport_scan solver entry" begin
    for name in (:solve_models_equilibrium, :solve_equilibrium_with_diagnostics,
        :build_scan_runtime, :execute_gap_transport_scan_point!, :build_scan_plan, :execute_scan_plan!)
        @test isdefined(Main, name)
    end
    @test :transport_point_api in Main.Models.list_precompile_capabilities()
    @test :transport_point_api in Main.Models.list_precompile_profile(:scan)
end

@testset "fixed-temperature aliases use the same sampling profile" begin
    include(joinpath(dirname(_SCAN_SCRIPT),"scan_transport_vs_xi_T150_muB800.jl"))
    opts=Main.FixedTemperatureTransportScan.parse_args(["--out","fixed.csv","--T-MeV","190",
        "--muB-MeV","450","--xi-list","-0.5,0,0.37"])
    @test (opts.tmin_mev,opts.tmax_mev)==(190.0,190.0)
    @test (opts.mubmin_mev,opts.mubmax_mev)==(450.0,450.0)
    @test opts.xi_values==[-0.5,0.0,0.37]
    @test opts.output=="fixed.csv"
    @test !opts.compute_bulk
    @test opts.sigma_cache_policy===:thermal_weighted
    @test Main.FixedTemperatureTransportScan.parse_args(["--xi0-only"]).xi_values==[0.0]
    @test_throws ArgumentError Main.FixedTemperatureTransportScan.parse_args(["--sigma-cache-dir","old"])
    @test_throws ArgumentError Main.FixedTemperatureTransportScan.parse_args(["--xi-min","-1"])
    @test_throws ArgumentError Main.GapTransportScanCLI.parse_args(["--max-iter","0"])
end

@testset "preloaded Models identity remains stable" begin
    original_models = Main.Models
    solver = Main.Models.NLsolveGapSolver(method=:trust_region, jacobian=:finite)
    include(_SCAN_SCRIPT)
    @test Main.Models === original_models
    @test solver isa Main.Models.AbstractGapSolver
end

@testset "gap transport scan provenance shell" begin
    opts = Main.parse_args(["--output", "tmp.csv", "--channel-diagnostics-output", "diag.csv",
        "--failed-points-output", "failed.csv", "--provenance-dir", "prov",
        "--xi-list", "-0.5,0", "--overwrite", "--no-resume", "--no-compute-bulk"])
    effective = Main.build_effective_config(opts)
    @test effective["output"] == "tmp.csv"
    @test effective["sigma_cache_policy"] == "thermal_weighted"
    @test effective["integration_mode"] == "semi_infinite"
    @test effective["sampling_profile"].id == "thermal_weighted_v1"
    @test !haskey(effective,"tau_threshold_subtraction")
    summary = Main.build_summary(3, 1, 2)
    @test summary["points_total"] == 4
    @test summary["success_count"] == 3
    @test summary["skipped_count"] == 2
    @test Main.collect_artifacts(opts) == ["tmp.csv", "diag.csv", "failed.csv"]
end

@testset "prefer Models.solve over solve_constraint" begin
    # Preserve this solver-only check's historical background discretization.
    opts = (p_num=12,t_num=6)
    eq = Main.solve_models_equilibrium(150 / Main.Constants_PNJL.ħc_MeV_fm,
        0.0,-0.5,Main.TransportWorkflow.PNJL.HADRON_SEED_5,opts)
    @test eq !== nothing
    @test eq.solver_backend == :models
end

@testset "failed point sidecar" begin
    opts = Main.parse_args(["--output", "tmp.csv", "--failed-points-output", "tmp_failed.csv",
        "--tmin", "150", "--tmax", "150", "--tstep", "1",
        "--mubmin", "0", "--mubmax", "0", "--mubstep", "1", "--xi-list", "-0.5"])
    @test opts.failed_points_output == "tmp_failed.csv"
    io = IOBuffer()
    diag = (seed_source="seed",phase_prev=:unknown,phase_curr=:hadron)
    Main.write_failed_point_row!(io,150.0,0.0,-0.5,diag,ArgumentError("boom"))
    line = String(take!(io))
    @test occursin("150.0,0.0,-0.5",line)
    @test occursin("seed",line)
    @test occursin("ArgumentError",line)
    @test occursin("boom",line)
end

@testset "scan runtime binds the qualified sampling profile" begin
    opts = Main.parse_args(["--output","tmp.csv","--propagator-xi-policy","isotropic"])
    runtime = Main.build_scan_runtime(opts)
    @test opts.propagator_xi_policy === :isotropic
    @test opts.sigma_cache_policy === :thermal_weighted
    @test runtime.tau_kwargs.p_nodes == 48
    @test runtime.tau_kwargs.sampling_config.rtol == 1e-3
    @test runtime.profile.id == "thermal_weighted_v1"
    @test runtime.model.params.thermal_p_max_inv_fm == 25
    @test_throws ArgumentError Main.parse_args(["--sigma-grid-n","17"])
    @test_throws ArgumentError Main.parse_args(["--sigma-cache-policy","validated_anchored"])
end

@testset "phase-reference runtime switch defaults to accepted with explicit strict" begin
    default_opts = Main.parse_args(String[])
    @test default_opts.phase_reference_mode === :runtime
    @test default_opts.phase_reference_layer === :accepted
    strict_opts = Main.parse_args(["--phase-reference-mode","strict","--phase-reference-layer","strict"])
    @test strict_opts.phase_reference_mode === :strict
    @test strict_opts.phase_reference_layer === :strict
    @test_throws ErrorException Main.parse_args(["--phase-reference-mode","legacy"])
    diagnostic_opts = Main.parse_args(["--phase-reference-root","candidate","--phase-reference-mode","diagnostic"])
    @test diagnostic_opts.phase_reference_mode === :diagnostic
end
