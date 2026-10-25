using Test, JSON3, TOML
if !isdefined(Main,:PhaseGuidedCLI)
    include(joinpath(@__DIR__,"..","..","..","scripts","relaxtime","run_phase_guided_transport_scan.jl"))
end

@testset "sampling profile survives both scan CLI layers" begin
    mktempdir() do outdir
        args=["--mode","fixed-T-sparse-muB","--xi-list","0","--muB-list","0",
            "--T-list","200","--sampling-profile","thermal_weighted_v1","--outdir",outdir,"--dry-run"]
        opts=Main.PhaseGuidedCLI.parse_args(args)
        scan=Main._build_runtime_scan_opts(joinpath(outdir,"result.csv"),opts)
        @test (opts.p_num,opts.t_num)==(scan.p_num,scan.t_num)==(96,24)
        @test (scan.tau_p_nodes,scan.tau_angle_nodes,scan.tau_phi_nodes)==(48,36,36)
        @test scan.sigma_cache_policy===:thermal_weighted
        @test !hasproperty(scan,:tau_interpolation_mode)
        @test !scan.resume
        @test Main._build_runtime_scan_opts(scan.output,Main.PhaseGuidedCLI.parse_args(vcat(args,["--resume"]))).resume
        runtime=Main.build_scan_runtime(scan)
        @test runtime.model.params.thermal_p_max_inv_fm==25.0
        @test runtime.profile.source_path==opts.sampling_profile==scan.sampling_profile
        cfg=Main.build_effective_config(scan)
        @test cfg["sampling_profile"].background.thermal_p_max_inv_fm==25.0
        @test cfg["tr_p_max_fm"]===nothing
        @test cfg["terminal_integration_mode"]=="semi_infinite"
        @test !haskey(cfg,"sigma_grid_n")
        for extra in (["--p-num","12"],["--tau-angle-nodes","20"],["--sigma-cache-policy","default"],
            ["--sigma-grid-n","560"],["--phase-anchor-policy","reference_interpolation"])
            @test_throws ArgumentError Main.PhaseGuidedCLI.parse_args(vcat(args,extra))
        end
        @test_throws ArgumentError Main.parse_args(["--sampling-profile","thermal_weighted_v1","--tr-p-max","8"])
        reference=Main._load_runtime_phase_reference(opts)
        ctx=Main.ProvenanceMetadata.new_run_context("test_transport_sampling_profile_cli",args)
        plan=Main.run_phase_guided_scan(opts,ctx;phase_reference=reference)
        @test plan.total==1
        config_files=filter(p->endswith(p,"effective_config.json"),readdir(outdir;join=true))
        @test length(config_files)==1
        saved=JSON3.read(read(only(config_files),String))["options"]
        @test saved["sampling_profile"]["algorithm"]=="thermal_weighted"
        @test saved["sampling_profile"]["terminal_integral"]["p_nodes"]==32
    end
end

@testset "profile resume cannot mix configurations" begin
    mktempdir() do dir
        helpers=Main.TransportSamplingScanConfig
        opts=(output=joinpath(dir,"result.csv"),sampling_profile="thermal_weighted_v1",
            propagator_xi_policy=:match_thermo,compute_bulk=true,overwrite=false,resume=true)
        contract=helpers.prepare_sampling_output!(opts)
        @test isfile(contract)
        write(opts.output,"T_MeV,muB_MeV,xi\n200,0,0\n")
        original=read(opts.output)
        @test helpers.prepare_sampling_output!(opts)==contract
        @test_throws ArgumentError helpers.prepare_sampling_output!(merge(opts,(compute_bulk=false,)))
        @test_throws ArgumentError helpers.prepare_sampling_output!(merge(opts,(propagator_xi_policy=:isotropic,)))
        @test_throws ArgumentError helpers.prepare_sampling_output!(merge(opts,(resume=false,)))
        @test_throws ArgumentError helpers.prepare_sampling_output!(merge(opts,(sampling_profile=nothing,)))
        @test read(opts.output)==original
        rm(contract)
        @test_throws ArgumentError helpers.prepare_sampling_output!(opts)
        write(helpers.sampling_diagnostics_path(opts),"old diagnostic")
        @test isfile(helpers.prepare_sampling_output!(merge(opts,(overwrite=true,))))
        @test !isfile(helpers.sampling_diagnostics_path(opts))
    end
end

@testset "sampling failure never produces a valid transport point" begin
    mktempdir() do dir
        data=TOML.parsefile(Main.Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1").source_path)
        data["sampling"]["max_events"]=1
        path=joinpath(dir,"exhausted.toml")
        open(io->TOML.print(io,data),path,"w")
        args=["--output",joinpath(dir,"result.csv"),"--failed-points-output",joinpath(dir,"failed.csv"),
            "--sampling-profile",path,"--tmin","200","--tmax","200","--mubmin","0","--mubmax","0"]
        opts=Main.parse_args(args)
        reference=Main._load_runtime_phase_reference(opts)
        ctx=Main.ProvenanceMetadata.new_run_context("test_transport_sampling_failure",args)
        Main.run_scan(opts,ctx;phase_reference=reference)
        stats=Main.ProvenanceMetadata.csv_stats(opts.output)
        @test stats.points_total==0
        @test Main.ProvenanceMetadata.csv_stats(opts.failed_points_output).points_total==1
        manifest=JSON3.read(read(joinpath(dir,"run_manifest.json"),String))
        @test manifest["summary"]["success_count"]==0
        @test manifest["summary"]["error_count"]==1
        @test occursin("event budget",read(opts.failed_points_output,String))
    end
end
