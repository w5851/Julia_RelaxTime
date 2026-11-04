using Test, JSON3, SHA

if !isdefined(Main,:PhaseGuidedCLI)
    include(joinpath(@__DIR__,"..","..","..","scripts","relaxtime","run_phase_guided_transport_scan.jl"))
end

@testset "phase guided transport dry-run artifacts" begin
    mktempdir() do outdir
        args=["--mode","fixed-T-sparse-muB","--outdir",outdir,"--case-name","smoke_case",
            "--xi-list","-0.2,0,0.2","--muB-list","150,400","--T-list","160,200",
            "--propagator-xi-policy","isotropic","--channel-diagnostics","--dry-run"]
        opts=Main.PhaseGuidedCLI.parse_args(args)
        ctx=Main.ProvenanceMetadata.new_run_context("scripts/relaxtime/run_phase_guided_transport_scan.jl",args)
        reference=Main._load_runtime_phase_reference(opts)
        plan=Main.run_phase_guided_scan(opts,ctx;phase_reference=reference)
        @test plan.total==12
        @test isfile(joinpath(outdir,"sampling_plan.csv"))
        @test isfile(joinpath(outdir,"README.md"))
        config=JSON3.read(read(joinpath(outdir,"effective_config.json"),String))["options"]
        @test config["mode"]=="mode_b_fixed_T_sparse_muB"
        @test config["propagator_xi_policy"]=="isotropic"
        @test config["sigma_cache_policy"]=="thermal_weighted"
        @test config["sampling_profile"]["tau_integral"]["p_nodes"]==48
        @test config["channel_diagnostics"]
        manifest=JSON3.read(read(joinpath(outdir,"run_manifest.json"),String))
        @test manifest["script"]=="scripts/relaxtime/run_phase_guided_transport_scan.jl"
        for entry in manifest["artifacts"]
            path=joinpath(outdir,basename(String(entry["path"])))
            @test isfile(path)
            @test entry["sha256"]==bytes2hex(sha256(read(path)))
        end
    end
end
