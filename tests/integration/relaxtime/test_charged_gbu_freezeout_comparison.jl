using Test, CSV, JSON3

isdefined(Main, :ChargedGBUFreezeoutComparison) || include(joinpath(@__DIR__, "..", "..", "..",
    "scripts", "analysis", "relaxtime", "run_charged_gbu_freezeout_comparison.jl"))
const FC_TEST = Main.ChargedGBUFreezeoutComparison

@testset "q0 chemical freezeout comparison inputs and dispatch" begin
    f = FC_TEST
    experimental_path = joinpath(f.ROOT, "data", "inputs", "meson_freezeout_experiment_20261017", "experimental_ratios.csv")
    @test all(r.sqrt_s_NN_GeV in f.DEFAULT_ENERGIES for r in CSV.File(experimental_path))
    pts = f.comparison_points(["cleymans2006", "andronic2018"], [7.7, 200.])
    @test [p.profile_name for p in pts] == ["cleymans2006", "cleymans2006", "andronic2018", "andronic2018"]
    @test [p.point.sqrt_s_NN_GeV for p in pts] == [200., 7.7, 200., 7.7]
    @test pts[3].point.T_MeV ≈ 158.384 atol=.001
    @test f.numerical_settings("screening") == f.Scan.screening_settings(route=f.ROUTE)
    @test f.numerical_settings("standard").q_nodes == 16
    @test f.numerical_settings("fine").q_nodes == 32
    @test f.numerical_settings("fine").omega_lower_inv_fm == 1e-5
    @test_throws ArgumentError f.numerical_settings("production")
    @test_throws ArgumentError f.comparison_points(["andronic2018", "andronic2018"], [200.])
    for invalid in ([], [7.7, 7.7], [NaN], [0.], [-1.])
        @test_throws ArgumentError f.comparison_points(["andronic2018"], invalid)
    end
    @test_throws ArgumentError f.comparison_points(["crossover"], [200.])
    @test f.parse_args(["--help"]) === nothing
    args = f.parse_args(["--output", "unused", "--profiles", "andronic2018", "--sqrts-list", "7.7,200",
        "--resolution", "fine", "--resume"])
    @test args[:profiles] == ["andronic2018"] && args[:energies] == [7.7, 200.]
    @test args[:resolution] == "fine" && args[:resume]
    @test args[:parameter_policy] == "published_central"
    band_args = f.parse_args(["--output", "unused", "--parameter-policy", "literature_envelope_exploration"])
    @test band_args[:parameter_policy] == "literature_envelope_exploration"
    @test_throws ArgumentError f.parse_args(["--output", "unused", "--parameter-policy", "confidence_band"])
    @test_throws ArgumentError f.parse_args(["--output"])
end

@testset "relative profile directory retains its source and snapshots inside output" begin
    f = FC_TEST
    mktempdir(f.ROOT) do dir
        profile_dir = joinpath(dir, "profiles")
        mkpath(profile_dir)
        source = joinpath(profile_dir, "andronic2018.toml")
        cp(joinpath(f.FP.freezeout_profile_dir(), "andronic2018.toml"), source)
        original_hash = f.Workflow.hashfile(source)
        relative_dir = relpath(profile_dir, pwd())
        hashes = f.source_hashes(f.Workflow.DEFAULT_CONFIG, ["andronic2018"], relative_dir)
        @test hashes == f.source_hashes(f.Workflow.DEFAULT_CONFIG, ["andronic2018"], profile_dir)
        @test all(!isabspath(p) for p in keys(hashes))
        out = joinpath(dir, "run")
        result = f.run_scan(output=out, profiles=["andronic2018"], energies=[200.],
            resolution="screening", profile_dir=relative_dir, parameter_policy="literature_envelope_exploration")
        @test result.successful_points == 1
        @test result.parameter_policy == "literature_envelope_exploration"
        @test f.Workflow.hashfile(source) == original_hash
        key = replace(relpath(source, f.ROOT), '\\'=>'/')
        @test f.Workflow.hashfile(joinpath(out, "source_snapshot", key)) == original_hash
    end
end

@testset "real q0 run, frozen source, CSV, resume and failure semantics" begin
    f = FC_TEST
    mktempdir() do dir
        out = joinpath(dir, "run")
        result = f.run_scan(output=out, energies=[200.], resolution="screening")
        @test result.point_count == result.successful_points == 2
        @test result.density_route == "q0_lambda_reference"
        @test result.parameter_policy == "published_central"
        @test !result.production_authorized && result.diagnostic_only
        rows = collect(CSV.File(joinpath(out, "freezeout_ratios.csv")))
        @test length(rows) == 2
        for row in rows
            @test row.plus_passed && row.background_residual <= 1e-7
            @test row.Kplus_over_pi_plus ≈ row.K_plus_density_inv_fm3/row.pi_plus_density_inv_fm3
            @test row.endpoint_policy == "finite_window"
            path = joinpath(out, "points", row.profile_name, "energy_001.json")
            point = f.Workflow.readjson(path)
            # Reevaluate one channel at its retained background: verifies the
            # new path runner uses the existing q0 reference, not the legacy BU.
            bg = f.Workflow.restore_background(point.background)
            direct = f.Scan._channel_record(bg, :K_plus, result.settings; route=f.ROUTE)
            @test direct.density_inv_fm3 ≈ row.K_plus_density_inv_fm3 rtol=1e-12
        end
        @test rows[1].T_MeV != rows[2].T_MeV
        frozen = copy(result.outputs)
        resumed = f.run_scan(output=out, energies=[200.], resolution="screening", resume=true)
        @test resumed.outputs == frozen
        @test_throws ArgumentError f.run_scan(output=out, energies=[200.])
        @test_throws ArgumentError f.run_scan(output=out, energies=[130.], resolution="screening", resume=true)
        @test_throws ArgumentError f.run_scan(output=out, energies=[200.], resolution="screening",
            resume=true, parameter_policy="literature_envelope_exploration")
        state = f.Workflow.readjson(joinpath(out, "run.json"))
        @test state.background_contract.nS_inv_fm3 == 0.
        @test state.parameter_policy == result.parameter_policy
        @test state.background_contract.nQ_over_nB == .4
        @test !state.background_contract.meson_feedback
        @test state.coordinate_contract.bose_frequency == "external omega"
        for (path, hash) in pairs(state.source_hashes)
            @test f.Workflow.hashfile(joinpath(out, "source_snapshot", String(path))) == hash
        end
        sample = f.Workflow.readjson(joinpath(out, "points", "andronic2018", "energy_001.json"))
        failed = merge(Dict(String(k)=>v for (k,v) in pairs(sample)),
            Dict("status"=>"gate_failed", "ratios"=>Dict("plus_passed"=>false, "Kplus_over_pi_plus"=>nothing)))
        parsed = JSON3.read(JSON3.write(failed))
        summary = f.summarize(parsed)
        @test !summary.plus_passed && isnan(summary.Kplus_over_pi_plus)
    end
end
