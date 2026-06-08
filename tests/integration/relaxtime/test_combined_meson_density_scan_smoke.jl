using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_combined_meson_density_scan.jl")
const CONTRACT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "combined_meson_density_scan_contract.jl")
const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const FIGDIR = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const OUTCSV = joinpath(OUTDIR, "combined_meson_density_scan.csv")
const OUTSVG = joinpath(FIGDIR, "combined_meson_density_scan.svg")
const PLOTMANIFEST = joinpath(FIGDIR, "plot_manifest.json")
const OUTREADME = joinpath(OUTDIR, "README.md")
const ASYM_PLUS_PROFILE = joinpath(PROJECT_ROOT, "config", "physics", "meson_chemical", "asymmetric_kplus_over_piplus_signed.toml")
const ASYM_MINUS_PROFILE = joinpath(PROJECT_ROOT, "config", "physics", "meson_chemical", "asymmetric_kminus_over_piminus_signed.toml")
const RUN_COMBINED_MESON_DENSITY_CLI = lowercase(get(ENV, "RUN_COMBINED_MESON_DENSITY_CLI_SMOKE", "0")) in ("1", "true", "yes")

if !isdefined(Main, :CombinedMesonDensityScanContract)
    include(CONTRACT)
end
const CMD_CONTRACT = Main.CombinedMesonDensityScanContract

@testset "combined meson density scan CLI contract" begin
    @test isfile(CONTRACT)
    tmu = CMD_CONTRACT.parse_args([
        "--output-dir", OUTDIR,
        "--figure-dir", FIGDIR,
        "--regimes", "stable,phase_shift_current",
        "--muq-values", "0,100",
        "--no-plot",
    ])
    @test tmu.path_strategy === :tmu
    @test tmu.output_dir == OUTDIR
    @test tmu.figure_dir == FIGDIR
    @test tmu.regimes == [:stable, :phase_shift_current]
    @test tmu.muq_values_MeV == [0.0, 100.0]
    @test tmu.plot == false

    asym = CMD_CONTRACT.parse_args([
        "--path", "trho_asymmetric",
        "--regimes", "stable,phase_shift_current",
        "--rho-values", "0.05,0.1",
        "--asym-ud-ratio-target", "0.876",
        "--asym-s-target", "0",
        "--no-trho-reverse-rho",
        "--meson-profile", "asymmetric_kplus_over_piplus_signed",
        "--density-policy", "x_min_cut",
        "--bose-x-min", "1e-6",
    ])
    @test asym.path_strategy === :trho_asymmetric
    @test asym.output_dir == CMD_CONTRACT.DEFAULT_TRHO_ASYMMETRIC_OUTPUT_DIR
    @test endswith(asym.figure_dir, joinpath("relaxtime", "meson_density", "combined_trho_asymmetric_smoke_scan"))
    @test asym.regimes == [:stable, :phase_shift_current]
    @test asym.rho_values == [0.05, 0.1]
    @test asym.trho_reverse_rho == false
    @test asym.asym_ud_ratio_target ≈ 0.876
    @test asym.meson_profile == "asymmetric_kplus_over_piplus_signed"
    @test asym.density_policy === :x_min_cut
    @test asym.bose_x_min ≈ 1e-6

    @test "constraint_mode" in CMD_CONTRACT.OUTPUT_COLUMNS
    @test "rho_u_over_rho_d" in CMD_CONTRACT.OUTPUT_COLUMNS
    @test "mu_pi_MeV" in CMD_CONTRACT.OUTPUT_COLUMNS
    @test "phase_display" in CMD_CONTRACT.OUTPUT_COLUMNS
    @test :phase_shift_gbu_reference in CMD_CONTRACT.DEFAULT_REGIMES
    @test_throws ArgumentError CMD_CONTRACT.parse_args(["--path", "tmu", "--rho-values", "0.05"])
    @test_throws ArgumentError CMD_CONTRACT.parse_args(["--path", "tmu", "--rhomin", "0", "--rhomax", "0.1", "--rhostep", "0.05"])
    @test_throws ArgumentError CMD_CONTRACT.parse_args(["--path", "trho_asymmetric", "--muq-values", "0,10"])
    @test_throws ArgumentError CMD_CONTRACT.parse_args(["--path", "trho_asymmetric", "--rhomin", "0.2", "--rhomax", "0.1", "--rhostep", "0.05"])
    @test_throws ArgumentError CMD_CONTRACT.parse_args(["--mumin", "0", "--mumax", "10", "--mustep", "0"])
    @test_throws ArgumentError CMD_CONTRACT.parse_args(["--bose-x-min", "-1e-6"])
end

if RUN_COMBINED_MESON_DENSITY_CLI
    @testset "combined meson density scan CLI smoke" begin
        isdir(OUTDIR) && rm(OUTDIR; recursive=true, force=true)
        isdir(FIGDIR) && rm(FIGDIR; recursive=true, force=true)
        mkpath(OUTDIR)

        cmd = `julia --project=. $SCRIPT --output-dir $OUTDIR --overwrite --tmin 150 --tmax 150 --tstep 10 --muq 0 --p-num 4 --t-num 2 --max-iter 6 --stable-q-nodes 6 --qmax 2 --q-nodes 3 --omega-max 2 --omega-nodes 3 --phase-display fold_0_pi`
        run(cmd)

        @test isfile(OUTCSV)
        @test isfile(OUTSVG)
        @test isfile(PLOTMANIFEST)
        @test isfile(OUTREADME)

        text = read(OUTCSV, String)
        @test occursin("# bridge: path_strategy x density_regime", text)
        @test occursin("path_strategy,path_point_index,T_MeV", text)
        @test occursin("# rho_values: not_applicable", text)
        @test occursin("# asym_ud_ratio_target: not_applicable", text)
        @test occursin("# asym_s_target: not_applicable", text)
        @test !occursin("# trho_seed_policy:", text)
        @test occursin("# bose_x_min: 0.0", text)
        @test occursin("# density_policy_scope: phase_shift_current,phase_shift_gbu_reference", text)
        @test occursin("phase_display", text)
        @test occursin("stable", text)
        @test occursin("strict_bw_stage1", text)
        @test occursin("phase_shift_current", text)
        @test occursin("phase_shift_gbu_reference", text)

        data_lines = [
            line for line in split(text, '\n')
            if !isempty(strip(line)) && !startswith(strip(line), "#") && !startswith(strip(line), "path_strategy,")
        ]
        @test length(data_lines) == 4

        summary = read(OUTREADME, String)
        @test occursin("Bridge-style composition", summary)
        @test occursin("figure directory", summary)
        @test occursin("strict_bw_omega_min=0.05", summary)
        @test occursin("Status Counts", summary)

        svg = read(OUTSVG, String)
        @test occursin("<svg", svg)
        @test occursin("Combined Meson Density Scan", svg)

        manifest = read(PLOTMANIFEST, String)
        @test occursin("combined_meson_density_plot_manifest_v1", manifest)
        @test occursin("combined_meson_density_scan.svg", manifest)
    end
else
    @info "Skipping slow combined meson density CLI smoke" env="RUN_COMBINED_MESON_DENSITY_CLI_SMOKE"
end

@testset "combined meson density FixedAsymmetricRho path contract" begin
    source = read(SCRIPT, String)
    @test occursin("combined_meson_density_scan_contract.jl", source)
    @test occursin("function _run_trho_asymmetric_scan", source)
    @test occursin("trho_reverse_rho", source)
    @test occursin("temperature_grouped_rho_continuity", source)
    @test occursin("Models.FixedAsymmetricRho", source)
    @test occursin("solve_meson_point_from_equilibrium", source)
    @test occursin("function _heatmap_axis_config", source)
    @test occursin("\"rho_target\", \"rho/rho0\", \"\"", source)
    @test occursin("Combined Meson Density Scan: T-rho heatmap", source)
    @test occursin("_write_svg_plot(plot_path, opts, rows)", source)

    @test isfile(ASYM_PLUS_PROFILE)
    @test isfile(ASYM_MINUS_PROFILE)
    plus = read(ASYM_PLUS_PROFILE, String)
    minus = read(ASYM_MINUS_PROFILE, String)
    @test occursin("mu_pi_rule = \"mu_u_minus_mu_d_signed\"", plus)
    @test occursin("mu_K_rule = \"mu_u_minus_mu_s_signed\"", plus)
    @test occursin("pi_label = \"pi_plus\"", plus)
    @test occursin("k_label = \"K_plus\"", plus)
    @test occursin("mu_pi_rule = \"mu_u_minus_mu_d_signed\"", minus)
    @test occursin("mu_K_rule = \"mu_u_minus_mu_s_signed\"", minus)
    @test occursin("pi_label = \"pi_minus\"", minus)
    @test occursin("k_label = \"K_minus\"", minus)
end
