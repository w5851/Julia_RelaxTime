using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_combined_meson_density_scan.jl")
const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const FIGDIR = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const OUTCSV = joinpath(OUTDIR, "combined_meson_density_scan.csv")
const OUTSVG = joinpath(FIGDIR, "combined_meson_density_scan.svg")
const PLOTMANIFEST = joinpath(FIGDIR, "plot_manifest.json")
const OUTREADME = joinpath(OUTDIR, "README.md")
const ASYM_PLUS_PROFILE = joinpath(PROJECT_ROOT, "config", "physics", "meson_chemical", "asymmetric_kplus_over_piplus_signed.toml")
const ASYM_MINUS_PROFILE = joinpath(PROJECT_ROOT, "config", "physics", "meson_chemical", "asymmetric_kminus_over_piminus_signed.toml")
const RUN_COMBINED_MESON_DENSITY_CLI = lowercase(get(ENV, "RUN_COMBINED_MESON_DENSITY_CLI_SMOKE", "0")) in ("1", "true", "yes")

@testset "combined meson density scan script contract" begin
    source = read(SCRIPT, String)
    @test occursin("const DEFAULT_REGIMES", source)
    @test occursin(":stable", source)
    @test occursin(":strict_bw_stage1", source)
    @test occursin(":phase_shift_current", source)
    @test occursin(":phase_shift_gbu_reference", source)
    @test occursin("\"path_strategy\", \"path_point_index\"", source)
    @test occursin("\"phase_display\"", source)
    @test occursin("function _write_csv", source)
    @test occursin("function _write_svg_plot", source)
    @test occursin("function _write_plot_manifest", source)
    @test occursin("Bridge-style composition", source)
    @test occursin("--figure-dir", source)
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
    @test occursin("--path <tmu|trho_asymmetric>", source)
    @test occursin("DEFAULT_TRHO_ASYMMETRIC_OUTPUT_DIR", source)
    @test occursin("function _run_trho_asymmetric_scan", source)
    @test occursin("Models.FixedAsymmetricRho", source)
    @test occursin("solve_meson_point_from_equilibrium", source)
    @test occursin("\"constraint_mode\"", source)
    @test occursin("\"rho_u_over_rho_d\"", source)
    @test occursin("\"mu_u_MeV\"", source)
    @test occursin("\"mu_pi_MeV\"", source)

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
