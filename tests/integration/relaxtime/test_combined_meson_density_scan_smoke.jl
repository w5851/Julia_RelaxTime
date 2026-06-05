using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_combined_meson_density_scan.jl")
const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const FIGDIR = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const OUTCSV = joinpath(OUTDIR, "combined_meson_density_scan.csv")
const OUTSVG = joinpath(FIGDIR, "combined_meson_density_scan.svg")
const PLOTMANIFEST = joinpath(FIGDIR, "plot_manifest.json")
const OUTREADME = joinpath(OUTDIR, "README.md")
const ASYM_OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "test_outputs", "combined_meson_density_asymmetric")
const ASYM_FIGDIR = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "scan", "test_outputs", "combined_meson_density_asymmetric")
const ASYM_OUTCSV = joinpath(ASYM_OUTDIR, "combined_meson_density_scan.csv")

@testset "combined meson density scan smoke" begin
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

@testset "combined meson density FixedAsymmetricRho path smoke" begin
    isdir(ASYM_OUTDIR) && rm(ASYM_OUTDIR; recursive=true, force=true)
    isdir(ASYM_FIGDIR) && rm(ASYM_FIGDIR; recursive=true, force=true)
    mkpath(ASYM_OUTDIR)

    cmd = `julia --project=. $SCRIPT --path trho_asymmetric --output-dir $ASYM_OUTDIR --figure-dir $ASYM_FIGDIR --overwrite --regimes stable,phase_shift_current --tmin 100 --tmax 100 --tstep 10 --rho-values 0.05 --asym-ud-ratio-target 0.876 --asym-s-target 0 --meson-profile asymmetric_kplus_over_piplus_signed --p-num 4 --t-num 2 --max-iter 8 --stable-q-nodes 6 --qmax 2 --q-nodes 3 --omega-max 2 --omega-nodes 3 --phase-display fold_0_pi`
    run(cmd)

    @test isfile(ASYM_OUTCSV)
    text = read(ASYM_OUTCSV, String)
    @test occursin("# path_strategy: trho_asymmetric", text)
    @test occursin("constraint_mode", text)
    @test occursin("rho_u_over_rho_d", text)
    @test occursin("mu_u_MeV", text)
    @test occursin("mu_pi_MeV", text)
    @test occursin("stable", text)
    @test occursin("phase_shift_current", text)

    lines = [
        split(line, ',')
        for line in split(text, '\n')
        if !isempty(strip(line)) && !startswith(strip(line), "#")
    ]
    header = lines[1]
    rows = lines[2:end]
    @test length(rows) == 2
    col(name) = findfirst(==(name), header)
    mu_u_idx = col("mu_u_MeV")
    mu_d_idx = col("mu_d_MeV")
    rho_ratio_idx = col("rho_u_over_rho_d")
    mu_pi_idx = col("mu_pi_MeV")
    @test mu_u_idx !== nothing
    @test mu_d_idx !== nothing
    @test rho_ratio_idx !== nothing
    @test mu_pi_idx !== nothing

    ok_rows = [row for row in rows if row[col("status")] == "ok"]
    @test !isempty(ok_rows)
    first_ok = ok_rows[1]
    @test parse(Float64, first_ok[mu_u_idx]) != parse(Float64, first_ok[mu_d_idx])
    @test parse(Float64, first_ok[rho_ratio_idx]) ≈ 0.876 atol=5e-3
    @test parse(Float64, first_ok[mu_pi_idx]) != 0.0
end
