using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_phase_shift_meson_thermo_scan.jl")

function _non_comment_data_lines(text::String)
    return [
        line for line in split(text, '\n')
        if !isempty(strip(line)) && !startswith(strip(line), "#") && !startswith(strip(line), "T_MeV,")
    ]
end

@testset "phase-shift meson thermo canonical scan smoke" begin
    @test isfile(SCRIPT)

    outdir = mktempdir()
    cmd = `julia --project=. $SCRIPT --outdir $outdir --overwrite --tmin 210 --tmax 210 --tstep 5 --q-nodes 8 --omega-nodes 8 --qmax 4 --omega-max 3 --p-num 8 --t-num 4 --max-iter 20`
    run(cmd)

    out_csv = joinpath(outdir, "scan.csv")
    readme = joinpath(outdir, "README.md")
    eff_cfg = joinpath(outdir, "effective_config.json")
    manifest = joinpath(outdir, "run_manifest.json")

    @test isfile(out_csv)
    @test isfile(readme)
    @test isfile(eff_cfg)
    @test isfile(manifest)

    text = read(out_csv, String)
    @test occursin("# workflow_entry: Models.solve_gap_and_phase_shift_meson_thermo_point", text)
    @test occursin("# canonical_case: muB=0", text)
    @test occursin("P_total", text)
    @test occursin("trace_anomaly", text)
    @test occursin("P_meson_qp", text)
    @test occursin("P_meson_ld", text)
    @test occursin("thermo_derivation_mode", text)
    @test occursin("primary_channel,secondary_channel", text)

    data_lines = _non_comment_data_lines(text)
    @test length(data_lines) == 1

    cols = split(data_lines[1], ',')
    @test length(cols) == 51
    @test cols[9] == "pi"
    @test cols[10] == "sigma_pi"
    @test cols[41] == "omega_total_ad"

    readme_text = read(readme, String)
    @test occursin("canonical mu_B = 0 phase-shift meson thermo case", readme_text)
    @test occursin("scan.csv", readme_text)
    @test occursin("run_manifest.json", readme_text)
    @test occursin("sigma_pi", readme_text)

    cfg_text = read(eff_cfg, String)
    @test occursin("\"workflow_entry\":\"Models.solve_gap_and_phase_shift_meson_thermo_point\"", cfg_text)
    @test occursin("\"channels\":[\"pi\",\"sigma_pi\"]", cfg_text)

    manifest_obj = JSON3.read(read(manifest, String))
    @test String(manifest_obj["script"]) == "scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl"
    @test length(manifest_obj["artifacts"]) == 2
end
