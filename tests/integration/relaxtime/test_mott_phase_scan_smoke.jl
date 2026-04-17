using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_mott_phase_scan.jl")

function _read_csv_header_and_rows(path::String)
    header = String[]
    rows = Vector{Dict{String,String}}()
    open(path, "r") do io
        header_seen = false
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !header_seen
                header = [strip(x) for x in split(s, ',')]
                header_seen = true
                continue
            end
            vals = split(s, ',')
            d = Dict{String,String}()
            for (i, k) in enumerate(header)
                d[k] = i <= length(vals) ? strip(vals[i]) : ""
            end
            push!(rows, d)
        end
    end
    return header, rows
end

@testset "mott phase scan smoke" begin
    @test isfile(SCRIPT_PATH)

    outdir = mktempdir()
    cmd = `julia --project=. $SCRIPT_PATH --outdir $outdir --overwrite --tmin 150 --tmax 150 --tstep 10 --p-num 8 --t-num 4 --max-iter 30`
    run(cmd)

    out_csv = joinpath(outdir, "mott_phase_scan.csv")
    eff_cfg = joinpath(outdir, "effective_config.json")
    manifest = joinpath(outdir, "run_manifest.json")

    @test isfile(out_csv)
    @test isfile(eff_cfg)
    @test isfile(manifest)

    header, rows = _read_csv_header_and_rows(out_csv)
    @test !isempty(rows)

    required_cols = Set([
        "T_MeV", "muB_MeV", "xi",
        "M_pi", "M_K", "Gamma_pi", "Gamma_K",
        "residual_pi", "residual_K",
        "root_quality_pi", "root_quality_K",
        "selected_method_pi", "selected_method_K",
        "m_u", "m_d", "m_s",
        "status", "error_code", "error_message", "timestamp_utc",
    ])
    @test required_cols ⊆ Set(header)

    xi_vals = sort(unique(parse.(Float64, [r["xi"] for r in rows])))
    @test xi_vals == [-0.3, -0.15, 0.0, 0.15, 0.3]

    muB_vals = unique(parse.(Float64, [r["muB_MeV"] for r in rows]))
    @test muB_vals == [0.0]

    cfg_text = read(eff_cfg, String)
    @test occursin("\"T_min_MeV\":150", cfg_text)
    @test occursin("\"T_max_MeV\":150", cfg_text)
    @test occursin("\"T_step_MeV\":10", cfg_text)
    @test occursin("\"xi_list\":[-0.3,-0.15,0.0,0.15,0.3]", cfg_text)

    manifest_obj = JSON3.read(read(manifest, String))
    config_path = String(manifest_obj["config_path"])
    output_csv_path = String(manifest_obj["output_csv"])
    @test !occursin('\\', config_path)
    @test !occursin('\\', output_csv_path)
    @test !occursin(r"^[A-Za-z]:", config_path)
    @test !occursin(r"^[A-Za-z]:", output_csv_path)
    @test occursin("config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml", config_path)
    @test endswith(output_csv_path, "mott_phase_scan.csv")
end
