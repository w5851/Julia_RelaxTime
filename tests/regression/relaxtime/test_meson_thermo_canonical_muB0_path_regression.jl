using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(
    PROJECT_ROOT,
    "tests", "baselines", "relaxtime", "baseline_meson_thermo_canonical_muB0_path_v1.csv",
)
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_phase_shift_meson_thermo_scan.jl")

function _read_csv_rows(path::String)
    isfile(path) || error("CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("CSV is empty: $path")

    header_line = nothing
    start_idx = 0
    for (idx, line) in enumerate(lines)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        header_line = s
        start_idx = idx + 1
        break
    end
    header_line === nothing && error("CSV has no header: $path")
    header = [strip(col) for col in split(header_line, ',')]

    rows = NamedTuple[]
    for line in lines[start_idx:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ','; keepempty=true)
        length(cols) == length(header) || error("invalid CSV row: $line")
        push!(rows, (; (Symbol(header[i]) => strip(cols[i]) for i in eachindex(header))...))
    end
    return rows
end

function _run_scan(outdir::String)
    cmd = `julia --project=. $SCRIPT --outdir $outdir --overwrite --tmin 170 --tmax 230 --tstep 20 --q-nodes 6 --omega-nodes 6 --qmax 4 --omega-max 3 --p-num 8 --t-num 4 --max-iter 20`
    run(cmd)
    return (
        scan_csv=joinpath(outdir, "scan.csv"),
        readme=joinpath(outdir, "README.md"),
        manifest=joinpath(outdir, "run_manifest.json"),
    )
end

function _as_float(row, key::Symbol)
    return parse(Float64, getproperty(row, key))
end

function _as_bool(row, key::Symbol)
    return lowercase(getproperty(row, key)) == "true"
end

@testset "Meson thermo canonical muB0 path regression" begin
    @test isfile(SCRIPT)
    baseline = _read_csv_rows(BASELINE_PATH)
    outdir = mktempdir()
    paths = _run_scan(outdir)
    rows = _read_csv_rows(paths.scan_csv)

    @test length(rows) == length(baseline) == 4

    for (expected, actual) in zip(baseline, rows)
        @test _as_float(actual, :T_MeV) == _as_float(expected, :T_MeV)
        @test _as_float(actual, :muB_MeV) == _as_float(expected, :muB_MeV)
        @test _as_float(actual, :xi) == _as_float(expected, :xi)
        @test actual.workflow == expected.workflow
        @test actual.channel_set == expected.channel_set
        @test actual.primary_channel == expected.primary_channel
        @test actual.secondary_channel == expected.secondary_channel
        @test actual.phase_shift_variant == expected.phase_shift_variant
        @test actual.thermo_derivation_mode == expected.thermo_derivation_mode
        @test _as_bool(actual, :equilibrium_converged) == _as_bool(expected, :equilibrium_converged)

        @test isapprox(_as_float(actual, :P_meson), _as_float(expected, :P_meson); rtol=5e-3, atol=1e-8)
        @test isapprox(_as_float(actual, :P_meson_qp), _as_float(expected, :P_meson_qp); rtol=5e-3, atol=1e-8)
        @test isapprox(_as_float(actual, :P_meson_ld), _as_float(expected, :P_meson_ld); rtol=5e-3, atol=1e-8)
        @test isapprox(_as_float(actual, :P_total), _as_float(expected, :P_total); rtol=5e-3, atol=1e-8)
        @test isapprox(_as_float(actual, :trace_anomaly), _as_float(expected, :trace_anomaly); rtol=1e-2, atol=1e-8)

        @test isapprox(
            _as_float(actual, :P_meson),
            _as_float(actual, :P_meson_qp) + _as_float(actual, :P_meson_ld);
            rtol=1e-10, atol=1e-10
        )
        @test actual.thermo_derivation_mode != "workflow_fd_legacy"
    end

    readme_text = read(paths.readme, String)
    @test occursin("canonical mu_B = 0 phase-shift meson thermo case", readme_text)
    @test occursin("points_total: `4`", readme_text)
    @test occursin("success_count: `4`", readme_text)
    @test occursin("error_count: `0`", readme_text)
    @test occursin("sigma_pi", readme_text)

    manifest_obj = JSON3.read(read(paths.manifest, String))
    @test String(manifest_obj["script"]) == "scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl"
    @test Int(manifest_obj["summary"]["points_total"]) == 4
    @test Int(manifest_obj["summary"]["success_count"]) == 4
    @test Int(manifest_obj["summary"]["error_count"]) == 0
end
