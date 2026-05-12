using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCAN_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_phase_shift_meson_thermo_scan.jl")
const PLOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "build_meson_thermo_plot_review_case.py")
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_meson_thermo_plot_review_case_v1.csv")

function _load_metric_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    metrics = Dict{String,NamedTuple{(:expected, :rtol, :atol),Tuple{Float64,Float64,Float64}}}()
    for line in readlines(path)[2:end]
        row = strip(line)
        isempty(row) && continue
        startswith(row, "#") && continue
        cols = split(row, ',')
        length(cols) == 4 || error("invalid baseline row: $line")
        metrics[strip(cols[1])] = (
            expected=parse(Float64, strip(cols[2])),
            rtol=parse(Float64, strip(cols[3])),
            atol=parse(Float64, strip(cols[4])),
        )
    end
    return metrics
end

function _read_csv_rows(path::String)
    isfile(path) || error("CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("CSV is empty: $path")

    header = nothing
    start_idx = 0
    for (idx, line) in enumerate(lines)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        header = split(s, ',')
        start_idx = idx + 1
        break
    end
    header === nothing && error("CSV has no header: $path")

    rows = NamedTuple[]
    for line in lines[start_idx:end]
        row = strip(line)
        isempty(row) && continue
        startswith(row, "#") && continue
        cols = split(row, ','; keepempty=true)
        length(cols) == length(header) || error("invalid CSV row: $line")
        push!(rows, (; (Symbol(header[i]) => cols[i] for i in eachindex(header))...))
    end
    return rows
end

function _run_scan(outdir::String)
    run(`julia --project=. $SCAN_SCRIPT --outdir $outdir --overwrite --tmin 170 --tmax 230 --tstep 20 --q-nodes 6 --omega-nodes 6 --qmax 4 --omega-max 3 --p-num 8 --t-num 4 --max-iter 20`)
end

function _run_plot_review(source_dir::String, outdir::String)
    python_cmd = Sys.which("python") !== nothing ? "python" : (Sys.which("python3") !== nothing ? "python3" : nothing)
    python_cmd === nothing && return false
    success(`$(python_cmd) -c "import matplotlib"`) || return false
    run(`$(python_cmd) $PLOT_SCRIPT --source-dir $source_dir --out-dir $outdir`)
    return true
end

function _metrics(rows)
    qp_share = [parse(Float64, row.qp_share) for row in rows]
    ld_share = [parse(Float64, row.ld_share) for row in rows]
    p_total = [parse(Float64, row.P_total) for row in rows]
    trace = [parse(Float64, row.trace_anomaly) for row in rows]
    return Dict(
        "row_count" => Float64(length(rows)),
        "first_T" => parse(Float64, first(rows).T_MeV),
        "last_T" => parse(Float64, last(rows).T_MeV),
        "max_P_total" => maximum(p_total),
        "min_trace_anomaly" => minimum(trace),
        "max_ld_share" => maximum(ld_share),
        "min_qp_share" => minimum(qp_share),
    )
end

@testset "Meson thermo plot-review canonical case regression" begin
    @test isfile(SCAN_SCRIPT)
    @test isfile(PLOT_SCRIPT)

    baseline = _load_metric_baseline(BASELINE_PATH)
    source_dir = mktempdir()
    outdir = mktempdir()
    _run_scan(source_dir)

    if !_run_plot_review(source_dir, outdir)
        @test_skip "python/matplotlib unavailable for plot-review regression"
        return
    end

    workflow_rows = _read_csv_rows(joinpath(outdir, "workflow_scan.csv"))
    summary_rows = _read_csv_rows(joinpath(outdir, "plot_review_summary.csv"))
    source_rows = _read_csv_rows(joinpath(source_dir, "scan.csv"))

    @test length(workflow_rows) == length(source_rows) == 4
    @test length(summary_rows) == 4

    for (source, copied) in zip(source_rows, workflow_rows)
        @test source.T_MeV == copied.T_MeV
        @test source.P_total == copied.P_total
        @test source.trace_anomaly == copied.trace_anomaly
        @test source.thermo_derivation_mode == copied.thermo_derivation_mode
    end

    metrics = _metrics(summary_rows)
    for metric in sort!(collect(keys(baseline)))
        @test haskey(metrics, metric)
        expected = baseline[metric]
        @test isapprox(metrics[metric], expected.expected; rtol=expected.rtol, atol=expected.atol)
    end

    readme = read(joinpath(outdir, "README.md"), String)
    @test occursin("manual trend inspection", readme)
    @test occursin("not an external validation gate", readme)
    @test occursin("max P_total", readme)
    @test occursin("min trace anomaly", readme)

    for fig in ("pressure_overlay.png", "trace_anomaly_overlay.png", "qp_ld_split.png")
        path = joinpath(outdir, fig)
        @test isfile(path)
        @test filesize(path) > 10_000
    end
end
