using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_meson_density_plot_review_case_v1.csv")
const SOURCE_COMPARISON_PATH = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "results", "relaxtime", "meson_density", "freezeout_validation",
    "blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default", "comparison_vs_target.csv",
)
const PLOT_REVIEW_DIR = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "results", "relaxtime", "meson_density", "plot_review",
    "freezeout_kminus_piminus_mu_pi_100",
)
const PLOT_REVIEW_FIGURE_DIR = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "figures", "relaxtime", "meson_density", "plot_review",
    "freezeout_kminus_piminus_mu_pi_100",
)
const PLOT_COMPARISON_PATH = joinpath(PLOT_REVIEW_DIR, "plot_review_comparison.csv")
const README_PATH = joinpath(PLOT_REVIEW_DIR, "README.md")
const OVERLAY_PNG_PATH = joinpath(PLOT_REVIEW_FIGURE_DIR, "overlay_kminus_piminus_mu_pi_100.png")
const RESIDUAL_PNG_PATH = joinpath(PLOT_REVIEW_FIGURE_DIR, "residual_kminus_piminus_mu_pi_100.png")
const PLOT_MANIFEST_PATH = joinpath(PLOT_REVIEW_FIGURE_DIR, "plot_manifest.json")

function _load_metric_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    metrics = Dict{String, NamedTuple{(:expected, :rtol, :atol), Tuple{Float64, Float64, Float64}}}()
    for line in readlines(path)[2:end]
        row = strip(line)
        isempty(row) && continue
        startswith(row, "#") && continue
        cols = split(row, ',')
        length(cols) == 4 || error("invalid baseline row: $line")
        metric = strip(cols[1])
        expected = parse(Float64, strip(cols[2]))
        rtol = parse(Float64, strip(cols[3]))
        atol = parse(Float64, strip(cols[4]))
        metrics[metric] = (expected=expected, rtol=rtol, atol=atol)
    end
    return metrics
end

function _read_csv_rows(path::String)
    isfile(path) || error("CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    rows = NamedTuple[]
    for line in lines[2:end]
        row = strip(line)
        isempty(row) && continue
        cols = split(row, ','; keepempty=true)
        length(cols) == length(header) || error("invalid CSV row: $line")
        push!(rows, (; (Symbol(header[i]) => cols[i] for i in eachindex(header))...))
    end
    return rows
end

_maybe_float(s::AbstractString) = isempty(strip(s)) ? nothing : parse(Float64, strip(s))

function _comparison_metrics(rows)
    abs_vals = Float64[parse(Float64, row.abs_diff) for row in rows]
    rel_vals = Float64[]
    for row in rows
        val = _maybe_float(row.rel_diff)
        isnothing(val) || push!(rel_vals, val)
    end
    return Dict(
        "row_count" => Float64(length(rows)),
        "first_x" => parse(Float64, first(rows).x_value),
        "last_x" => parse(Float64, last(rows).x_value),
        "first_model_y" => parse(Float64, first(rows).model_y),
        "last_model_y" => parse(Float64, last(rows).model_y),
        "max_abs_diff" => maximum(abs_vals),
        "max_rel_diff" => maximum(rel_vals),
    )
end

function _compare_source_and_plot_review(source_rows, plot_rows)
    length(source_rows) == length(plot_rows) || return false
    for (source, plot) in zip(source_rows, plot_rows)
        parse(Float64, source.x_value) == parse(Float64, plot.x_value) || return false
        parse(Float64, source.target_y) == parse(Float64, plot.target_y) || return false
        parse(Float64, source.model_y) == parse(Float64, plot.model_y) || return false
        parse(Float64, source.abs_diff) == parse(Float64, plot.abs_diff) || return false

        source_rel_raw = strip(source.rel_diff)
        plot_rel_raw = strip(plot.rel_diff)
        if source_rel_raw == "NaN" && isempty(plot_rel_raw)
            # plot-review CSV deliberately normalizes undefined relative diff to empty string.
        else
            source_rel = _maybe_float(source.rel_diff)
            plot_rel = _maybe_float(plot.rel_diff)
            source_rel == plot_rel || return false
        end

        parse(Float64, source.T_MeV) == parse(Float64, plot.T_MeV) || return false
        parse(Float64, source.muB_MeV) == parse(Float64, plot.muB_MeV) || return false
    end
    return true
end

@testset "Meson density plot-review canonical case regression" begin
    baseline = _load_metric_baseline(BASELINE_PATH)
    source_rows = _read_csv_rows(SOURCE_COMPARISON_PATH)
    plot_rows = _read_csv_rows(PLOT_COMPARISON_PATH)
    metrics = _comparison_metrics(plot_rows)

    @test _compare_source_and_plot_review(source_rows, plot_rows)

    for metric in sort!(collect(keys(baseline)))
        expected = baseline[metric]
        @test haskey(metrics, metric)
        @test isapprox(metrics[metric], expected.expected; rtol=expected.rtol, atol=expected.atol)
    end

    readme = read(README_PATH, String)
    @test occursin("manual trend inspection", readme)
    @test occursin("not an external validation gate", readme)
    @test occursin("points: `48`", readme)
    @test occursin("max abs diff: `0.147787`", readme)
    @test occursin("max rel diff (finite only): `0.873012`", readme)

    @test isfile(OVERLAY_PNG_PATH)
    @test isfile(RESIDUAL_PNG_PATH)
    @test isfile(PLOT_MANIFEST_PATH)
    @test filesize(OVERLAY_PNG_PATH) > 10_000
    @test filesize(RESIDUAL_PNG_PATH) > 10_000
end
