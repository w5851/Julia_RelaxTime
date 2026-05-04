using Test

const _LITERATURE_TARGETS_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "literature", "meson_density_targets",
))

const _EXPECTED_LITERATURE_TARGETS = [
    "blaschke2019col_kminus_piminus_mu_pi_100_fig4_right",
    "blaschke2019col_kminus_piminus_mu_pi_134p5_fig4_right",
    "blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_with_anomalous",
    "blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_no_anomalous",
]

function _read_xy_csv(path::String)
    isfile(path) || error("literature target CSV missing: $path")
    lines = readlines(path)
    isempty(lines) && error("literature target CSV is empty: $path")
    strip(lines[1]) == "x_value,y_value" || error("unexpected CSV header in $path")

    xs = Float64[]
    ys = Float64[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) == 2 || error("invalid literature target row in $path: $line")
        push!(xs, parse(Float64, strip(cols[1])))
        push!(ys, parse(Float64, strip(cols[2])))
    end
    return xs, ys
end

@testset "Meson density literature targets smoke" begin
    @test isdir(_LITERATURE_TARGETS_DIR)

    manifest = joinpath(_LITERATURE_TARGETS_DIR, "manifest_wpd_import_report.csv")
    @test isfile(manifest)
    manifest_text = read(manifest, String)
    @test occursin("clamped_negative_to_zero", manifest_text)

    for slug in _EXPECTED_LITERATURE_TARGETS
        csv_path = joinpath(_LITERATURE_TARGETS_DIR, slug * ".csv")
        meta_path = joinpath(_LITERATURE_TARGETS_DIR, slug * ".meta.md")

        @test isfile(csv_path)
        @test isfile(meta_path)

        xs, ys = _read_xy_csv(csv_path)
        @test !isempty(xs)
        @test length(xs) == length(ys)
        @test issorted(xs)
        @test all(diff(xs) .> 0.0)
        @test all(isfinite, xs)
        @test all(isfinite, ys)
        @test all(y -> y >= 0.0, ys)

        meta_text = read(meta_path, String)
        @test occursin("source_paper: `Blaschke:2019col`", meta_text)
        @test occursin("path_context:", meta_text)
        @test occursin("extraction_method:", meta_text)
    end
end
