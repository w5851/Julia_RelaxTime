using Test
using SHA

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

function _read_text(path::String)
    isfile(path) || error("regression fixture not found: $path")
    return read(path, String)
end

function _sha1_hex(path::String)
    return bytes2hex(sha1(_read_text(path)))
end

function _parse_probe_summary(path::String)
    lines = split(chomp(_read_text(path)), '\n')
    length(lines) >= 2 || error("probe summary is empty: $path")
    header = split(lines[1], ',')
    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) == length(header) || error("invalid probe summary row: $line")
        push!(rows, (; (Symbol(header[i]) => cols[i] for i in eachindex(header))...))
    end
    return rows
end

function _find_row(rows, center_xi::String, sigma_grid_n::String, adaptive_refinement::String)
    for row in rows
        if row.center_xi == center_xi && row.sigma_grid_n == sigma_grid_n && row.adaptive_refinement == adaptive_refinement && row.scenario == "baseline"
            return row
        end
    end
    error("probe summary row not found: center_xi=$center_xi sigma_grid_n=$sigma_grid_n adaptive_refinement=$adaptive_refinement")
end

@testset "tau xi probe regression fixtures" begin
    t190_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "_xi_probe_T190_summary.csv")
    t200_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "_xi_probe_T200_summary.csv")

    @test _sha1_hex(t190_path) == "127841e4bd0e5ea0d488271527073963da96fbe1"
    @test _sha1_hex(t200_path) == "4114daeded42f46774063ee9bc5089619f20a0e0"

    t190_rows = _parse_probe_summary(t190_path)
    t200_rows = _parse_probe_summary(t200_path)

    t190_row = _find_row(t190_rows, "-0.44", "128", "true")
    t200_row = _find_row(t200_rows, "-0.2", "256", "true")

    @test parse(Float64, t190_row.kink_metric_tauinv) < 0.2
    @test parse(Float64, t190_row.kink_metric_tau) < 0.35
    @test parse(Float64, t200_row.kink_metric_tauinv) < 0.02
    @test parse(Float64, t200_row.kink_metric_tau) < 0.02
end
