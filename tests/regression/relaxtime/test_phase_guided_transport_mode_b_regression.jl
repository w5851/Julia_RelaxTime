using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_phase_guided_transport_mode_b_v1.csv")
const SCRIPT_PATH = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")

if !isdefined(Main, :run_phase_guided_scan)
    include(SCRIPT_PATH)
end

function _parse_csv_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    in_quotes = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            next_i = nextind(line, i)
            if in_quotes && next_i <= lastindex(line) && line[next_i] == '"'
                print(buf, '"')
                i = next_i
            else
                in_quotes = !in_quotes
            end
        elseif c == ',' && !in_quotes
            push!(fields, String(take!(buf)))
        else
            print(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end

function _read_csv_rows(path::String)
    header = String[]
    rows = Dict{String,String}[]
    open(path, "r") do io
        for raw in eachline(io)
            line = strip(raw)
            isempty(line) && continue
            startswith(line, "#") && continue
            if isempty(header)
                header = _parse_csv_line(line)
                continue
            end
            values = _parse_csv_line(line)
            length(values) == length(header) || error("invalid csv row in $(path): $(line)")
            row = Dict{String,String}()
            for (k, v) in zip(header, values)
                row[k] = v
            end
            push!(rows, row)
        end
    end
    return rows
end

function _baseline_map(path::String)
    rows = _read_csv_rows(path)
    out = Dict{Tuple{Float64,Float64,Float64}, NamedTuple}()
    for row in rows
        key = (
            parse(Float64, row["T_MeV"]),
            parse(Float64, row["muB_MeV"]),
            parse(Float64, row["xi"]),
        )
        out[key] = (
            tau_u=parse(Float64, row["tau_u"]),
            tau_s=parse(Float64, row["tau_s"]),
            eta=parse(Float64, row["eta"]),
            sigma=parse(Float64, row["sigma"]),
            eta_over_s=parse(Float64, row["eta_over_s"]),
            sigma_over_T=parse(Float64, row["sigma_over_T"]),
        )
    end
    return out
end

@testset "phase-guided transport mode b regression" begin
    tmp = mktempdir()
    opts = Main.PhaseGuidedTransportScanCLI.PhaseGuidedScanOptions(
        :mode_b_fixed_T_sparse_muB,
        tmp,
        "regression_mode_b_v1",
        [-0.5, -0.2],
        [0.0, 260.0],
        Float64[1.0],
        [120.0, 138.0],
        false,
        true,
        false,
    )
    ctx = Main.ProvenanceMetadata.new_run_context("tests/regression/relaxtime/test_phase_guided_transport_mode_b_regression.jl", String[])
    Main.run_phase_guided_scan(opts, ctx)

    result_csv = joinpath(tmp, "phase_guided_transport_scan.csv")
    @test isfile(result_csv)

    rows = _read_csv_rows(result_csv)
    observed = Dict{Tuple{Float64,Float64,Float64}, Dict{String,String}}()
    for row in rows
        key = (
            parse(Float64, row["T_MeV"]),
            parse(Float64, row["muB_MeV"]),
            parse(Float64, row["xi"]),
        )
        observed[key] = row
    end

    baseline = _baseline_map(BASELINE_PATH)
    rtol = 1e-6
    atol = 1e-8
    for (key, expected) in baseline
        @test haskey(observed, key)
        row = observed[key]
        @test isapprox(parse(Float64, row["tau_u"]), expected.tau_u; rtol=rtol, atol=atol)
        @test isapprox(parse(Float64, row["tau_s"]), expected.tau_s; rtol=rtol, atol=atol)
        @test isapprox(parse(Float64, row["eta"]), expected.eta; rtol=rtol, atol=atol)
        @test isapprox(parse(Float64, row["sigma"]), expected.sigma; rtol=rtol, atol=atol)
        @test isapprox(parse(Float64, row["eta_over_s"]), expected.eta_over_s; rtol=rtol, atol=atol)
        @test isapprox(parse(Float64, row["sigma_over_T"]), expected.sigma_over_T; rtol=rtol, atol=atol)
    end
end
