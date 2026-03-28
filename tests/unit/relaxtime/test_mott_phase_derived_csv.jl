using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_mott_phase_derived_csv.jl")

function _read_csv_rows(path::String)
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

@testset "mott derived csv builder" begin
    @test isfile(SCRIPT_PATH)

    tmp = mktempdir()
    input_csv = joinpath(tmp, "input.csv")
    output_csv = joinpath(tmp, "derived.csv")

    open(input_csv, "w") do io
        write(io, "# format: scan_csv_v1\n")
        write(io, "T_MeV,muB_MeV,xi,m_u,m_d,m_s,status,error_code,error_message,timestamp_utc\n")
        write(io, "150.0,0.0,0.0,1.1,1.2,1.3,ok,,,2026-03-28T00:00:00Z\n")
        write(io, "150.0,0.0,0.15,NaN,1.0,1.5,error,E_SOLVE,solver failed,2026-03-28T00:00:01Z\n")
    end

    run(`julia --project=. $SCRIPT_PATH --in $input_csv --out $output_csv`)

    @test isfile(output_csv)
    header, rows = _read_csv_rows(output_csv)
    @test "M_u_plus_M_d" in header
    @test "M_u_plus_M_s" in header
    @test length(rows) == 2

    @test parse(Float64, rows[1]["M_u_plus_M_d"]) ≈ 2.3 atol = 1e-12
    @test parse(Float64, rows[1]["M_u_plus_M_s"]) ≈ 2.4 atol = 1e-12

    @test rows[2]["status"] == "error"
    @test rows[2]["error_code"] == "E_SOLVE"
    @test rows[2]["M_u_plus_M_d"] == "NaN"
    @test rows[2]["M_u_plus_M_s"] == "NaN"
end
