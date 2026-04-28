using Test

const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _SCRIPT_PATH = joinpath(_REPO_ROOT, "scripts", "relaxtime", "merge_xi_smoothness_review.jl")

function _write_csv(path::String, header::String, rows::Vector{String})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, row)
        end
    end
end

function _parse_csv_line(raw::AbstractString)
    vals = String[]
    buf = IOBuffer()
    in_quotes = false
    i = firstindex(raw)
    while i <= lastindex(raw)
        c = raw[i]
        if c == '"'
            if in_quotes
                ni = nextind(raw, i)
                if ni <= lastindex(raw) && raw[ni] == '"'
                    write(buf, '"')
                    i = ni
                else
                    in_quotes = false
                end
            else
                in_quotes = true
            end
        elseif c == ',' && !in_quotes
            push!(vals, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(raw, i)
    end
    push!(vals, String(take!(buf)))
    return vals
end

function _read_csv_rows(path::String)
    header = String[]
    rows = Vector{Dict{String, String}}()
    open(path, "r") do io
        header_seen = false
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !header_seen
                header = [strip(x) for x in _parse_csv_line(s)]
                header_seen = true
                continue
            end
            vals = [strip(x) for x in _parse_csv_line(line)]
            row = Dict{String, String}()
            for (i, k) in enumerate(header)
                row[k] = i <= length(vals) ? vals[i] : ""
            end
            push!(rows, row)
        end
    end
    return header, rows
end

function _find_row(rows::Vector{Dict{String, String}}, sample_id::String, field::String)
    matched = [r for r in rows if r["sample_id"] == sample_id && r["field"] == field]
    length(matched) == 1 || error("expected exactly one row for $(sample_id),$(field)")
    return matched[1]
end

@testset "merge xi smoothness review cli" begin
    @test isfile(_SCRIPT_PATH)

    tmp = mktempdir()
    flags_csv = joinpath(tmp, "flags.csv")
    manual_csv = joinpath(tmp, "manual_review.csv")
    out_csv = joinpath(tmp, "smoothness_final.csv")

    _write_csv(
        flags_csv,
        "sample_id,field,label,reason,S2,S1jump,N_spike",
        [
            "S001,tau_u,suspect,S2_above_smooth,0.091,0.21,1",
            "S001,tau_s,smooth,all_metrics_within_smooth_thresholds,0.01,0.03,0",
            "S002,tau_u,not_smooth,S2_above_suspect,0.42,0.63,3",
        ],
    )
    _write_csv(
        manual_csv,
        "sample_id,field,manual_label,reason",
        [
            "S001,tau_u,confirm_not_smooth,\"manual_visual_spike, near xi=\"\"-0.34\"\"\"",
            "S001,tau_s,confirm_smooth,manual_curve_consistent",
        ],
    )

    run(`julia --project=. $_SCRIPT_PATH --flags $flags_csv --manual-review $manual_csv --out $out_csv`)
    @test isfile(out_csv)

    header, rows = _read_csv_rows(out_csv)
    @test header == ["sample_id", "field", "auto_label", "manual_label", "final_label", "reason"]
    @test length(rows) == 3

    r_override = _find_row(rows, "S001", "tau_u")
    @test r_override["auto_label"] == "suspect"
    @test r_override["manual_label"] == "confirm_not_smooth"
    @test r_override["final_label"] == "confirm_not_smooth"
    @test r_override["reason"] == "manual_visual_spike, near xi=\"-0.34\""

    r_auto = _find_row(rows, "S002", "tau_u")
    @test r_auto["auto_label"] == "not_smooth"
    @test r_auto["manual_label"] == ""
    @test r_auto["final_label"] == "not_smooth"

    invalid_manual_csv = joinpath(tmp, "manual_review_invalid.csv")
    _write_csv(
        invalid_manual_csv,
        "sample_id,field,manual_label,reason",
        ["S001,tau_u,recheck_later,bad_label"],
    )

    err = IOBuffer()
    proc = run(pipeline(
        ignorestatus(`julia --project=. $_SCRIPT_PATH --flags $flags_csv --manual-review $invalid_manual_csv --out $(joinpath(tmp, "should_not_exist.csv"))`);
        stdout=devnull,
        stderr=err,
    ))
    @test proc.exitcode != 0
    msg = String(take!(err))
    @test occursin("invalid manual_label", msg)
end
