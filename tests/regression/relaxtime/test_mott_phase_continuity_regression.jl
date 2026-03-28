using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_mott_phase_scan.jl")

function _load_rows(path::String)
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
    return rows
end

function _count_large_jumps(rows::Vector{Dict{String,String}}, field::String; threshold::Float64)
    by_xi = Dict{Float64,Vector{Tuple{Float64,Float64}}}()
    for r in rows
        r["status"] == "ok" || continue
        xi = parse(Float64, r["xi"])
        T = parse(Float64, r["T_MeV"])
        v = tryparse(Float64, r[field])
        v === nothing && continue
        isfinite(v) || continue
        get!(by_xi, xi, Tuple{Float64,Float64}[])
        push!(by_xi[xi], (T, v))
    end

    count = 0
    for (_, seq0) in by_xi
        seq = sort(seq0; by=x -> x[1])
        for i in 2:length(seq)
            if abs(seq[i][2] - seq[i - 1][2]) > threshold
                count += 1
            end
        end
    end
    return count
end

@testset "mott phase continuity regression" begin
    @test isfile(SCRIPT_PATH)

    outdir = mktempdir()
    run(`julia --project=. $SCRIPT_PATH --outdir $outdir --overwrite --tmin 120 --tmax 260 --tstep 10 --p-num 8 --t-num 4 --max-iter 20`)

    out_csv = joinpath(outdir, "mott_phase_scan.csv")
    @test isfile(out_csv)
    rows = _load_rows(out_csv)
    @test !isempty(rows)

    # Continuation should suppress most branch-switch artifacts.
    # Keep a moderate threshold to avoid over-constraining true physical transitions.
    n_jump_mk = _count_large_jumps(rows, "M_K"; threshold=0.8)
    n_jump_mpi = _count_large_jumps(rows, "M_pi"; threshold=0.8)

    @test n_jump_mk <= 2
    @test n_jump_mpi <= 2
end
