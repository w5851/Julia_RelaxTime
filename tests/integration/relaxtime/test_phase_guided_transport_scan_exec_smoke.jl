using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")

function _non_comment_lines(path::String)
    lines = String[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            push!(lines, s)
        end
    end
    return lines
end

@testset "phase guided transport one-point execution smoke" begin
    outdir = mktempdir()
    cmd = `julia --project=. $SCRIPT_PATH --mode fixed-T-sparse-muB --outdir $outdir --case-name exec_smoke --xi-list 0.0 --muB-list 0 --T-list 120 --overwrite`
    run(cmd)

    result_csv = joinpath(outdir, "phase_guided_transport_scan.csv")
    failed_csv = joinpath(outdir, "failed_points.csv")

    @test isfile(result_csv)
    @test isfile(failed_csv)

    result_lines = _non_comment_lines(result_csv)
    @test length(result_lines) == 2
    @test occursin("mode_b_fixed_T_sparse_muB", result_lines[2])
    @test occursin(",crossover,", result_lines[2])

    failed_lines = _non_comment_lines(failed_csv)
    @test length(failed_lines) == 1
    @test startswith(failed_lines[1], "T_MeV,muB_MeV,xi,")
end
