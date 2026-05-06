using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_external_path_meson_density_scan.jl")
const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "meson_density", "external_path", "test_outputs")
const OUTFILE = joinpath(OUTDIR, "external_path_meson_density_scan_smoke.csv")
const INFILE = joinpath(OUTDIR, "external_path_points_smoke.csv")

@testset "external-path meson density scan smoke" begin
    mkpath(OUTDIR)
    open(INFILE, "w") do io
        println(io, "source_fig,case_id,line_style,point_index,muB_MeV,T_MeV")
        println(io, "smoke_fig,smoke_case,dashed,0,30.0,170.0")
    end
    isfile(OUTFILE) && rm(OUTFILE)

    cmd = `julia --project=. $SCRIPT --input-csv $INFILE --source-fig smoke_fig --case-id smoke_case --line-style dashed --output $OUTFILE --overwrite --p-num 8 --t-num 4 --max-iter 20`
    run(cmd)

    @test isfile(OUTFILE)
    text = read(OUTFILE, String)
    @test occursin("path_source,path_case_id,path_line_style,path_point_index,path_order_key", text)
    @test occursin("smoke_fig,smoke_case,dashed,0", text)
    @test occursin("kpi_ratio", text)
end
