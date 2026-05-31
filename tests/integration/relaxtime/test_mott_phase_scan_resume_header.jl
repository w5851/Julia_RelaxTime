using Test

const REPO_ROOT_MOTT_RESUME = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MOTT_PHASE_SCAN_SCRIPT = joinpath(REPO_ROOT_MOTT_RESUME, "scripts", "relaxtime", "run_mott_phase_scan.jl")

@testset "Mott phase scan resume header contract" begin
    tmp = mktempdir()
    outdir = joinpath(tmp, "out")
    mkpath(outdir)
    csv_path = joinpath(outdir, "mott_phase_scan.csv")
    write(csv_path, join([
        "# format: scan_csv_v1",
        "run_id,T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,status",
        "old,100,0,0,1,1,0,0,ok",
    ], "\n") * "\n")

    cmd = ignorestatus(`$(Base.julia_cmd()) --project=$REPO_ROOT_MOTT_RESUME $MOTT_PHASE_SCAN_SCRIPT --outdir $outdir --tmin 100 --tmax 100 --xi-list 0 --resume`)
    log_path = joinpath(tmp, "run.log")
    open(log_path, "w") do io
        run(pipeline(cmd; stdout=io, stderr=io))
    end
    result = read(log_path, String)
    @test occursin("existing output CSV header is incompatible with current schema", result)
    @test occursin("--overwrite", result)
end
