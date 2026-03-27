using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ORCH_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_relaxtime_orchestrator.jl")

function _count_data_rows(path::String)
    count = 0
    open(path, "r") do io
        header_seen = false
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !header_seen
                header_seen = true
                continue
            end
            count += 1
        end
    end
    return count
end

@testset "cross-section orchestrated smoke" begin
    @test isfile(ORCH_PATH)

    outdir = mktempdir()
    cfg = joinpath(outdir, "cfg.toml")
    open(cfg, "w") do io
        write(io, """
schema_version = "v1"

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0]
xi_list = [0.0, 0.1]
processes = ["ud_to_ud", "us_to_us"]

[scan.cross_section.energy]
mode = "list"
sqrt_s_list_MeV = [500.0, 600.0]
""")
    end

    cmd = `julia --project=. $ORCH_PATH cross-section --config $cfg --outdir $outdir`
    run(cmd)

    out_csv = joinpath(outdir, "cross_section_orchestrated.csv")
    @test isfile(out_csv)
    @test _count_data_rows(out_csv) == 8
end
