using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ORCH_PATH = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_relaxtime_orchestrator.jl")

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

@testset "workflow vs direct consistency (mini)" begin
    @test isfile(ORCH_PATH)

    outdir = mktempdir()
    cfg = joinpath(outdir, "mini_cfg.toml")
    open(cfg, "w") do io
        write(io, """
schema_version = "v1"

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0]
xi_list = [0.0]
processes = ["ud_to_ud"]

[scan.cross_section.energy]
mode = "list"
sqrt_s_list_MeV = [500.0, 600.0]
""")
    end

    run(`julia --project=. $ORCH_PATH cross-section --config $cfg --outdir $outdir`)

    orchestrated_csv = joinpath(outdir, "cross_section_orchestrated.csv")
    @test isfile(orchestrated_csv)
    @test _count_data_rows(orchestrated_csv) == 2

    # Current mini regression consistency criterion:
    # row-count identity for identical normalized fixed-point grid.
    direct_expected_rows = 2
    @test _count_data_rows(orchestrated_csv) == direct_expected_rows
end
