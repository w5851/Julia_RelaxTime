using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","compute_charged_gbu_bqs_point.jl"))
const BQSPointAdapter=ChargedGBUBQSPoint
@testset "Fixed-point adapter preserves existing BQS/GBU contracts" begin
    opts=BQSPointAdapter.parse_args(String[])
    @test opts[:T_MeV]==50.0
    @test opts[:muB_MeV]==620.0
    @test BQSPointAdapter.parse_args(["--workers","1"])[:workers]==1
    for args in (["--T-MeV","0"],["--T-MeV","NaN"],["--muB-MeV","Inf"],
                 ["--workers","5"],["--bad","x"],["--output"])
        @test_throws ArgumentError BQSPointAdapter.parse_args(args)
    end
    mktempdir() do dir
        @test_throws ErrorException BQSPointAdapter.run_point(output=dir)
        for ch in BQSPointAdapter.W.R.CHANNELS
            BQSPointAdapter.W.writejson(joinpath(dir,"$(ch).json"),
                (passed=true,status="accepted",density=startswith(String(ch),"K") ? .05 : 1.))
        end
        r=BQSPointAdapter.summarize(dir,(T_MeV=50.,muB_MeV=620.))
        @test r.Kplus_over_pi_plus==.05 && r.plus_passed
        @test !hasproperty(r,:sqrt_s_NN_GeV)
        BQSPointAdapter.W.writejson(joinpath(dir,"pi_plus.json"),
            (passed=false,status="local_gate_failed",density=NaN))
        r=BQSPointAdapter.summarize(dir,(T_MeV=50.,muB_MeV=620.))
        @test isnan(r.Kplus_over_pi_plus) && !r.plus_passed
        @test r.Kminus_over_pi_minus==.05 && r.minus_passed
    end
end
