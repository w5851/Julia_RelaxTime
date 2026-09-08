using Test
const _DF_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_DF_ROOT,"scripts","analysis","relaxtime","audit_causal_gbu_dense_freezeout.jl"))
const testcausalgbudensefreezeout_DF=CausalGBUDenseFreezeout
@testset "Dense freezeout preserves failures and ratio semantics" begin
    @test testcausalgbudensefreezeout_DF.energy_grid("200,3,7.7")==[3.,7.7,200.]
    for bad in ("0,3","NaN,3","3,3","-1,5","broken")
        @test_throws ArgumentError testcausalgbudensefreezeout_DF.energy_grid(bad)
    end
    input=[(sqrt_s_NN_GeV=3.,channel=c,density_inv_fm3=n,passed=true) for
        (c,n) in (("pi_plus",2.),("K_plus",1.),("pi_minus",4.),("K_minus",0.5))]
    r=testcausalgbudensefreezeout_DF.ratio_table(input,input,[3.,5.])
    @test r[1].direct_Kplus_over_pi_plus==0.5
    @test r[1].reference_Kminus_over_pi_minus==0.125
    @test r[1].direct_plus_passed && r[1].direct_minus_passed
    @test isnan(r[2].direct_Kplus_over_pi_plus) && !r[2].direct_plus_passed
    failed=copy(input);failed[2]=merge(input[2],(passed=false,))
    @test !testcausalgbudensefreezeout_DF.ratio_table(failed,input,[3.])[1].direct_plus_passed
    negative=copy(input);negative[2]=merge(input[2],(density_inv_fm3=-1.,))
    n=testcausalgbudensefreezeout_DF.ratio_table(negative,input,[3.])[1]
    @test n.direct_Kplus_over_pi_plus==-0.5 && !n.direct_plus_passed
    @test_throws ErrorException testcausalgbudensefreezeout_DF.ratio_table(vcat(input,input),input,[3.])
    @test all(!x.production_authorized for x in r)
end
