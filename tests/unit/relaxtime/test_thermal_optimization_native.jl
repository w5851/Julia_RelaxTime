using Test
if !isdefined(Main, :ThresholdCacheCandidates)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "ThresholdCacheCandidates.jl"))
end

using .ThresholdUpgradeCommon
using .ThresholdCacheCandidates

@testset "Medium-frame Møller speed keeps the event energies" begin
    # Equal masses and perpendicular incoming momenta with |p_i|=|p_j|=m=1:
    # s=6, v_cm=2/sqrt(3), while the medium-frame Møller definition is
    # sqrt(3)/2.  This is the algebraic counterexample used in the convention
    # audit; it must not be reduced to a function of s alone.
    mi=mj=1.0
    ei=ej=sqrt(2.0)
    s=6.0
    vcm=AS._relative_velocity_for_mode(:cm,s,mi,mj,ei,ej)
    vm=AS._relative_velocity_for_mode(:medium_moller,s,mi,mj,ei,ej)
    @test isapprox(vcm,2/sqrt(3);rtol=1e-14,atol=0.0)
    @test isapprox(vm,sqrt(3)/2;rtol=1e-14,atol=0.0)
    @test isapprox(vcm/vm,4/3;rtol=1e-14,atol=0.0)
    @test AS._relative_velocity_for_mode("cm",s,mi,mj,ei,ej)==vcm
    @test_throws ArgumentError AS._relative_velocity_for_mode(:unknown,s,mi,mj,ei,ej)
end

@testset "Factored equal-mass relative velocity and threshold" begin
    for m in (0.0,0.01,1.0,3.0), x in (1e-12,1e-8,1e-3,1.0,100.0)
        s=4m^2+x
        expected=2sqrt((s-4m^2)/s)
        @test AS._relative_velocity_cm(s,m,m)==expected
        @test 0<=AS._relative_velocity_cm(s,m,m)<=2
    end
    for (mi,mj) in ((0.1,2.0),(1.0,1.0+1e-10)), x in (1e-6,0.1,10.0)
        s=(mi+mj)^2+x
        invariant=2s*sqrt((s-(mi+mj)^2)*(s-(mi-mj)^2))/(s^2-(mi^2-mj^2)^2)
        @test isapprox(AS._relative_velocity_cm(s,mi,mj),invariant;rtol=1e-8)
    end
end
