using Test
if !isdefined(Main,:full_rs_identity)
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","run_transport_full_rs_layer_validation.jl"))
end

@testset "full-RS layer identity separates rounding from physics" begin
    state=(quark_params=(A=(u=-3.0,d=-3.0,s=-2.0),),K_coeffs=(K0_plus=0.5,K0_minus=0.75))
    @test full_rs_identity(state,state).passed
    rounded=merge(state,(K_coeffs=(K0_plus=nextfloat(0.5),K0_minus=0.75),))
    identity=full_rs_identity(state,rounded)
    @test identity.passed
    @test identity.A_exact
    @test 0<identity.coupling_max_relative<=identity.coupling_roundoff_limit<1e-12
    changed=merge(state,(K_coeffs=(K0_plus=0.5*(1+1e-12),K0_minus=0.75),))
    @test !full_rs_identity(state,changed).passed
    changed_a=merge(state,(quark_params=(A=(u=nextfloat(-3.0),d=-3.0,s=-2.0),),))
    @test !full_rs_identity(state,changed_a).passed
    for invalid in (NaN,Inf,-Inf)
        @test !full_rs_identity(state,merge(state,(K_coeffs=(K0_plus=invalid,K0_minus=0.75),))).passed
    end
    @test_throws ErrorException full_rs_identity(state,merge(state,(K_coeffs=(K0_plus=0.5,),)))
end
