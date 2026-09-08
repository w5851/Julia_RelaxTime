using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","audit_causal_gbu_infinite_readiness.jl"))
const testcausalgbuinfinitereadiness_IRG=CausalGBUInfiniteReadiness

@testset "Weak-limit reducer rejects incomplete and nonconvergent evidence" begin
    rows=[(channel="toy",q_inv_fm=1.,eta_inv_fm=e,lower_inv_fm=l,
        density=1-e,pv_density=1.,relative_difference=e) for e in (0.003,0.001,0.0003) for l in (0.001,0.0003)]
    @test testcausalgbuinfinitereadiness_IRG.eta_gate(rows)
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate(rows[1:0])
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate(rows[1:end-1])
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate(vcat(rows,rows[1:1]))
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate([merge(r,(relative_difference=0.02,)) for r in rows])
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate([merge(r,(relative_difference=0.0001/r.eta_inv_fm,)) for r in rows])
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate([merge(r,(density=r.lower_inv_fm==0.001 ? 2. : 1.,)) for r in rows])
    @test !testcausalgbuinfinitereadiness_IRG.eta_gate([merge(r,(relative_difference=NaN,)) for r in rows])
end
