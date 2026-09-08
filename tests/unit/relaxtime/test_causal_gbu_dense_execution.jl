using Test
const _DE_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_DE_ROOT,"scripts","analysis","relaxtime","causal_gbu_dense_execution.jl"))
const testcausalgbudenseexecution_DE=CausalGBUDenseExecution
@testset "Reference reuse matches slow finite-q-discarding oracle" begin
    bg=(m=(u=1.0,d=1.1,s=1.4),mu=(u=0.15,d=0.1,s=0.05),T=0.8,Phi=0.3,PhiBar=0.4,
        coupling=Dict(c=>0.20 for c in testcausalgbudenseexecution_DE.R.CHANNELS),vacuum=3.0)
    s=testcausalgbudenseexecution_DE.R.Settings(mesh=16,np=16,nx=8,ne=16,nw=80,nq=8,thermal=6.0)
    for channel in testcausalgbudenseexecution_DE.R.CHANNELS
        b=testcausalgbudenseexecution_DE.R.bubble_at(bg,channel,0.0,s)
        ref=(bubble=b,gap=testcausalgbudenseexecution_DE.R.gap_audit(b,s))
        for q in (0.0,0.8,3.0)
            fast=testcausalgbudenseexecution_DE.reference_shell(bg,channel,q,s,ref)
            slow=testcausalgbudenseexecution_DE.R.shell(bg,channel,q,s;reference=ref)
            @test isequal(fast,slow)
        end
    end
    a=testcausalgbudenseexecution_DE.reference_density(bg,:K_plus,s)
    b=testcausalgbudenseexecution_DE.R.density(bg,:K_plus,s;reference=true)
    @test isequal(a,b)
end
