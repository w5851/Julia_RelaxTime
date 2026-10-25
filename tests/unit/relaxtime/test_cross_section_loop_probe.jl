using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","CrossSectionLoopProbe.jl"))
using .CrossSectionLoopProbe, .ThresholdUpgradeCommon

@testset "Explicit loop probe covers all channels and restores the production cache" begin
    state=retained_state((m_u=1.4,m_d=1.4,m_s=2.2,muq_fm=0.13,T_fm=0.8,Phi=0.4,Phibar=0.6,xi=0.2))
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    pc=Main.RelaxTime.PolarizationCache
    backend=Main.RelaxTime.PolarizationAniso.polarization_aniso
    probe=LoopProbe(backend)
    for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
        particles=CS.parse_particles_from_process(process)
        m=map(p->AS.get_mass(p,qp),particles)
        s=max((m[1]+m[2])^2,(m[3]+m[4])^2)+2.3
        pc.reset_cache!()
        baseline=CS.total_cross_section(process,s,qp,tp,kc;n_points=4)
        saved=copy(pc.POLARIZATION_CACHE)
        actual=probe_cross_section(probe,process,s,state;n_points=4)
        @test isapprox(actual,baseline;rtol=2e-12,atol=1e-14)
        @test pc.POLARIZATION_CACHE==saved
        failed=LoopProbe((args...)->(NaN,0.0))
        @test_throws DomainError probe_cross_section(failed,process,s,state;n_points=4)
        @test pc.POLARIZATION_CACHE==saved
    end
    @test_throws ArgumentError LoopProbe(backend;channels=:unknown)
end
