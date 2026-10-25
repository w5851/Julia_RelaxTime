using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","TerminalTransportQuadrature.jl"))
using .TerminalTransportQuadrature, .ThresholdUpgradeCommon

@testset "Gram moments reproduce the actual RTA kernels" begin
    qp=(m=(u=0.73,d=0.81,s=1.41),μ=(u=0.29,d=0.25,s=0.57))
    tau=(u=0.7,d=0.8,s=1.4,ubar=0.9,dbar=0.6,sbar=1.1)
    coeffs=(masses=[qp.m.u,qp.m.d,qp.m.s],v_n_sq=0.23,dμB_dT_sigma=-0.39,
        dM_dT=[-1.1,-0.9,-0.5],dM_dμB=[-0.17,-0.13,-0.04])
    p,pw=Main.RelaxTime.GaussLegendre.gauleg(0.0,15.0,32)
    ct,cw=Main.RelaxTime.GaussLegendre.gauleg(-1.0,1.0,16)
    config=Main.RelaxTime.TransportCoefficients.TransportIntegrationConfig(p_grid=p,p_w=pw,cos_grid=ct,cos_w=cw)
    TC=Main.RelaxTime.TransportCoefficients
    for xi in (0.0,-0.49,0.37)
        tp=(T=0.82,Φ=0.38,Φbar=0.47,ξ=xi)
        grams=terminal_grams(qp,tp,p,pw,ct,cw)
        result=transport_from_grams(grams,qp,tp,tau;bulk_coeffs=coeffs)
        @test isapprox(result.eta,TC.shear_viscosity(qp,tp;tau=tau,config=config);rtol=2e-13)
        @test isapprox(result.sigma,TC.electric_conductivity(qp,tp;tau=tau,config=config);rtol=2e-13)
        @test isapprox(result.zeta,TC.bulk_viscosity_isentropic(qp,tp;tau=tau,config=config,bulk_coeffs_isentropic=coeffs);rtol=2e-13)
        fp,fw=Main.RelaxTime.GaussLegendre.gauleg(0.0,20.0,64)
        fine=terminal_grams(qp,tp,fp,fw,ct,cw)
        for s in SPECIES
            G=getproperty(grams,s);R=getproperty(fine,s)
            bound=gram_relative_change(G,R)
            for b in ([0.2,1.4,-0.3],[-0.13,-0.6,1.1],[0.0,1.0,0.0])
                @test abs(dot(b,(G-R)*b))/dot(b,R*b)<=bound+1e-12
            end
        end
    end
end
