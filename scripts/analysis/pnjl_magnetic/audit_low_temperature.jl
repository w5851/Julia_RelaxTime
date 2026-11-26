#!/usr/bin/env julia
#= Compare the magnetic production kernel with independent low-T formulas.

Run under the root environment, preferably on quote-remote. Results remain
diagnostic; no production profile or acceptance target is changed.
=#
using JSON
using LinearAlgebra: norm
using StaticArrays
using ForwardDiff
include(joinpath(@__DIR__, "MagneticReference.jl"))
using .MagneticReference
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main,:Models)
    include(joinpath(ROOT,"src","models","Models.jl"))
end
const MT = Models.MagneticThermodynamics
const MI = MT.MagneticIntegrals
const REF_P = MagneticReference.Parameters()
const H = MagneticReference.HBARC

function records(path)
    data=JSON.parsefile(path)
    values=Any[]
    for c in data["cases"]
        haskey(c,"targets") || continue
        for target in c["targets"]
            target["muB_MeV"] in (900.0,1050.0,1150.0,1500.0) || continue
            for (idx,r) in enumerate(target["catalog_roots"])
                push!(values,(;eB=c["eB_GeV2"],T=c["T_MeV"],muB=target["muB_MeV"],
                    index=idx,x=r["x_state"]))
            end
        end
    end
    return values
end

function audit_state(r)
    eB=r.eB*(1000/H)^2
    T=r.T/H
    mu=fill(r.muB/(3H),3)
    x=SVector{5,Float64}(r.x)
    ref16=MagneticReference.evaluate(r.x,T,mu,eB,REF_P,
        MagneticReference.Controls(order=16))
    ref32=MagneticReference.evaluate(r.x,T,mu,eB,REF_P,
        MagneticReference.Controls(order=32))
    comparisons=Any[]
    for count in (32,64,128,256,512)
        conf=MT.MagneticConfig(eB_fm2=eB,n_max=32,p_num=count,pz_max=40.0)
        comp=MT.calculate_magnetic_omega_components(x,SVector{3,Float64}(mu),T,conf)
        # Independent loop gap evaluated on the same unsplit quadrature:
        # first derivative w.r.t mass plus the analytic MFIR loggamma term.
        loops=zeros(3)
        for f in 1:3
            _,lv=MagneticReference.vacuum(comp.masses[f],eB*MagneticReference.CHARGES[f],
                REF_P,MagneticReference.Controls())
            thermal_derivative=ForwardDiff.derivative(comp.masses[f]) do M
                MI.omegat_flavor_landau(M,mu[f],T,x[4],x[5],MagneticReference.CHARGES[f],eB;
                    n_max=32,p_num=count,pz_max=40.0)
            end
            loops[f]=lv+thermal_derivative
        end
        model=Models.PNJLMagneticModel(Models.PNJLModel(),conf)
        residual=Models.magnetic_gap_residual(model,x,T,SVector{3,Float64}(mu);
            p_num=count,n_max=32,pz_max=40.0,t_num=4)
        rho=MT.calculate_magnetic_rho(x,SVector{3,Float64}(mu),T,conf)
        push!(comparisons,(;p_num=count,omega=comp.omega,vac=comp.vac,therm=comp.therm,
            rho=collect(rho),stationarity_norm=norm(residual),
            condensate_gap_norm=norm(loops-r.x[1:3]),
            omega_reference_delta=comp.omega-ref32.omega,
            rho_reference_delta=norm(rho-ref32.rho)))
    end
    # Same unsplit nodes: true self consistency vs Omega stationarity can be
    # distinguished even if quadrature is poorly resolved.
    return merge(r,(;reference=ref32,
        reference_order_rho_delta=norm(ref16.rho-ref32.rho),
        reference_order_gap_delta=norm(ref16.condensate_gap-ref32.condensate_gap),
        reference_order_omega_delta=abs(ref16.omega-ref32.omega),comparisons))
end

function fixed_kernel_checks()
    values=Any[]
    for x in (0.01,0.1,0.5,1.0,2.0,8.0,25.0)
        val=MagneticReference.zeta_derivative(x)
        actual=MI.zeta_prime_minus_one(x;zeta_num=256)
        push!(values,(;kind="zeta",x,reference=val,production=actual,error=actual-val))
    end
    for eBG in (0.05,0.155,0.4), f in (1,2), MMeV in (10.0,350.0)
        eB=eBG*(1000/H)^2
        q=MagneticReference.CHARGES[f]
        M=MMeV/H
        conf=MagneticReference.Controls()
        _,loop=MagneticReference.vacuum(M,q*eB,REF_P,conf)
        production=ForwardDiff.derivative(M) do mass
            -2REF_P.colors*Models.PNJLCore.vacuum_integral_with_cutoff(mass,convert(typeof(mass),REF_P.cutoff))+
                MI.omega_magnetic_mfir(mass,q,eB;zeta_num=256)
        end
        push!(values,(;kind="vacuum_condensate",eB_GeV2=eBG,flavor=f,M_MeV=MMeV,
            reference=loop,production,error=production-loop))
    end
    for Phi in (0.0,0.25), delta in (-20.0,-1.0,0.0,1.0,20.0)
        E,mu,T=1.0,1.0+0.05delta,0.05
        ref=MagneticReference.occupation((mu-E)/T,Phi,Phi)
        r2=MagneticReference.occupation((-mu-E)/T,Phi,Phi)
        prod_log,prod_net=MI._polyakov_log_and_net_density(E,mu,T,Phi,Phi)
        push!(values,(;kind="occupation",Phi,delta,
            log_error=prod_log-ref[1]-r2[1],net_error=prod_net-3(ref[2]-r2[2])))
    end
    return values
end

function main(args=ARGS)
    length(args)==2 || throw(ArgumentError("input PALC JSON and output JSON required"))
    checks=fixed_kernel_checks()
    states=Any[]
    for r in records(args[1])
        println("audit eB=",r.eB," muB=",r.muB," root=",r.index);flush(stdout)
        push!(states,audit_state(r))
    end
    mkpath(dirname(args[2]))
    open(args[2],"w") do io
        JSON.print(io,(;julia_version=string(VERSION),checks,states),2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
