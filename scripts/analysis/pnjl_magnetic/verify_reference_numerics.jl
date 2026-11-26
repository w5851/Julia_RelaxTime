#!/usr/bin/env julia
# Independent derivative, quadrature and Landau-budget checks for saved PALC
# branches. These are diagnostics, not production-baseline tolerances.
using JSON, LinearAlgebra
using StaticArrays, ForwardDiff
include(joinpath(@__DIR__,"MagneticReference.jl"))
using .MagneticReference
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
if !isdefined(Main,:Models)
    include(joinpath(ROOT,"src","models","Models.jl"))
end
const H=MagneticReference.HBARC
const P=MagneticReference.Parameters()
const MT=Models.MagneticThermodynamics
const MI=MT.MagneticIntegrals

function _state_checks(x,T,mu,eB,route)
    controls=MagneticReference.Controls(order=48,route=route)
    base=MagneticReference.evaluate(x,T,mu,eB,P,controls)
    variants=Any[]
    for order in (24,32,64,96)
        ev=MagneticReference.evaluate(x,T,mu,eB,P,MagneticReference.Controls(order=order,route=route))
        push!(variants,(;kind="order",value=order,omega_delta=ev.omega-base.omega,
            gap_delta=norm(ev.condensate_gap-base.condensate_gap),rho_delta=norm(ev.rho-base.rho),
            residual_norm=norm(ev.residual)))
    end
    for nmax in (16,48,64)
        ev=MagneticReference.evaluate(x,T,mu,eB,P,MagneticReference.Controls(order=64,route=route,n_max=nmax))
        push!(variants,(;kind="n_max",value=nmax,omega_delta=ev.omega-base.omega,
            gap_delta=norm(ev.condensate_gap-base.condensate_gap),rho_delta=norm(ev.rho-base.rho),
            residual_norm=norm(ev.residual)))
    end
    for pz in (10.0,20.0,60.0)
        ev=MagneticReference.evaluate(x,T,mu,eB,P,MagneticReference.Controls(order=64,route=route,pz_max=pz))
        push!(variants,(;kind="pz_max",value=pz,omega_delta=ev.omega-base.omega,
            gap_delta=norm(ev.condensate_gap-base.condensate_gap),rho_delta=norm(ev.rho-base.rho),
            residual_norm=norm(ev.residual)))
    end
    derivative_errors=Any[]
    for h in (1e-4,1e-5)
        grad=zeros(5)
        rho=zeros(3)
        for j in 1:5
            xp,xm=copy(x),copy(x)
            xp[j]+=h;xm[j]-=h
            grad[j]=(MagneticReference.evaluate(xp,T,mu,eB,P,controls).omega-
                MagneticReference.evaluate(xm,T,mu,eB,P,controls).omega)/(2h)
        end
        for j in 1:3
            mp,mm=copy(mu),copy(mu)
            mp[j]+=h;mm[j]-=h
            rho[j]=-(MagneticReference.evaluate(x,T,mp,eB,P,controls).omega-
                MagneticReference.evaluate(x,T,mm,eB,P,controls).omega)/(2h)
        end
        push!(derivative_errors,(;h,omega_gradient_delta=norm(grad-base.stationarity),
            rho_derivative_delta=norm(rho-base.rho)))
    end
    production=Any[]
    for count in (128,512,1024)
        conf=MT.MagneticConfig(eB_fm2=eB,n_max=32,p_num=count,pz_max=40.0,route=route)
        st=SVector{5,Float64}(x)
        muv=SVector{3,Float64}(mu)
        comp=MT.calculate_magnetic_omega_components(st,muv,T,conf)
        dens=MT.calculate_magnetic_rho(st,muv,T,conf)
        model=Models.PNJLMagneticModel(Models.PNJLModel(),conf)
        res=Models.magnetic_gap_residual(model,st,T,muv;p_num=count,n_max=32,pz_max=40.0,t_num=4)
        push!(production,(;p_num=count,omega_delta=comp.omega-base.omega,
            rho_delta=norm(dens-base.rho),stationarity_norm=norm(res),
            rhoB_over_rho0=sum(dens)/(3*.16)))
    end
    return (;base_residual_norm=norm(base.residual),base_gap_norm=norm(base.condensate_gap),
        mass_jacobian_sigma_min=base.mass_jacobian_sigma_min,variants,derivative_errors,production)
end

function vacuum_checks()
    out=Any[]
    for ebG in (0.05,0.155,0.4), M in (10/H,350/H,550/H), q in (1/3,2/3)
        eB=ebG*(1000/H)^2
        _,expected=MagneticReference.vacuum(M,q*eB,P,MagneticReference.Controls())
        prod=zeros(4)
        for (i,z) in enumerate((64,256,512,1024))
            prod[i]=ForwardDiff.derivative(M) do m
                -2P.colors*Models.PNJLCore.vacuum_integral_with_cutoff(m,convert(typeof(m),P.cutoff))+
                    MI.omega_magnetic_mfir(m,q,eB;zeta_num=z)
            end
        end
        push!(out,(;eB_GeV2=ebG,M_MeV=M*H,q,expected,zeta_counts=[64,256,512,1024],errors=prod.-expected))
    end
    return out
end

function main(args=ARGS)
    length(args)==3 || throw(ArgumentError("MFIR file, smooth file and output required"))
    out=Any[]
    for fn in args[1:2]
        for c in JSON.parsefile(fn)["cases"]
            haskey(c,"points") || continue
            route=Symbol(c["route"])
            T=c["T_MeV"]/H
            eB=c["eB_GeV2"]*(1000/H)^2
            # Vacuum, two competing low-T roots, and a high-density strange
            # threshold root are chosen before examining numerical errors.
            selected=Any[c["points"][1]]
            for t in c["targets"]
                t["muB_MeV"] in (1050.0,1500.0) || continue
                append!(selected,[merge(r,Dict("muB_MeV"=>t["muB_MeV"])) for r in t["palc_roots"] if r["physical"]])
            end
            for r in selected
                println("verify ",route," eb=",c["eB_GeV2"]," muB=",r["muB_MeV"]);flush(stdout)
                mu=fill(r["muB_MeV"]/(3H),3)
                value=_state_checks(r["x_state"],T,mu,eB,route)
                push!(out,(;route=String(route),eB_GeV2=c["eB_GeV2"],muB_MeV=r["muB_MeV"],
                    x_state=r["x_state"],checks=value))
            end
        end
    end
    open(args[3],"w") do io
        JSON.print(io,(;states=out,vacuum_checks=vacuum_checks()),2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
