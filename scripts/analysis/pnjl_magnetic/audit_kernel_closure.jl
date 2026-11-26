#!/usr/bin/env julia

#= Fixed-state closure of complete production Omega derivatives against the
independent reference, after explicit p and zeta refinement. Input selection
reuses the worst previous p1024 density cases, one vacuum anchor, and one
multi-root low-T case in each route. All computations run on quote-remote.
=#
using JSON, LinearAlgebra, StaticArrays, ForwardDiff
include(joinpath(@__DIR__,"MagneticReference.jl"))
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
if !isdefined(Main,:Models)
    include(joinpath(ROOT,"src","models","Models.jl"))
end
const MT=Models.MagneticThermodynamics
const MI=MT.MagneticIntegrals
const P=MagneticReference.Parameters()
const H=MagneticReference.HBARC

function production_at(x,T,mu,eB,route,p_num,zeta_num,nmax)
    conf=MT.MagneticConfig(eB_fm2=eB,n_max=nmax,p_num=p_num,
        pz_max=40.0,zeta_num=zeta_num,route=route)
    st=SVector{5,Float64}(x)
    muv=SVector{3,Float64}(mu)
    ev=MT.calculate_magnetic_omega_components(st,muv,T,conf)
    loop=zeros(3)
    for f in 1:3
        q=MagneticReference.CHARGES[f]
        loop[f]=ForwardDiff.derivative(ev.masses[f]) do m
            vacuum=if route===:mfir
                -2P.colors*Models.PNJLCore.vacuum_integral_with_cutoff(m,convert(typeof(m),P.cutoff))+
                    MI.omega_magnetic_mfir(m,q,eB;zeta_num=zeta_num)
            else
                MI.omega0_flavor_landau(m,q,eB;n_max=nmax,p_num=p_num,pz_max=40.0)
            end
            vacuum+MI.omegat_flavor_landau(m,mu[f],T,x[4],x[5],q,eB;
                n_max=nmax,p_num=p_num,pz_max=40.0)
        end
    end
    model=Models.PNJLMagneticModel(Models.PNJLModel(),conf)
    gradient=Models.magnetic_gap_residual(model,st,T,muv;p_num=p_num,n_max=nmax)
    rho=MT.calculate_magnetic_rho(st,muv,T,conf)
    J=ForwardDiff.jacobian(v->Models.calculate_mass_vec(model,v),x[1:3])
    gap=loop-x[1:3]
    return (;omega=ev.omega,rho=collect(rho),gap,gradient=collect(gradient),
        chain_rule_error=norm(gradient[1:3]-J'*gap),
        mass_jacobian_sigma_min=minimum(svdvals(J)))
end

function main(args=ARGS)
    length(args)==3 || throw(ArgumentError("representative states, refined folds, output"))
    input,foldpath,output=args
    states=JSON.parsefile(input)["states"]
    folds=JSON.parsefile(foldpath)["folds"]
    for f in folds
        r=only(v["result"] for v in f["numerical"] if v["id"]=="o128_n96")
        push!(states,Dict("id"=>f["id"],"route"=>"landau_legacy",
            "eB_GeV2"=>.05,"T_MeV"=>5.0,"muB_MeV"=>r["muB_MeV"],"x_state"=>r["x_state"]))
    end
    out=Dict{String,Any}("julia_version"=>string(VERSION),"input"=>input,
        "fold_source"=>foldpath,"states"=>Any[])
    mkpath(dirname(output))
    for s in states
        x=Float64.(s["x_state"]); T=s["T_MeV"]/H
        mu=fill(s["muB_MeV"]/(3H),3);eB=s["eB_GeV2"]*(1000/H)^2
        route=Symbol(s["route"])
        ctl=MagneticReference.Controls(order=128,n_max=96,route=route)
        reference=MagneticReference.evaluate(x,T,mu,eB,P,ctl)
        checks=Any[]
        counts=s["id"] in ("first_minimum","first_maximum","second_minimum",
            "second_maximum","main_minimum") ? (512,1024,2048) : (128,512,1024,2048,4096)
        for p in counts
            pr=production_at(x,T,mu,eB,route,p,512,96)
            push!(checks,(;p_num=p,zeta_num=512,n_max=96,
                omega_delta=pr.omega-reference.omega,rho_delta=norm(pr.rho-reference.rho),
                gap_delta=norm(pr.gap-reference.condensate_gap),
                gradient_delta=norm(pr.gradient-reference.stationarity),
                chain_rule_error=pr.chain_rule_error))
            println(s["id"]," p=",p," density delta=",checks[end].rho_delta,
                " gap delta=",checks[end].gap_delta);flush(stdout)
        end
        high=MagneticReference.evaluate(x,T,mu,eB,P,
            MagneticReference.Controls(order=192,n_max=128,route=route))
        tail=MagneticReference.evaluate(x,T,mu,eB,P,
            MagneticReference.Controls(order=128,n_max=96,route=route,
                thermal_width=3.,thermal_tail=50.,pz_max=60.))
        push!(out["states"],(;input=s,reference_gap=norm(reference.residual),checks,
            reference_high_gap_delta=norm(high.condensate_gap-reference.condensate_gap),
            reference_high_rho_delta=norm(high.rho-reference.rho),
            reference_tail_gap_delta=norm(tail.condensate_gap-reference.condensate_gap),
            reference_tail_rho_delta=norm(tail.rho-reference.rho)))
        open(io->JSON.print(io,out,2),output,"w")
    end
end

function audit_zeta_floor(input,output)
    d=JSON.parsefile(input)
    s=only(s for s in d["states"] if s["id"]=="mfir_worst_previous_p1024")
    x=Float64.(s["x_state"])
    eB=s["eB_GeV2"]*(1000/H)^2
    masses=MagneticReference.masses(x,eB,P)
    checks=Any[]
    for f in 1:3
        mass=masses[f]; q=MagneticReference.CHARGES[f]
        vref,gref=MagneticReference.vacuum(mass,q*eB,P,MagneticReference.Controls())
        variants=Any[]
        for z in (128,256,512,1024,2048)
            vf=m->-2P.colors*Models.PNJLCore.vacuum_integral_with_cutoff(m,convert(typeof(m),P.cutoff))+
                MI.omega_magnetic_mfir(m,q,eB;zeta_num=z)
            push!(variants,(;zeta_num=z,omega_delta=vf(mass)-vref,
                condensate_delta=ForwardDiff.derivative(vf,mass)-gref))
        end
        push!(checks,(;flavor=f,mass_MeV=mass*H,zeta_argument=mass^2/(2q*eB),variants))
    end
    out=JSON.parsefile(output)
    mu=fill(s["muB_MeV"]/(3H),3)
    reference=MagneticReference.evaluate(x,s["T_MeV"]/H,mu,eB,P,
        MagneticReference.Controls(order=128,n_max=96))
    pr=production_at(x,s["T_MeV"]/H,mu,eB,:mfir,4096,1024,96)
    combined=(;p_num=4096,zeta_num=1024,n_max=96,
        omega_delta=pr.omega-reference.omega,rho_delta=norm(pr.rho-reference.rho),
        gap_delta=norm(pr.gap-reference.condensate_gap),
        gradient_delta=norm(pr.gradient-reference.stationarity),
        chain_rule_error=pr.chain_rule_error)
    out["small_mass_zeta_checks"]=(;input=s,checks,combined)
    open(io->JSON.print(io,out,2),output,"w")
    println(JSON.json(checks))
end

function audit_configuration(output)
    x=SVector{5,Float64}(-1.8,-1.8,-2.2,.001,.001)
    T=5/H; mu=SVector{3,Float64}(fill(500/(3H),3));eB=.05*(1000/H)^2
    probes=Any[]
    for route in (:mfir,:landau_legacy)
        model=Models.PNJLMagneticModel(eB_fm2=eB,profile="unittest_lambda",
            route=route,p_num=1024,n_max=96,pz_max=40.,zeta_num=1024)
        pa=model.magnetic.params
        ref=MagneticReference.Parameters(cutoff=pa.Λ_inv_fm,G0=pa.G_fm2,
            K=pa.K_fm5,ml=pa.m_ud0_inv_fm,ms=pa.m_s0_inv_fm,
            T0=pa.T0_inv_fm)
        legacy_default_cutoff=MagneticReference.Parameters(cutoff=P.cutoff,
            G0=pa.G_fm2,K=pa.K_fm5,ml=pa.m_ud0_inv_fm,ms=pa.m_s0_inv_fm,
            T0=pa.T0_inv_fm)
        ctl=MagneticReference.Controls(order=128,n_max=96,route=route)
        prod=MT.calculate_magnetic_omega_components(x,mu,T,model.magnetic)
        expected=MagneticReference.evaluate(x,T,mu,eB,ref,ctl)
        frozen=MagneticReference.evaluate(x,T,mu,eB,legacy_default_cutoff,ctl)
        push!(probes,(;route=String(route),requested_cutoff_MeV=pa.Λ_inv_fm*H,
            mass_delta=norm(prod.masses-expected.masses),production_vacuum=prod.vac,
            configured_reference_vacuum=expected.vac,
            default_cutoff_reference_vacuum=frozen.vac,
            configured_vacuum_delta=prod.vac-expected.vac,
            frozen_cutoff_vacuum_delta=prod.vac-frozen.vac))
    end
    budgets=Any[]
    for field in (.05,.155,.4)
        conf=MT.MagneticConfig(eB_fm2=field*(1000/H)^2,route=:landau_legacy)
        push!(budgets,(;eB_GeV2=field,n_max=MT.resolve_magnetic_nmax(T,mu,conf)))
    end
    out=JSON.parsefile(output)
    out["configuration_checks"]=(;state=collect(x),T_MeV=5.,muB_MeV=500.,probes,
        legacy_automatic_thermal_budgets=budgets)
    open(io->JSON.print(io,out,2),output,"w")
    println(JSON.json(out["configuration_checks"]))
end

if abspath(PROGRAM_FILE)==@__FILE__
    if length(ARGS)==3 && ARGS[3]=="--zeta-only"
        audit_zeta_floor(ARGS[1],ARGS[2])
    elseif length(ARGS)==2 && ARGS[2]=="--configuration-only"
        audit_configuration(ARGS[1])
    else
        main()
    end
end
