#!/usr/bin/env julia

"""Refine the five preselected Fig. 5-3 folds and audit explicit conventions.

This isolated diagnostic uses the independent MagneticReference equations.
It changes no production defaults or numerical baselines. A fold is solved as
F(x,mu)=0 together with det(F_x)=0; the singular-value scale is fixed from
the initial estimate. Separate finite-difference steps and neighboring points
check that the result is a turning point, not an undersampled PALC vertex.
"""

using JSON, NLsolve, LinearAlgebra, Dates
include(joinpath(@__DIR__, "MagneticReference.jl"))
using .MagneticReference

const H = MagneticReference.HBARC
const P0 = MagneticReference.Parameters()
const FEATURE_IDS = ("first_minimum", "first_maximum", "second_minimum",
    "second_maximum", "main_minimum")

function params(; a=P0.imc_a, K_lambda5=12.36)
    return MagneticReference.Parameters(imc_a=a, K=K_lambda5/P0.cutoff^5)
end

function controls(; order=96, n_max=64, pz_max=40.0)
    return MagneticReference.Controls(order=order, n_max=n_max,
        pz_max=pz_max, route=:landau_legacy)
end

function evaluate_variant(x, muq, eB, p, ctl;
    polyakov_zero=false, vacuum_level_cutoff=false, vacuum_extra_levels=0)
    ev = MagneticReference.evaluate(x, 5/H, fill(muq,3), eB, p, ctl)
    gap = copy(ev.condensate_gap)
    omega = ev.omega
    if vacuum_level_cutoff
        # A controlled finite-sum convention, not the infinite smooth sum
        # printed in Eq. (5-7). Matter levels remain uncut.
        for f in 1:3
            qB = MagneticReference.CHARGES[f]*eB
            ncut = floor(Int,p.cutoff^2/(2qB))+vacuum_extra_levels
            cutctl = controls(order=ctl.order, n_max=ncut, pz_max=ctl.pz_max)
            v0,g0 = MagneticReference.vacuum(ev.masses[f],qB,p,ctl)
            v1,g1 = MagneticReference.vacuum(ev.masses[f],qB,p,cutctl)
            gap[f] += g1-g0
            omega += v1-v0
        end
    end
    residual = vcat(gap,polyakov_zero ? x[4:5] : ev.residual[4:5])
    return merge(ev,(;residual,condensate_gap=gap,omega))
end

function fd_jacobian(f,x,h)
    out = Matrix{Float64}(undef,length(f(x)),length(x))
    for j in eachindex(x)
        step=h*max(1.0,abs(x[j]))
        xp,xm=copy(x),copy(x)
        xp[j]+=step; xm[j]-=step
        out[:,j]=(f(xp)-f(xm))/(2step)
    end
    return out
end

function solve_residual(f,seed; h=5e-5, iterations=90)
    f!(r,x)=(r.=f(x))
    j!(J,x)=(J.=fd_jacobian(f,x,h))
    sol=nlsolve(f!,j!,copy(seed);method=:trust_region,ftol=2e-11,
        xtol=1e-12,iterations=iterations)
    residual=norm(f(sol.zero))
    return (;x=collect(sol.zero),residual,iterations=sol.iterations,
        converged=isfinite(residual)&&residual<2e-7)
end

function refine_fold(seed,eB,p,ctl; fd_step=2e-5,
    polyakov_zero=false,vacuum_level_cutoff=false,vacuum_extra_levels=0,
    check_neighbors=true)
    eval=(x,muq)->evaluate_variant(x,muq,eB,p,ctl;
        polyakov_zero,vacuum_level_cutoff,vacuum_extra_levels)
    F=(x,muq)->eval(x,muq).residual
    J=(x,muq)->fd_jacobian(y->F(y,muq),x,fd_step)
    s=svdvals(J(seed[1:5],seed[6]))
    determinant_scale=prod(s[1:4])
    augmented=z->vcat(F(z[1:5],z[6]),det(J(z[1:5],z[6]))/determinant_scale)
    sol=solve_residual(augmented,seed)
    x,muq=sol.x[1:5],sol.x[6]
    ev=eval(x,muq)
    jac=J(x,muq)
    decomp=svd(jac)
    neighbors=Any[]
    if sol.converged && check_neighbors
        direction=decomp.V[:,end]
        index=argmax(abs.(direction))
        for delta in (-1e-3,1e-3)
            val=x[index]+delta
            neighbor_eq=z->vcat(F(z[1:5],z[6]),z[index]-val)
            initial=copy(sol.x)
            initial[1:5].+=delta.*direction./direction[index]
            r=solve_residual(neighbor_eq,initial;h=2e-5)
            en=eval(r.x[1:5],r.x[6])
            push!(neighbors,(;delta_coordinate=delta,index,converged=r.converged,
                muB_delta_MeV=3H*(r.x[6]-muq),
                rhoB_delta=sum(en.rho-ev.rho)/(3*.16),residual=r.residual))
        end
    end
    return (;converged=sol.converged,augmented_residual=sol.residual,
        gap_residual=norm(ev.residual),gap_sigma_min=minimum(decomp.S),
        mass_jacobian_sigma_min=ev.mass_jacobian_sigma_min,
        muB_MeV=3H*muq,rhoB_over_rho0=sum(ev.rho)/(3*.16),
        masses_MeV=ev.masses.*H,x_state=x,neighbors)
end

fold_seed(r)=vcat(r.x_state,r.muB_MeV/(3H))

function parameter_transition(seed,eB,p0,p1,ctl;steps=6,kwargs...)
    last=nothing
    z=copy(seed)
    for k in 1:steps
        t=k/steps
        p=params(a=(1-t)*p0.imc_a+t*p1.imc_a,
            K_lambda5=((1-t)*p0.K+t*p1.K)*p0.cutoff^5)
        last=refine_fold(z,eB,p,ctl;check_neighbors=k==steps,kwargs...)
        last.converged || return last
        z=fold_seed(last)
    end
    return last
end

function vacuum_anchor(c,p,ctl;vacuum_level_cutoff=false,vacuum_extra_levels=0)
    eB=c["eB_GeV2"]*(1000/H)^2
    point=c["points"][argmin([abs(v["muB_MeV"]-500) for v in c["points"]])]
    x=Float64.(point["x_state"])
    muq=500/(3H)
    function eq(v)
        evaluate_variant(v,muq,eB,p,ctl;vacuum_level_cutoff,vacuum_extra_levels).residual
    end
    sol=solve_residual(eq,x;iterations=160)
    ev=evaluate_variant(sol.x,muq,eB,p,ctl;vacuum_level_cutoff,vacuum_extra_levels)
    return (;eB_GeV2=c["eB_GeV2"],converged=sol.converged,
        residual=sol.residual,masses_MeV=ev.masses.*H,
        G_over_G0=MagneticReference.coupling(eB,p)/p.G0,
        x_state=sol.x)
end

function main(args=ARGS)
    length(args)==2 || throw(ArgumentError("smooth PALC input and output required"))
    input,output=args
    cases=JSON.parsefile(input)["cases"]
    c=only(c for c in cases if c["eB_GeV2"]==0.05)
    eB=.05*(1000/H)^2
    selected=c["folds"][2:6]
    variants=[
        (;id="o48_n32",ctl=controls(order=48,n_max=32),h=2e-5),
        (;id="o96_n32",ctl=controls(order=96,n_max=32),h=2e-5),
        (;id="o96_n64",ctl=controls(order=96,n_max=64),h=2e-5),
        (;id="o128_n96",ctl=controls(order=128,n_max=96),h=2e-5),
        (;id="o128_n96_pz60",ctl=controls(order=128,n_max=96,pz_max=60),h=2e-5),
        (;id="o128_n96_fd5e6",ctl=controls(order=128,n_max=96),h=5e-6)]
    out=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "source"=>input,"T_MeV"=>5.0,"eB_GeV2"=>.05,
        "method"=>"F=0 and det(F_x)=0; finite-difference steps and two-sided neighbors",
        "folds"=>Any[],"anchors"=>Any[])
    persist()=open(io->JSON.print(io,out,2),output,"w")
    mkpath(dirname(output))
    for (id,point) in zip(FEATURE_IDS,selected)
        println("fold ",id);flush(stdout)
        seed=vcat(Float64.(point["x_state"]),point["muB_MeV"]/(3H))
        numerical=Any[]
        for v in variants
            r=refine_fold(seed,eB,P0,v.ctl;fd_step=v.h)
            push!(numerical,(;id=v.id,order=v.ctl.order,n_max=v.ctl.n_max,
                pz_max=v.ctl.pz_max,fd_step=v.h,result=r))
            println(" ",v.id," ",r.converged," muB=",r.muB_MeV,
                " residual=",r.augmented_residual);flush(stdout)
            r.converged && (seed=fold_seed(r))
        end
        convention=Any[]
        ctl=controls(order=96,n_max=64)
        for (name,p) in (("thesis_K12p56",params(K_lambda5=12.56)),
            ("printed_a0p108805",params(a=.108805)))
            r=parameter_transition(seed,eB,P0,p,ctl)
            push!(convention,(;id=name,a=p.imc_a,K_lambda5=p.K*p.cutoff^5,result=r))
            println(" ",name," ",r.converged," muB=",r.muB_MeV);flush(stdout)
        end
        for (name,pz,vc) in (("polyakov_zero",true,false),
            ("vacuum_levels_below_Lambda",false,true))
            r=refine_fold(seed,eB,P0,ctl;polyakov_zero=pz,vacuum_level_cutoff=vc)
            push!(convention,(;id=name,result=r))
            println(" ",name," ",r.converged," muB=",r.muB_MeV);flush(stdout)
        end
        push!(out["folds"],(;id,initial_muB_MeV=point["muB_MeV"],numerical,convention))
        persist()
    end
    for c in cases, (name,p) in (("source_parameters",P0),
        ("thesis_K12p56",params(K_lambda5=12.56)),("printed_a0p108805",params(a=.108805)))
        println("anchor ",c["eB_GeV2"]," ",name);flush(stdout)
        r=vacuum_anchor(c,p,controls(order=128,n_max=96))
        push!(out["anchors"],(;id=name,a=p.imc_a,K_lambda5=p.K*p.cutoff^5,result=r))
        persist()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
