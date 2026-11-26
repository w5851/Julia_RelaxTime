#!/usr/bin/env julia

#= Targeted, independent literature checks for magnetized three-flavor matter.

Grunfeld et al., PRD 90, 044024 (2014), Eqs. (3)--(7), Table I:
exact T=0 MFIR loops and pressure equality, including weak cascaded transitions.
Costa et al., PRD 92, 036012 (2015), Sec. III B: two T=0 coexistences.
Ferreira et al., PRD 97, 014014 (2018), Tables I--II: finite-T CEPs.
No production model/kernel is used by the reference calculations.
=#
using JSON, Dates, LinearAlgebra, NLsolve
include(joinpath(@__DIR__, "MagneticReference.jl"))
const MR = MagneticReference
const H = MR.HBARC

function literature_params(; set=2, imc=false, T0=270.0)
    L = (set == 1 ? 631.4 : 602.3)/H
    return MR.Parameters(cutoff=L, G0=1.835/L^2,
        K=(set == 1 ? 9.29 : 12.36)/L^5, ms=(set == 1 ? 135.7 : 140.7)/H,
        T0=T0/H, imc_a=imc ? .0108805 : 0., imc_b=imc ? -1.0133e-4 : 0.,
        imc_c=imc ? .02228 : 0., imc_d=imc ? 1.84558e-4 : 0.)
end

function jacobian_fd(f,x;h=2e-5)
    J=Matrix{Float64}(undef,length(f(x)),length(x))
    for j in eachindex(x)
        d=h*max(1.,abs(x[j])); xp=copy(x); xm=copy(x)
        xp[j]+=d; xm[j]-=d
        J[:,j]=(f(xp)-f(xm))/(2d)
    end
    J
end

function solve_fd(f,seed;h=2e-5,ftol=1e-10,iterations=120)
    f!(r,x)=(r.=f(x))
    j!(J,x)=(J.=jacobian_fd(f,x;h))
    sol=nlsolve(f!,j!,Float64.(seed);method=:trust_region,
        ftol,xtol=1e-12,iterations)
    err=norm(f(sol.zero))
    (;x=collect(sol.zero),residual=err,converged=isfinite(err)&&err<1e-7,
        iterations=sol.iterations)
end

"""Exact T=0 medium terms; chemical potential is quark mu, units fm^-1."""
function cold_eval(M,mu,eB,p; with_omega=true)
    loop=zeros(3);rho=zeros(3); gas=0.;levels=Int[]
    ctl=MR.Controls()
    for f in 1:3
        qB=MR.CHARGES[f]*eB
        if eB==0
            m=abs(M[f]);L=p.cutoff;E=sqrt(L^2+m^2)
            loop[f]=-p.colors*M[f]/(2pi^2)*(L*E-m^2*asinh(L/m))
            gas-=p.colors/(8pi^2)*(L*E*(2L^2+m^2)-m^4*asinh(L/m))
            if abs(mu)>m
                k=sqrt(mu^2-m^2);a=asinh(k/m)
                loop[f]+=p.colors*M[f]/(2pi^2)*(abs(mu)*k-m^2*a)
                rho[f]=p.colors*sign(mu)*k^3/(3pi^2)
                gas-=p.colors/(24pi^2)*(abs(mu)*k*(2mu^2-5m^2)+3m^4*a)
            end
            push!(levels,-1)
            continue
        end
        if with_omega
            vf,lf=MR.vacuum(M[f],qB,p,ctl)
            gas+=vf;loop[f]=lf
        else
            m=abs(M[f]);L=p.cutoff;x=m^2/(2qB)
            loop[f]=-p.colors*M[f]/(2pi^2)*(L*sqrt(L^2+m^2)-m^2*asinh(L/m))-
                p.colors*M[f]*qB/(2pi^2)*(MR.loggamma(x)-log(2pi)/2+x-(x-.5)*log(x))
        end
        nmax=mu^2>M[f]^2 ? floor(Int,(mu^2-M[f]^2)/(2qB)) : -1
        push!(levels,nmax)
        for n in 0:nmax
            mn2=M[f]^2+2n*qB;k=sqrt(max(mu^2-mn2,0.));a=asinh(k/sqrt(mn2))
            pref=p.colors*qB/(2pi^2)*(n==0 ? 1 : 2)
            loop[f]+=pref*M[f]*a;rho[f]+=pref*sign(mu)*k
            gas-=pref/2*(abs(mu)*k-mn2*a)
        end
    end
    G=MR.coupling(eB,p)
    target=[p.ml-4G*loop[1]+2p.K*loop[2]*loop[3],
        p.ml-4G*loop[2]+2p.K*loop[1]*loop[3],
        p.ms-4G*loop[3]+2p.K*loop[1]*loop[2]]
    omega=gas+2G*sum(abs2,loop)-4p.K*prod(loop)
    (;gap=M-target,omega,condensates=loop,rho,levels,masses=M)
end

function cold_catalogue(mu,eB,p;step_MeV=.2)
    # Grunfeld's reduction: solve fd=fs=0 at fixed Mu, then bracket all fu zeros.
    # The reduction is checked again with a finer Mu grid for selected cases.
    grid=collect(2.0:step_MeV:800.)/H
    records=Any[];roots=Any[]; seed=[350.,550.]/H; failures=0
    last=nothing
    for u in grid
        f=z->cold_eval([u,z[1],z[2]],mu,eB,p;with_omega=false).gap[2:3]
        r=solve_fd(f,seed;iterations=80,h=2e-6)
        if !r.converged || minimum(r.x)<=0
            failures+=1;last=nothing;continue
        end
        seed=r.x
        M=[u;seed];v=cold_eval(M,mu,eB,p;with_omega=false).gap[1]
        if last!==nothing && signbit(v)!=signbit(last.v)
            full=x->cold_eval(x,mu,eB,p;with_omega=false).gap
            rr=solve_fd(full,(M+last.M)/2;h=2e-7)
            if rr.converged && minimum(rr.x)>0 && !any(norm(rr.x-q.masses)<1e-6 for q in roots)
                ev=cold_eval(rr.x,mu,eB,p)
                push!(roots,merge(ev,(;residual=rr.residual)))
            end
        end
        last=(;M,v)
    end
    sort!(roots;by=r->r.masses[1])
    (;roots,failed_reductions=failures,step_MeV)
end

function coexistence(M1,M2,mu,eB,p;h=2e-6)
    function f(z)
        a=cold_eval(z[1:3],z[7],eB,p);b=cold_eval(z[4:6],z[7],eB,p)
        vcat(a.gap,b.gap,a.omega-b.omega)
    end
    r=solve_fd(f,vcat(M1,M2,mu);h,ftol=2e-11,iterations=150)
    a=cold_eval(r.x[1:3],r.x[7],eB,p);b=cold_eval(r.x[4:6],r.x[7],eB,p)
    (;converged=r.converged,residual=r.residual,muq_MeV=r.x[7]*H,
        muB_MeV=3r.x[7]*H,mass_separation_MeV=norm(a.masses-b.masses)*H,
        phases=[(;masses_MeV=e.masses.*H,condensates=e.condensates,
            rhoB_over_rho0=sum(e.rho)/(3*.16),omega=e.omega,levels=e.levels) for e in (a,b)])
end

function run_cold(output)
    out=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "method"=>"exact T=0 MFIR, independent mass gap and pressure equality",
        "cases"=>Any[],"vacuum"=>Any[])
    persist()=open(io->JSON.print(io,out,2),output,"w")
    for set in (1,2)
        p=literature_params(;set)
        r=solve_fd(M->cold_eval(M,0.,0.,p).gap,[350.,350.,550.]/H)
        push!(out["vacuum"],(;set,converged=r.converged,residual=r.residual,masses_MeV=r.x.*H))
    end
    jobs=[(;id="Grunfeld2014_set1_B001",set=1,eB_GeV2=.01,
        literature_muqs=[335.3,342.3,345.4]),
        (;id="Costa2015_set2_B009",set=2,eB_GeV2=.09,
        literature_muqs=[1061/3,1083/3])]
    for job in jobs
        p=literature_params(set=job.set);eB=job.eB_GeV2*(1000/H)^2
        targets=Any[]
        for muMeV in job.literature_muqs
            println(job.id," catalogue mu=",muMeV);flush(stdout)
            cat=cold_catalogue(muMeV/H,eB,p)
            println(" roots=",length(cat.roots)," masses=",[r.masses.*H for r in cat.roots]);flush(stdout)
            pairs=Any[]
            # Adjacent minima in mass order are separated by one maximum.
            for i in 1:max(length(cat.roots)-2,0)
                co=coexistence(cat.roots[i].masses,cat.roots[i+2].masses,muMeV/H,eB,p)
                push!(pairs,co)
                println(" coexist mu=",co.muq_MeV," sep=",co.mass_separation_MeV," ok=",co.converged);flush(stdout)
            end
            push!(targets,(;literature_muq_MeV=muMeV,catalogue=cat,coexistences=pairs))
            # Persist after each useful target, rather than lose earlier results.
            push!(out["cases"],merge(job,(;target=targets[end])))
            persist()
        end
    end
    persist()
end

function finite_equilibrium(T,mu,eB,p,ctl,seed)
    f=x->begin
        try
            MR.evaluate(x,T,fill(mu,3),eB,p,ctl).residual
        catch err
            err isa DomainError || rethrow()
            fill(1e6,5)
        end
    end
    solve_fd(f,seed;h=2e-5,iterations=150)
end

function run_cep_preview(output)
    targets=[(;id="light_B03_constant",sector="light",eB=.3,imc=false,T=192.,muB=674.,rho=3.54),
        (;id="light_B06_constant",sector="light",eB=.6,imc=false,T=214.,muB=692.,rho=7.22),
        (;id="light_B03_imc",sector="light",eB=.3,imc=true,T=177.,muB=627.,rho=2.63),
        (;id="light_B06_imc",sector="light",eB=.6,imc=true,T=171.,muB=535.,rho=3.90),
        (;id="strange_B03_constant_1",sector="strange",eB=.3,imc=false,T=62.,muB=1330.,rho=6.67),
        (;id="strange_B06_constant",sector="strange",eB=.6,imc=false,T=124.,muB=1234.,rho=12.10),
        (;id="strange_B03_imc_1",sector="strange",eB=.3,imc=true,T=48.,muB=1193.,rho=5.40),
        (;id="strange_B03_imc_2",sector="strange",eB=.3,imc=true,T=18.,muB=1539.,rho=10.75)]
    out=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "source"=>"Ferreira2018 Tables I-II; printed rounded values, not exact roots",
        "cases"=>Any[])
    mkpath(dirname(output))
    for t in targets
        eB=t.eB*(1000/H)^2;p=literature_params(imc=t.imc)
        ctl=MR.Controls(order=32,n_max=120,pz_max=40.)
        seed=t.sector=="light" ? [-1.,-.9,-2.3,.4,.4] : [-.04,-.04,-1.3,.15,.15]
        r=finite_equilibrium(t.T/H,t.muB/(3H),eB,p,ctl,seed)
        ev=MR.evaluate(r.x,t.T/H,fill(t.muB/(3H),3),eB,p,ctl)
        push!(out["cases"],(;target=t,converged=r.converged,residual=r.residual,
            x_state=r.x,masses_MeV=ev.masses.*H,rhoB_over_rho0=sum(ev.rho)/(3*.16),
            gap_sigma_min=minimum(svdvals(jacobian_fd(x->MR.evaluate(x,t.T/H,
                fill(t.muB/(3H),3),eB,p,ctl).residual,r.x)))))
        println(t.id," ok=",r.converged," masses=",ev.masses.*H," rho=",sum(ev.rho)/(3*.16));flush(stdout)
        open(io->JSON.print(io,out,2),output,"w")
    end
end

function cep_inverse_curve(s,T,eB,p,ctl,k,seed;h=2e-5)
    # Parameterize the regular solution manifold by one condensate instead of
    # mu. This stays regular at a fold and at the cusp endpoint.
    f=z->begin
        try
            vcat(MR.evaluate(z[1:5],T,fill(z[6],3),eB,p,ctl).residual,z[k]-s)
        catch err
            err isa DomainError || rethrow()
            fill(1e6,6)
        end
    end
    r=solve_fd(f,seed;h,ftol=4e-13,iterations=90)
    r.residual<1e-9 || error("inverse curve failed: $(r.residual), s=$s T=$(T*H)")
    r
end

function cep_derivatives(s,T,eB,p,ctl,k,seed;delta=.008)
    center=cep_inverse_curve(s,T,eB,p,ctl,k,seed)
    a=cep_inverse_curve(s-delta,T,eB,p,ctl,k,center.x)
    b=cep_inverse_curve(s+delta,T,eB,p,ctl,k,center.x)
    dmu=(b.x[6]-a.x[6])/(2delta)
    d2mu=(b.x[6]+a.x[6]-2center.x[6])/delta^2
    (;dmu,d2mu,center)
end

function locate_cep(c;order=32,n_max=120,delta=.008,seed=nothing,temperature=nothing)
    t=c["target"];eB=t["eB"]*(1000/H)^2;p=literature_params(imc=t["imc"])
    ctl=MR.Controls(;order,n_max,pz_max=40.)
    k=t["sector"]=="light" ? 1 : 3
    x=seed===nothing ? vcat(Float64.(c["x_state"]),t["muB"]/(3H)) : copy(seed)
    T0=temperature===nothing ? t["T"]/H : temperature
    function f(y)
        ev=cep_derivatives(y[1],y[2],eB,p,ctl,k,x;delta)
        [ev.dmu,ev.d2mu]
    end
    r=solve_fd(f,[x[k],T0];h=1e-3,ftol=1e-9,iterations=50)
    d=cep_derivatives(r.x[1],r.x[2],eB,p,ctl,k,x;delta)
    z=d.center.x
    ev=MR.evaluate(z[1:5],r.x[2],fill(z[6],3),eB,p,ctl)
    (;converged=r.converged,residual=r.residual,order,n_max,delta,
        T_MeV=r.x[2]*H,muB_MeV=3H*z[6],rhoB_over_rho0=sum(ev.rho)/(3*.16),
        x_state=z[1:5],muq=z[6],masses_MeV=ev.masses.*H,
        gap_residual=norm(ev.residual),dmu=d.dmu,d2mu=d.d2mu,
        mass_jacobian_sigma_min=ev.mass_jacobian_sigma_min)
end

function run_cep(output)
    input=joinpath(dirname(output),"literature_cep_preview.json")
    source=JSON.parsefile(input)
    out=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "source"=>"Ferreira2018 Tables I-II", "method"=>"inverse-condensate curve cusp: dmu/dphi=d2mu/dphi2=0",
        "cases"=>Any[])
    for c in source["cases"]
        t=c["target"];println("CEP ",t["id"]);flush(stdout)
        rows=Any[]
        seed=nothing;temp=nothing
        for v in ((32,120,.012),(48,180,.006),(64,240,.003))
            try
                r=locate_cep(c;order=v[1],n_max=v[2],delta=v[3],seed,temperature=temp)
                push!(rows,r)
                println(" o=",v[1]," delta=",v[3]," CEP=",r.T_MeV," ",r.muB_MeV,
                    " rho=",r.rhoB_over_rho0," residual=",r.residual);flush(stdout)
                seed=vcat(r.x_state,r.muq);temp=r.T_MeV/H
            catch err
                push!(rows,(;order=v[1],n_max=v[2],delta=v[3],error=sprint(showerror,err)))
                println(" error ",sprint(showerror,err));flush(stdout)
                break
            end
        end
        push!(out["cases"],(;target=t,checks=rows))
        open(io->JSON.print(io,out,2),output,"w")
    end
end

function run_cold_refine(output)
    input=JSON.parsefile(joinpath(dirname(output),"literature_cold.json"))
    firstcase=first(input["cases"])
    p=literature_params(set=1);eb=.01*(1000/H)^2
    firstco=only(firstcase["target"]["coexistences"])
    m=Float64.(firstco["phases"][1]["masses_MeV"])/H
    z0=[m[2],m[3],firstco["muq_MeV"]/H]
    points=Any[];failures=Any[]
    for mU in range(m[1]*H,300.;step=-.02)
        u=mU/H
        f=z->cold_eval([u,z[1],z[2]],z[3],eb,p;with_omega=false).gap
        r=solve_fd(f,z0;h=2e-7,ftol=2e-12,iterations=160)
        if !r.converged
            push!(failures,(;Mu_MeV=mU,residual=r.residual));continue
        end
        z0=r.x
        e=cold_eval([u,z0[1],z0[2]],z0[3],eb,p)
        push!(points,(;Mu_MeV=mU,muq_MeV=z0[3]*H,masses_MeV=e.masses.*H,
            omega=e.omega,levels=e.levels,gap=norm(e.gap)))
    end
    folds=[points[i] for i in 2:length(points)-1 if
        (points[i+1].muq_MeV-points[i].muq_MeV)*(points[i].muq_MeV-points[i-1].muq_MeV)<0]
    println("cold refined folds=",folds," failures=",length(failures));flush(stdout)
    crossings=Any[]
    for mu in (341.5,342.,342.2,342.3,342.4,342.5,343.)
        roots=Any[]
        for i in 2:length(points)
            a,b=points[i-1],points[i]
            min(a.muq_MeV,b.muq_MeV)<mu<max(a.muq_MeV,b.muq_MeV) || continue
            w=(mu-a.muq_MeV)/(b.muq_MeV-a.muq_MeV)
            seed=((1-w).*a.masses_MeV+w.*b.masses_MeV)/H
            r=solve_fd(x->cold_eval(x,mu/H,eb,p;with_omega=false).gap,seed;h=2e-7)
            r.converged && push!(roots,r.x)
        end
        sort!(roots;by=first)
        println("cold mu=",mu," high roots=",[q.*H for q in roots]);flush(stdout)
        if length(roots)>=3
            co=coexistence(first(roots),last(roots),mu/H,eb,p;h=2e-7)
            push!(crossings,co)
            println(" C0-C1 coexist ",co);flush(stdout)
        end
    end
    # Retain the local path needed to expose the narrow loop, not broad scans.
    open(io->JSON.print(io,(;generated_at=string(now()),julia_version=string(VERSION),
        source="Grunfeld2014 set1 eB=.01",folds,failures,crossings,points),2),output,"w")
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: crosscheck_literature.jl cold|cep-preview|cep OUTPUT.json")
    mkpath(dirname(ARGS[2]))
    if ARGS[1]=="cold"
        run_cold(ARGS[2])
    elseif ARGS[1]=="cep-preview"
        run_cep_preview(ARGS[2])
    elseif ARGS[1]=="cep"
        run_cep(ARGS[2])
    elseif ARGS[1]=="cold-refine"
        run_cold_refine(ARGS[2])
    else
        error("unknown mode")
    end
end
