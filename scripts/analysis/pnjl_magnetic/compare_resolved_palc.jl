#!/usr/bin/env julia
#= A single anchor PALC probe with independent Fermi-resolved quadrature.
   equation=gap traces condensate self consistency; equation=stationarity
   traces the Omega gradient. The two agree when dM/dphi is nonsingular.
   No production model, default or acceptance baseline is changed. =#

using Dates, JSON
using LinearAlgebra: norm
using NLsolve
import BifurcationKit
import BifurcationKit: @optic
include(joinpath(@__DIR__,"MagneticReference.jl"))
using .MagneticReference
const H = MagneticReference.HBARC

function jacobian_fd(f,x,p)
    J=Matrix{Float64}(undef,length(x),length(x))
    for j in eachindex(x)
        h=2e-5*max(1.0,abs(x[j]))
        xp,xm=copy(x),copy(x)
        xp[j]+=h;xm[j]-=h
        J[:,j]=(f(xp,p)-f(xm,p))/(2h)
    end
    return J
end

function correct(f,x,p)
    f!(r,v)=(r.=f(v,p))
    j!(J,v)=(J.=jacobian_fd(f,v,p))
    r=nlsolve(f!,j!,copy(Float64.(x));method=:trust_region,ftol=1e-11,xtol=1e-12,iterations=120)
    state=Float64.(r.zero)
    residual=norm(f(state,p))
    return (;x=state,converged=residual<1e-7,residual)
end

function branch_seeds(points,muq)
    seeds=Vector{Float64}[]
    for pt in points
        abs(pt.param-muq)<1e-11 && push!(seeds,copy(pt.solution))
    end
    for i in 2:length(points)
        a,b=points[i-1],points[i]
        min(a.param,b.param)<muq<max(a.param,b.param) || continue
        w=(muq-a.param)/(b.param-a.param)
        push!(seeds,(1-w).*a.solution+w.*b.solution)
    end
    return seeds
end

function root_records(f,seeds,p,T,eB,params,ctl)
    roots=Any[]
    failures=Any[]
    for seed in seeds
        try
            r=correct(f,seed,p)
            r.converged || (push!(failures,(;residual=r.residual));continue)
            any(norm(r.x-item.x_state)<1e-5 for item in roots) && continue
            ev=MagneticReference.evaluate(r.x,T,fill(p.muq,3),eB,params,ctl)
            push!(roots,(;x_state=r.x,residual_norm=r.residual,
                condensate_gap_norm=norm(ev.condensate_gap),stationarity_norm=norm(ev.stationarity),
                omega=ev.omega,masses_MeV=ev.masses.*H,rho=ev.rho,
                rhoB_over_rho0=sum(ev.rho)/(3*0.16),
                physical=all(>(0),ev.masses)&&all(v->-1e-10<=v<=1,r.x[4:5]),
                mass_jacobian_sigma_min=ev.mass_jacobian_sigma_min))
        catch err
            push!(failures,(;error=sprint(showerror,err)))
        end
    end
    return roots,failures
end

function run_case(ebG,opt)
    T=parse(Float64,get(opt,"T","5"))/H
    eb=ebG*(1000/H)^2
    route=Symbol(get(opt,"route","mfir"))
    equation=Symbol(get(opt,"equation","gap"))
    equation in (:gap,:stationarity) || throw(ArgumentError("unknown equation"))
    params=MagneticReference.Parameters(imc_a=parse(Float64,get(opt,"a","0.0108805")))
    ctl=MagneticReference.Controls(order=parse(Int,get(opt,"order","24")),route=route,
        n_max=parse(Int,get(opt,"n-max","32")))
    f=(x,p)->begin
        try
            ev=MagneticReference.evaluate(x,T,fill(p.muq,3),eb,params,ctl)
            equation===:gap ? ev.residual : ev.stationarity
        catch err
            # A Newton trial outside the logarithmic potential's domain is
            # rejected by the residual/line search. No outside-domain state
            # can pass the final independent evaluate() and residual check.
            err isa DomainError || rethrow()
            fill(1e6,5)
        end
    end
    side=get(opt,"anchor-side","high")
    anchor_mu=parse(Float64,get(opt,"anchor-muB",side=="high" ? "1700" : "500"))
    p0=(;muq=anchor_mu/(3H))
    seed=side=="high" ? [-0.02,-0.02,-0.30,0.01,0.01] : [-1.8,-1.8,-2.2,1e-7,1e-7]
    t0=time()
    anchor=correct(f,seed,p0)
    anchor.converged || error("single seed anchor failed: $(anchor.residual)")
    minimum(MagneticReference.masses(anchor.x,eb,params))>0 || error("anchor mass negative")
    println(" anchor ",anchor_mu," x=",anchor.x);flush(stdout)
    prob=BifurcationKit.BifurcationProblem(f,anchor.x,p0,(@optic _.muq);
        J=(x,p)->jacobian_fd(f,x,p),
        record_from_solution=(x,p;kwargs...)->(;solution=copy(Float64.(x))))
    opts=BifurcationKit.ContinuationPar(p_min=500/(3H),p_max=1700/(3H),
        ds=parse(Float64,get(opt,"ds","0.002")),dsmax=parse(Float64,get(opt,"dsmax","0.006")),
        max_steps=parse(Int,get(opt,"max-steps","1400")),detect_bifurcation=0,
        newton_options=BifurcationKit.NewtonPar(tol=1e-11,max_iterations=24,verbose=false,linesearch=true))
    br=BifurcationKit.continuation(prob,BifurcationKit.PALC(),opts;bothside=true,verbosity=0)
    seconds=time()-t0
    println(" continued ",length(br.branch)," in ",seconds," sec");flush(stdout)
    points=Any[]
    for (i,pt) in enumerate(br.branch)
        ev=MagneticReference.evaluate(pt.solution,T,fill(pt.param,3),eb,params,ctl)
        push!(points,(;index=i,muB_MeV=pt.param*3H,x_state=pt.solution,
            omega=ev.omega,rho=ev.rho,rhoB_over_rho0=sum(ev.rho)/(3*0.16),
            level_rho=ev.level_rho,masses_MeV=ev.masses.*H,
            condensate_gap_norm=norm(ev.condensate_gap),stationarity_norm=norm(ev.stationarity),
            mass_jacobian_sigma_min=ev.mass_jacobian_sigma_min,
            physical=all(>(0),ev.masses)&&all(v->-1e-10<=v<=1,pt.solution[4:5])))
    end
    folds=Any[]
    for i in 2:length(points)-1
        (points[i].muB_MeV-points[i-1].muB_MeV)*(points[i+1].muB_MeV-points[i].muB_MeV)<0 || continue
        push!(folds,points[i])
    end
    targets=Any[]
    target_values=parse.(Float64,split(get(opt,"targets","900,1000,1050,1100,1150,1200,1500,1650"),','))
    for muB in target_values
        p=(;muq=muB/(3H))
        seeds=branch_seeds(br.branch,p.muq)
        proots,pfail=root_records(f,seeds,p,T,eb,params,ctl)
        catalogue=[[-1.8,-1.8,-2.2,1e-6,1e-6],[-0.02,-0.02,-0.3,0.01,0.01],
            [-1.1,-1.1,-2.1,0.003,0.003],[-0.3,-0.5,-2.05,0.004,0.004],
            [-0.02,-0.02,-1.2,0.008,0.008],[-0.02,-0.02,-0.7,0.008,0.008]]
        croots,cfail=root_records(f,catalogue,p,T,eb,params,ctl)
        distances=[isempty(proots) ? Inf : minimum(norm(c.x_state-r.x_state) for r in proots) for c in croots]
        push!(targets,(;muB_MeV=muB,palc_roots=proots,catalog_roots=croots,
            palc_failures=pfail,catalog_failures=cfail,catalog_min_state_distances=distances,
            matched_catalog_roots=count(d->d<1e-5,distances)))
        println(" target ",muB," roots ",length(proots),"/",length(croots));flush(stdout)
    end
    return (;eB_GeV2=ebG,T_MeV=T*H,route=String(route),equation=String(equation),
        imc_a=params.imc_a,order=ctl.order,n_max=ctl.n_max,pz_max=ctl.pz_max,
        anchor_muB_MeV=anchor_mu,anchor_state=anchor.x,palc_seconds=seconds,total_seconds=time()-t0,
        points,folds,targets)
end

function main(args=ARGS)
    opt=Dict{String,String}()
    output="results/resolved_palc.json"
    for a in args
        if startswith(a,"--")
            k,v=split(a[3:end],'=';limit=2)
            opt[k]=v
        else
            output=a
        end
    end
    cases=Any[]
    mkpath(dirname(output))
    for eb in parse.(Float64,split(get(opt,"eB","0.05,0.155,0.4"),','))
        println("running resolved eB=",eb," route=",get(opt,"route","mfir"));flush(stdout)
        try
            push!(cases,run_case(eb,opt))
        catch err
            err isa InterruptException && rethrow()
            push!(cases,(;eB_GeV2=eb,error=sprint(showerror,err)))
            println("error ",sprint(showerror,err));flush(stdout)
        end
        open(output,"w") do io
            JSON.print(io,(;generated_at=string(Dates.now()),julia_version=string(VERSION),cases),2)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
