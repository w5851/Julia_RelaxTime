#!/usr/bin/env julia
# Targeted closure after the correctness fixes; independent literature targets
# remain outputs of crosscheck_literature.jl / crosscheck_fu2013.jl.
using JSON, Dates, LinearAlgebra, StaticArrays, ForwardDiff
include(joinpath(@__DIR__,"crosscheck_fu2013.jl"))
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(ROOT,"src","models","Models.jl"))
const MT=Models.MagneticThermodynamics
const MI=MT.MagneticIntegrals

function production_params(p;polynomial=false)
    Models.PNJLCore.PNJLParams(Λ_inv_fm=p.cutoff,m_ud0_inv_fm=p.ml,m_s0_inv_fm=p.ms,
        G_fm2=p.G0,K_fm5=p.K,T0_inv_fm=p.T0,
        polyakov_scheme=polynomial ? :poly : :log,
        a0=polynomial ? 6.75 : p.a0,a1=polynomial ? -1.95 : p.a1,
        a2=polynomial ? 2.625 : p.a2,a3=polynomial ? -7.44 : 0.,
        b3=polynomial ? .75 : p.b3,b4=polynomial ? 7.5 : 7.555)
end

function model_for(eB,p;polynomial=false,route=:mfir,p_num=128,nmax=240,zeta_num=256)
    pa=production_params(p;polynomial)
    imc=MT.MagneticIMCParams(a=p.imc_a,b=p.imc_b,c=p.imc_c,d=p.imc_d,
        Λ_QCD_MeV=p.qcd_scale*H)
    conf=MT.MagneticConfig(;eB_fm2=eB,n_max=nmax,p_num=p_num,pz_max=40.,zeta_num,
        route,imc,params=pa)
    Models.PNJLMagneticModel(Models.PNJLModel(),conf)
end

function compare_at(x,T,mu,eB,p,ctl;polynomial=false,p_num=128,zeta_num=256)
    model=model_for(eB,p;polynomial,route=ctl.route,p_num,nmax=ctl.n_max,zeta_num)
    st=SVector{5,Float64}(x);mv=SVector{3,Float64}(mu)
    prod=MT.calculate_magnetic_omega_components(st,mv,T,model.magnetic)
    gap=Models.magnetic_condensate_residual(model,st,T,mv)
    gradient=Models.magnetic_gap_residual(model,st,T,mv)
    rho=MT.calculate_magnetic_rho(st,mv,T,model.magnetic)
    if polynomial
        ref=evaluate_fu(x,T,eB,p,ctl)
        rgap=ref.residual[1:3];rgrad=vcat(MR.mass_jacobian(x,eB,p)'*rgap,ref.residual[4:5])
        rrho=zeros(3)
    else
        ref=MR.evaluate(x,T,mu,eB,p,ctl)
        rgap=ref.condensate_gap;rgrad=ref.stationarity;rrho=ref.rho
    end
    (;p_num,zeta_num,n_max=ctl.n_max,omega_delta=prod.omega-ref.omega,
        gap_delta=norm(gap-rgap),gradient_delta=norm(gradient-rgrad),rho_delta=norm(rho-rrho),
        gap_norm=norm(gap),reference_gap_norm=norm(rgap),
        chain_rule_error=norm(gradient[1:3]-MR.mass_jacobian(x,eB,p)'*gap))
end

function main_corrected(folder,output)
    out=Dict{String,Any}("generated_at"=>string(now()),"julia_version"=>string(VERSION),
        "states"=>Any[],"roots"=>Any[],"zeta"=>Any[])
    persist()=open(io->JSON.print(io,out,2),output,"w")
    old=JSON.parsefile(joinpath(folder,"kernel_closure_audit.json"))
    for item in old["states"][1:6]
        s=item["input"];x=Float64.(s["x_state"]);T=s["T_MeV"]/H
        mu=fill(s["muB_MeV"]/(3H),3);eb=s["eB_GeV2"]*(1000/H)^2
        route=Symbol(s["route"]);p=MR.Parameters()
        ctl=MR.Controls(order=128,n_max=96,route=route)
        rows=Any[]
        for order in (32,64,128,256)
            r=compare_at(x,T,mu,eb,p,ctl;p_num=order)
            push!(rows,r)
            println(s["id"]," p=",order," gap=",r.gap_delta," rho=",r.rho_delta);flush(stdout)
        end
        push!(out["states"],(;source=s,checks=rows));persist()
    end
    ceps=JSON.parsefile(joinpath(folder,"literature_cep.json"))
    for c in ceps["cases"]
        t=c["target"];r=last(c["checks"]);p=literature_params(imc=t["imc"])
        ctl=MR.Controls(order=96,n_max=300)
        x=Float64.(r["x_state"]);eb=t["eB"]*(1000/H)^2
        rows=[compare_at(x,r["T_MeV"]/H,fill(r["muq"],3),eb,p,ctl;p_num=n) for n in (64,128)]
        push!(out["states"],(;source=t,id="CEP_"*t["id"],checks=rows))
        println("CEP ",t["id"]," deltas=",rows[end]);flush(stdout);persist()
    end
    fu=JSON.parsefile(joinpath(folder,"literature_fu2013.json"))
    c=last(fu["cases"]);r=last(c["checks"]);p=literature_params()
    eb=c["eB_GeV2"]*(1000/H)^2;x=Float64.(r["x_state"])
    ctl=MR.Controls(order=96,n_max=300,route=:landau_legacy)
    rows=[compare_at(x,r["T_MeV"]/H,zeros(3),eb,p,ctl;polynomial=true,p_num=n) for n in (128,256,512)]
    push!(out["states"],(;source="Fu2013 20mpi2 peak",checks=rows));persist()
    qB=.4*(1000/H)^2*2/3
    for xz in (1e-8,7.60072e-5,.01,.9999,1.0001,10.)
        mass=sqrt(2qB*xz)
        expected=-3mass*qB/(2pi^2)*(MR.loggamma(xz)-log(2pi)/2+xz-(xz-.5)*log(xz))
        values=Any[]
        for n in (32,64,128,256)
            g=ForwardDiff.derivative(m->MI.omega_magnetic_mfir(m,2/3,qB/(2/3);zeta_num=n),mass)
            push!(values,(;zeta_num=n,derivative_delta=g-expected))
        end
        push!(out["zeta"],(;x=xz,mass_MeV=mass*H,values));persist()
    end
    # Reuse the complete previous target root sets. Correct with the actual
    # public solver and independently check every returned state and density.
    for (file,ebG) in (("resolved_mfir_o48_fine_v2.json",.155),)
        ref=JSON.parsefile(joinpath(folder,file))
        for c in ref["cases"]
            c["eB_GeV2"] in (.155,.4) || continue
            eb=c["eB_GeV2"]*(1000/H)^2;p=MR.Parameters()
            model=model_for(eb,p;p_num=64,nmax=96)
            for target in c["targets"]
                muB=target["muB_MeV"]
                muB in (900.,1000.,1050.,1100.,1150.,1200.,1500.) || continue
                mu=fill(muB/(3H),3)
                checks=Any[]
                for root in target["palc_roots"]
                    seed=Float64.(root["x_state"])
                    sol=Models.solve_magnetic_gap(model,5/H,mu;initial_guess=seed,
                        include_default_seeds=false,fallback_method=nothing,ftol=1e-10,
                        residual_norm_max=1e-8)
                    if sol.converged
                        state=sol.candidates[sol.selected_index].x_state
                        e=MR.evaluate(collect(state),5/H,mu,eb,p,MR.Controls(order=96,n_max=128))
                        push!(checks,(;converged=true,state_delta=norm(state-seed),
                            independent_gap_norm=norm(e.residual),x_state=collect(state)))
                    else
                        push!(checks,(;converged=false))
                    end
                end
                push!(out["roots"],(;eB_GeV2=c["eB_GeV2"],muB_MeV=muB,checks))
                println("roots ",c["eB_GeV2"]," ",muB," success=",count(q->q.converged,checks),"/",length(checks));flush(stdout)
                persist()
            end
        end
    end
end

function check_default_solver(folder,output)
    result=JSON.parsefile(output)
    source=JSON.parsefile(joinpath(folder,"corrected_palc_p128.json"))
    jobs=Any[]
    # Three fields, the seven-root slice, and high-density light/strange roots.
    # Perturb every input so success is not a zero-iteration reference replay.
    for c in source["cases"]
        targets=c["eB_GeV2"]==.155 ? (1000.,1500.) : (1050.,1500.)
        for t in c["targets"]
            t["muB_MeV"] in targets || continue
            for r in t["palc_roots"]
                push!(jobs,(;eB_GeV2=c["eB_GeV2"],T_MeV=c["T_MeV"],
                    muB_MeV=t["muB_MeV"],x=Float64.(r["x_state"]),
                    params=MR.Parameters(),id="cold_branch"))
            end
        end
    end
    # Stable high-T equilibrium, separately from the singular CEP equations.
    for c in JSON.parsefile(joinpath(folder,"external_anchor_drift.json"))["anchors"]
        suffix=last(split(c["record_id"],'_'))
        ebG=Dict("eb02"=>.2,"eb04"=>.4,"eb08"=>.8)[suffix]
        # External hbarc differs; these old states are used only as starting
        # guesses. Re-solve with the project's own units and configuration.
        push!(jobs,(;eB_GeV2=ebG,T_MeV=240.,muB_MeV=0.,
            x=Float64.(c["reference_state"]),params=MR.Parameters(),id="hot_anchor"))
    end
    rows=Any[]
    for job in jobs
        eb=job.eB_GeV2*(1000/H)^2;T=job.T_MeV/H;mu=fill(job.muB_MeV/(3H),3)
        model=model_for(eb,job.params;p_num=128,nmax=nothing)
        seed=job.x+[.002,-.0015,.001,.00002,-.00001]
        seed[4:5]=max.(seed[4:5],1e-8)
        sol=Models.solve_magnetic_gap(model,T,mu;initial_guess=seed,
            include_default_seeds=false,ftol=1e-11,xtol=1e-12,iterations=150,
            residual_norm_max=1e-8)
        candidate=sol.candidates[sol.selected_index]
        state=collect(candidate.x_state)
        reference=MR.evaluate(state,T,mu,eb,job.params,
            MR.Controls(order=96,n_max=max(128,candidate.n_max+32)))
        gap=Models.magnetic_condensate_residual(model,state,T,mu;n_max=candidate.n_max)
        prho=MT.calculate_magnetic_rho(state,mu,T,model.magnetic)
        row=(;id=job.id,eB_GeV2=job.eB_GeV2,T_MeV=job.T_MeV,muB_MeV=job.muB_MeV,
            n_max=candidate.n_max,converged=sol.converged,physical=candidate.physical,
            iterations=candidate.iterations,seed_distance=norm(state-seed),
            previous_state_distance=norm(state-job.x),
            condensate_gap=norm(gap),independent_residual=norm(reference.residual),
            density_delta=norm(prho-reference.rho),x_state=state)
        push!(rows,row)
        println("default solver ",job.id," eB=",job.eB_GeV2," muB=",job.muB_MeV,
            " nmax=",row.n_max," steps=",row.iterations," gap=",row.independent_residual);flush(stdout)
    end
    result["default_solver"]=rows
    open(io->JSON.print(io,result,2),output,"w")
end

if abspath(PROGRAM_FILE)==@__FILE__
    if length(ARGS)==3 && ARGS[3]=="--consumer"
        check_default_solver(ARGS[1],ARGS[2])
    else
        main_corrected(ARGS[1],ARGS[2])
    end
end
