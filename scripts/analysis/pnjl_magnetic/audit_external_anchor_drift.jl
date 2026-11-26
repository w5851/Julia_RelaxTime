#!/usr/bin/env julia
include(joinpath(@__DIR__,"check_corrected_entrypoint.jl"))

function read_rows(path)
    lines=readlines(path);header=split(lines[1],',')
    [Dict(zip(header,split(line,','))) for line in lines[2:end] if !isempty(strip(line))]
end

function external_params(row)
    hc=parse(Float64,row["external_hc_MeV_fm"]);L=602.3/hc
    MR.Parameters(cutoff=L,G0=1.835/L^2,K=12.36/L^5,ml=5.5/hc,ms=140.7/hc,
        T0=210/hc,qcd_scale=300/hc)
end

function old_zeta(x,n)
    nodes,weights=MI.zeta_nodes(n)
    scalar=x*(-x+2(x-1)*log(x))/4
    acc=zero(x)
    for (u,w) in zip(nodes,weights)
        t=-log1p(-u)
        acc+=w*(2x*MI._atan_nonnegative(t/x)+t*log(x*x+t*t))/
            ((1-u)*2*(exp(2pi*t)-1))
    end
    scalar+2acc
end

function historical_omega(x,T,mu,eB,p,pnum,znum,nmax)
    M=MR.masses(x,eB,p);G=MR.coupling(eB,p)
    vac=0.;therm=0.
    for f in 1:3
        m=abs(M[f])+1e-12;qB=MR.CHARGES[f]*eB;xx=m*m/(2qB)
        vac+=-2p.colors*Models.PNJLCore.vacuum_integral_with_cutoff(M[f],convert(typeof(M[f]),p.cutoff))-
            p.colors*qB^2/(2pi^2)*(old_zeta(xx,znum)-(xx*xx-xx)*log(xx)/2+xx*xx/4)
        nodes,weights=MI.pz_nodes(pnum)
        for n in 0:nmax
            therm-=T*qB/(2pi^2)*(n==0 ? 1 : 2)*40*sum(
                weights[j]*MI._log_polyakov_pair(sqrt(M[f]^2+2n*qB+(40nodes[j])^2),
                    mu[f],T,x[4],x[5]) for j in eachindex(nodes))
        end
    end
    poly=MR.polyakov(x[4],x[5],T,p)[1]
    2G*sum(abs2,x[1:3])-4p.K*prod(x[1:3])+poly+vac+therm
end

function main_external(output)
    folder=joinpath(ROOT,"tests","validation","data","targets","pnjl","reference")
    kernel=only(read_rows(joinpath(folder,"pnjl_magnetic_external_kernel_targets_v1.csv")))
    anchors=read_rows(joinpath(folder,"pnjl_magnetic_external_equilibrium_anchor_targets_v1.csv"))
    out=Dict{String,Any}("generated_at"=>string(now()),"source"=>"unchanged external pnjl_mag fixtures; old discretization reconstructed only for attribution","kernel"=>nothing,"anchors"=>Any[])
    for (i,row) in enumerate(vcat([kernel],anchors))
        p=external_params(row);hc=parse(Float64,row["external_hc_MeV_fm"])
        T=parse(Float64,row["T_MeV"])/hc;mu=fill(parse(Float64,row["muB_MeV"])/(3hc),3)
        eB=parse(Float64,row["eB_GeV2"])*(1000/hc)^2
        x=[parse(Float64,row[s]) for s in ("phi_u","phi_d","phi_s","Phi","PhiBar")]
        nmax=i==1 ? parse(Int,row["landau_levels"])-1 : parse(Int,row["n_max"])
        expected=parse(Float64,row[i==1 ? "expected_omega_fm4" : "omega_fm4"])
        n=parse(Int,row["p_num"]);z=parse(Int,row["zeta_num"])
        old=historical_omega(x,T,mu,eB,p,n,z,nmax)
        ref=MR.evaluate(x,T,mu,eB,p,MR.Controls(order=128,n_max=nmax))
        rows=Any[]
        for order in (16,32,64,128,256,512)
            oo=historical_omega(x,T,mu,eB,p,order,order<128 ? z : 1024,nmax)
            push!(rows,(;p_num=order,old_omega=oo,delta_to_reference=oo-ref.omega))
        end
        if i==1
            out["kernel"]=(;record_id=row["record_id"],expected,old_omega=old,
                reference_omega=ref.omega,old_minus_expected=old-expected,
                reference_minus_expected=ref.omega-expected,checks=rows)
        else
            eq=finite_equilibrium(T,mu[1],eB,p,MR.Controls(order=96,n_max=nmax),x)
            eqhi=finite_equilibrium(T,mu[1],eB,p,MR.Controls(order=160,n_max=nmax),eq.x)
            refinement=Any[]
            for (pn,zn) in ((128,256),(256,256),(128,1024),(256,1024))
                fn=y->ForwardDiff.gradient(v->historical_omega(v,T,mu,eB,p,pn,zn,nmax),y)
                rr=solve_fd(fn,x;h=2e-5,ftol=2e-12)
                push!(refinement,(;p_num=pn,zeta_num=zn,converged=rr.converged,
                    residual=rr.residual,state_delta=norm(rr.x-x),
                    reference_delta=norm(rr.x-eq.x)))
            end
            model=model_for(eB,p;p_num=128,nmax)
            # Map qcd_scale to the project's global hbarc used by coupling_GB.
            conf=model.magnetic
            imc=MT.MagneticIMCParams(a=p.imc_a,b=p.imc_b,c=p.imc_c,d=p.imc_d,
                Λ_QCD_MeV=300*Main.Constants_PNJL.ħc_MeV_fm/hc)
            conf=MT.MagneticConfig(eB_fm2=eB,n_max=nmax,p_num=128,pz_max=40.,
                zeta_num=256,params=conf.params,imc=imc)
            model=Models.PNJLMagneticModel(Models.PNJLModel(),conf)
            solved=Models.solve_magnetic_gap(model,T,mu;initial_guess=x,include_default_seeds=false,
                fallback_method=nothing,ftol=1e-11,residual_norm_max=1e-8)
            state=collect(solved.candidates[solved.selected_index].x_state)
            push!(out["anchors"],(;record_id=row["record_id"],expected,old_omega=old,
                source_independent_gap=norm(ref.residual),reference_state=eq.x,
                reference_state_delta=norm(eq.x-x),reference_refinement_delta=norm(eq.x-eqhi.x),
                production_state=state,production_state_delta=norm(state-x),
                production_reference_delta=norm(state-eq.x),checks=rows,refinement))
        end
        println(row["record_id"]," old delta=",old-expected," ref delta=",ref.omega-expected);flush(stdout)
        open(io->JSON.print(io,out,2),output,"w")
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    main_external(only(ARGS))
end
