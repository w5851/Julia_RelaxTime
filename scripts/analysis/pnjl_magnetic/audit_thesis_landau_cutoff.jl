#!/usr/bin/env julia
#= Test a cutoff-derived finite Landau sum against the printed infinite sum.

ceil(Lambda^2/(2qB)) is a plausible implementation choice suggested by a
sharp cutoff; it is not stated by the thesis and not adopted as the model.
The extra level is fixed before comparing its effect at the five features.
=#
include(joinpath(@__DIR__,"audit_thesis_equations.jl"))

function audit_cutoff(args=ARGS)
    length(args)==3 || throw(ArgumentError("smooth input, refined folds, output"))
    input,auditfile,output=args
    cases=JSON.parsefile(input)["cases"]
    prior=JSON.parsefile(auditfile)
    out=Dict{String,Any}("julia_version"=>string(VERSION),
        "convention"=>"vacuum n_max per flavor=floor(Lambda^2/(2qB))+1; thermal unchanged",
        "status"=>"controlled unstated convention; not evidence that author used it",
        "features"=>Any[],"anchors"=>Any[])
    ctl=controls(order=96,n_max=96)
    persist()=open(io->JSON.print(io,out,2),output,"w")
    for f in prior["folds"]
        point=only(v["result"] for v in f["numerical"] if v["id"]=="o128_n96")
        seed=vcat(Float64.(point["x_state"]),point["muB_MeV"]/(3H))
        r=refine_fold(seed,.05*(1000/H)^2,P0,ctl;
            vacuum_level_cutoff=true,vacuum_extra_levels=1)
        push!(out["features"],(;id=f["id"],result=r))
        println(f["id"]," muB=",r.muB_MeV," converged=",r.converged);flush(stdout)
        persist()
    end
    for c in cases
        r=vacuum_anchor(c,P0,ctl;vacuum_level_cutoff=true,vacuum_extra_levels=1)
        push!(out["anchors"],r)
        println("anchor ",c["eB_GeV2"]," masses=",r.masses_MeV);flush(stdout)
        persist()
    end
end

function audit_common_cutoff(input,auditfile,output)
    cases=JSON.parsefile(input)["cases"]
    prior=JSON.parsefile(auditfile)
    out=JSON.parsefile(output)
    result=Dict{String,Any}("convention"=>
        "common n_max=max_f ceil(Lambda^2/(2q_f B)); smooth vacuum tail omitted",
        "status"=>"controlled unstated convention; no parameter fitting",
        "features"=>Any[],"anchors"=>Any[])
    # The common budget is derived from Lambda and charge, not selected by
    # the degree of agreement with a published curve.
    budget(eB)=maximum(ceil(Int,P0.cutoff^2/(2q*eB)) for q in MagneticReference.CHARGES)
    eB=.05*(1000/H)^2
    ctl=controls(order=128,n_max=budget(eB))
    for f in prior["folds"]
        point=only(v["result"] for v in f["numerical"] if v["id"]=="o128_n96")
        seed=vcat(Float64.(point["x_state"]),point["muB_MeV"]/(3H))
        r=refine_fold(seed,eB,P0,ctl)
        push!(result["features"],(;id=f["id"],n_max=ctl.n_max,result=r))
        println("common ",f["id"]," muB=",r.muB_MeV);flush(stdout)
    end
    for c in cases
        ctl=controls(order=128,n_max=budget(c["eB_GeV2"]*(1000/H)^2))
        r=vacuum_anchor(c,P0,ctl)
        push!(result["anchors"],(;n_max=ctl.n_max,result=r))
        println("common anchor ",c["eB_GeV2"]," ",r.masses_MeV);flush(stdout)
    end
    out["common_cutoff_budget"]=result
    open(io->JSON.print(io,out,2),output,"w")
end

if abspath(PROGRAM_FILE)==@__FILE__
    if length(ARGS)==4 && ARGS[4]=="--common"
        audit_common_cutoff(ARGS[1:3]...)
    else
        audit_cutoff()
    end
end
