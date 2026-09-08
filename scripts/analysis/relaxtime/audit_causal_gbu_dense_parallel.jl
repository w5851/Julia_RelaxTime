"""Checkpointed dense extension reusing hash-bound accepted three-point evidence.

Run with --threads=4. Channel tasks own their results; only the coordinator
writes CSVs. This does not change quadrature, physical gates, or the primary route.
"""
module CausalGBUDenseParallel
include("audit_causal_gbu_dense_freezeout.jl")
include("causal_gbu_dense_execution.jl")
const D=CausalGBUDenseFreezeout
const E=CausalGBUDenseExecution
const R=E.R
using CSV,JSON3,LinearAlgebra
using Main.Models

function main()
    R.method_contract()
    BLAS.set_num_threads(1)
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    output=get(ENV,"GBU_DENSE_OUTPUT",joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3"))
    energies=D.energy_grid(get(ENV,"GBU_DENSE_ENERGIES",join(D.ENERGIES,',')))
    original=joinpath(base,"method_v1_freezeout_v3")
    dd=joinpath(base,"method_v1_saved_sparse")
    rd=joinpath(base,"method_v1_freezeout_q0_reference_v3")
    inputs=[joinpath(original,"backgrounds.csv"),joinpath(dd,"manifest.json"),joinpath(rd,"manifest.json"),
        joinpath(dd,"refined_densities.csv"),joinpath(dd,"refined_shells.csv"),
        joinpath(rd,"reference_densities.csv"),joinpath(rd,"reference_shells.csv")]
    input_hashes=Dict(p=>R.hashfile(p) for p in inputs)
    dm=JSON3.read(read(joinpath(dd,"manifest.json"),String))
    for (p,h) in pairs(dm.source_hashes)
        path=String(p)
        startswith(path,"src") || continue
        R.hashfile(joinpath(R.ROOT,path))==String(h) || error("retained numeric source mismatch: $(path)")
    end
    settings=R.Settings(mesh=512,ne=128,nw=4800,nq=24,lower=1e-5,upper=56.,thermal=24.)
    for (k,v) in pairs(dm.settings)
        getproperty(settings,Symbol(k))==v || error("retained settings mismatch: $(k)")
    end
    hashes=R.start_output(output)
    direct,reference,ds,rs,bgs,failures,timings=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    for (target,path,route) in ((direct,inputs[4],"direct_finite_q"),(reference,inputs[6],"q0_lambda_reference"))
        for r in CSV.File(path)
            r.sqrt_s_NN_GeV in energies || continue
            push!(target,(sqrt_s_NN_GeV=r.sqrt_s_NN_GeV,T_MeV=r.T_MeV,muB_MeV=r.muB_MeV,background_p=48,
                channel=String(r.channel),route=route,density_inv_fm3=r.density_inv_fm3,passed=Bool(r.passed),
                failed_shells=r.failed_shells,provenance="retained_hash_bound",production_authorized=false))
        end
    end
    for (target,path) in ((ds,inputs[5]),(rs,inputs[7]))
        append!(target,[NamedTuple(r) for r in CSV.File(path) if r.sqrt_s_NN_GeV in energies])
    end
    profile=Models.load_freezeout_profile(profile="default")
    points=Models.build_freezeout_scan_points(energies;profile=profile,traversal=:sqrts_descending)
    retained_bgs=collect(CSV.File(inputs[1]))
    manifest=Dict{String,Any}("schema"=>"charged_gbu_dense_freezeout_v3","status"=>"running",
        "background"=>"FixedMuBConservedCharges quark-only BQS; rhoQ/rhoB=0.4; rhoS=0",
        "input_hashes"=>input_hashes,"settings"=>settings,"energies_GeV"=>energies,
        "production_authorized"=>false,"complete_curve_authorized"=>false,"threads"=>Threads.nthreads(),
        "reference_execution"=>"reuse_q0_spectrum; bitwise-parity unit tested", "new_background_solves"=>0)
    function checkpoint()
        for (name,data) in (("backgrounds",bgs),("direct_densities",direct),("reference_densities",reference),
                ("direct_shells",ds),("reference_shells",rs),("failures",failures),("timings",timings))
            isempty(data) || CSV.write(joinpath(output,name*".csv"),data)
        end
        CSV.write(joinpath(output,"ratio_comparison.csv"),D.ratio_table(direct,reference,energies))
        manifest["completed_density_evaluations"]=length(direct)+length(reference)
        manifest["failed_evaluations"]=length(failures)
        open(joinpath(output,"progress.json"),"w") do io
            JSON3.write(io,manifest)
        end
    end
    checkpoint()
    model=Models.create_model(:PNJL)
    seed=nothing
    for pt in points
        old=filter(r->r.T_MeV==pt.T_MeV && r.p_nodes==48,retained_bgs)
        if !isempty(old)
            push!(bgs,NamedTuple(only(old)));checkpoint();continue
        end
        started=time()
        solved=try
            D.V.background(model,pt.T_MeV,pt.muB_MeV,48;seed=seed)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,channel="all",route="background",reason=sprint(showerror,err)))
            checkpoint();continue
        end
        seed=solved.seed
        push!(bgs,solved.row)
        manifest["new_background_solves"]+=1
        println("[gbu-dense-parallel-bg] s=$(pt.sqrt_s_NN_GeV) residual=$(solved.bg.residual)")
        checkpoint();flush(stdout)
        bg=solved.bg
        tasks=map(R.CHANNELS) do channel
            Threads.@spawn begin
                results=NamedTuple[]
                for ref in (false,true)
                    start=time()
                    item=try
                        value=ref ? E.reference_density(bg,channel,settings) : R.density(bg,channel,settings)
                        (result=value,error="")
                    catch err
                        err isa InterruptException && rethrow()
                        (result=nothing,error=sprint(showerror,err))
                    end
                    push!(results,(reference=ref,channel=channel,elapsed_s=time()-start,item=item))
                end
                results
            end
        end
        for task in tasks,result in fetch(task)
            channel=result.channel
            route=result.reference ? "q0_lambda_reference" : "direct_finite_q"
            if result.item.result===nothing
                push!(failures,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,channel=String(channel),route=route,reason=result.item.error))
            else
                d=result.item.result
                target,shells=result.reference ? (reference,rs) : (direct,ds)
                push!(target,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,T_MeV=pt.T_MeV,muB_MeV=pt.muB_MeV,background_p=48,
                    channel=String(channel),route=route,density_inv_fm3=d.density,passed=d.passed,
                    failed_shells=count(r->!r.passed,d.rows),provenance="new_solve",production_authorized=false))
                append!(shells,[merge((sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,background_p=48),r) for r in d.rows])
                println("[gbu-dense-parallel] s=$(pt.sqrt_s_NN_GeV) $(channel) $(route) n=$(d.density) pass=$(d.passed)")
            end
            push!(timings,(sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,channel=String(channel),route=route,elapsed_s=result.elapsed_s))
            checkpoint();flush(stdout)
        end
        println("[gbu-dense-point] s=$(pt.sqrt_s_NN_GeV) elapsed_s=$(time()-started)")
    end
    all(R.hashfile(p)==h for (p,h) in input_hashes) || error("input changed during run")
    manifest["status"]="diagnostic_completed"
    manifest["all_conditional_gates"]=length(direct)==length(reference)==4length(energies) && isempty(failures) && all(r.passed for r in vcat(direct,reference))
    checkpoint()
    R.finish_output(output,hashes,manifest)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
