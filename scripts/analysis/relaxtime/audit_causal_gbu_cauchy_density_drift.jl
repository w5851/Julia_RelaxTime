"""Replay retained backgrounds after a numerical Cauchy fix; never re-solve BQS."""
module CausalGBUCauchyDensityDrift
include("causal_gbu_dense_execution.jl")
const E=CausalGBUDenseExecution
const R=E.R
using CSV,JSON3,LinearAlgebra

function main()
    R.method_contract()
    BLAS.set_num_threads(1)
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3")
    output=get(ENV,"GBU_DRIFT_OUTPUT",joinpath(base,"cauchy_stability_density_drift_20260905"))
    energies=parse.(Float64,split(get(ENV,"GBU_DRIFT_ENERGIES","200,7.7,3"),','))
    all(isfinite,energies) && all(>(0),energies) && length(unique(energies))==length(energies) || error("invalid energy grid")
    im=JSON3.read(read(joinpath(input,"manifest.json"),String))
    for (f,h) in pairs(im.output_hashes)
        R.hashfile(joinpath(input,String(f)))==String(h) || error("input hash mismatch: $(f)")
    end
    old=vcat(collect(CSV.File(joinpath(input,"direct_densities.csv"))),collect(CSV.File(joinpath(input,"reference_densities.csv"))))
    s=R.Settings(mesh=512,ne=128,nw=4800,nq=24,lower=1e-5,upper=56.,thermal=24.)
    all(getproperty(s,Symbol(k))==v for (k,v) in pairs(im.settings)) || error("retained settings mismatch")
    input_hashes=Dict(f=>R.hashfile(joinpath(input,f)) for f in ("backgrounds.csv","manifest.json","direct_densities.csv","reference_densities.csv"))
    hashes=R.start_output(output)
    rows,shells,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    function checkpoint()
        for (name,values) in (("density_drift",rows),("shells",shells),("failures",failures))
            isempty(values) || CSV.write(joinpath(output,name*".csv"),values)
        end
    end
    for energy in energies
        bg=R.saved_background(input,first(filter(r->r.sqrt_s_NN_GeV==energy,old)).T_MeV,48)
        tasks=map(R.CHANNELS) do channel
            Threads.@spawn begin
                records=NamedTuple[]
                for ref in (false,true)
                    started=time()
                    d=try
                        (value=ref ? E.reference_density(bg,channel,s) : R.density(bg,channel,s),reason="")
                    catch err
                        err isa InterruptException && rethrow()
                        (value=nothing,reason=sprint(showerror,err))
                    end
                    push!(records,(channel=channel,reference=ref,result=d,elapsed_s=time()-started))
                end
                records
            end
        end
        for task in tasks,r in fetch(task)
            route=r.reference ? "q0_lambda_reference" : "direct_finite_q"
            if r.result.value===nothing
                push!(failures,(sqrt_s_NN_GeV=energy,channel=String(r.channel),route=route,reason=r.result.reason))
            else
                d=r.result.value
                before=only(filter(x->x.sqrt_s_NN_GeV==energy && x.channel==String(r.channel) && x.route==route,old))
                change=abs(d.density-before.density_inv_fm3)/abs(before.density_inv_fm3)
                push!(rows,(sqrt_s_NN_GeV=energy,channel=String(r.channel),route=route,
                    old_density_inv_fm3=before.density_inv_fm3,density_inv_fm3=d.density,
                    relative_change=change,passed=d.passed,drift_passed=change<0.01,
                    failed_shells=count(x->!x.passed,d.rows),elapsed_s=r.elapsed_s,production_authorized=false))
                append!(shells,[merge((sqrt_s_NN_GeV=energy,),x) for x in d.rows])
                println("[cauchy-density-drift] s=$(energy) $(r.channel) $(route) change=$(change) passed=$(d.passed)")
            end
            checkpoint();flush(stdout)
        end
    end
    all(R.hashfile(joinpath(input,f))==h for (f,h) in input_hashes) || error("input changed")
    R.finish_output(output,hashes,Dict("status"=>"diagnostic_density_replay_completed","input_directory"=>input,
        "input_hashes"=>input_hashes,"energies_GeV"=>energies,"settings"=>s,"solver_called"=>false,
        "expected_density_rows"=>8length(energies),"failed_evaluations"=>length(failures),
        "all_conditional_gates"=>length(rows)==8length(energies) && isempty(failures) && all(r.passed && r.drift_passed for r in rows),
        "old_plot_preserved"=>true,"full_spectrum_certified"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
