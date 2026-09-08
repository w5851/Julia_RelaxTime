"""Compute the q=0 lambda-invariant reference on retained freezeout backgrounds.

Diagnostic only: this does not alter production defaults or historical outputs.
"""
module CausalGBUFreezeoutReference
include("causal_gbu_research_utils.jl")
using .CausalGBUResearch
using CSV, JSON3
const R = CausalGBUResearch

function main()
    base = joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input = get(ENV,"GBU_RESEARCH_INPUT",joinpath(base,"method_v1_freezeout_v3"))
    output = get(ENV,"GBU_RESEARCH_OUTPUT",joinpath(base,"method_v1_freezeout_q0_reference"))
    ispath(output) && error("refusing to overwrite $(output)")
    mkpath(output)
    input_files=("manifest.json","backgrounds.csv","sparse_densities.csv","sparse_shells.csv")
    input_hashes=Dict(f=>R.hashfile(joinpath(input,f)) for f in input_files)
    sparse=collect(CSV.File(joinpath(input,"sparse_densities.csv")))
    rows=NamedTuple[]
    shells=NamedTuple[]
    settings=R.Settings(mesh=512,ne=128,nw=4800,nq=24,lower=1e-5,upper=56.0,thermal=24.0)
    for energy in (200.0,7.7,3.0)
        bgrow=first(filter(r->r.sqrt_s_NN_GeV==energy && r.background_p==48,sparse))
        bg=R.saved_background(input,bgrow.T_MeV,48)
        for channel in R.CHANNELS
            d=R.density(bg,channel,settings;reference=true)
            append!(shells,[merge((sqrt_s_NN_GeV=energy,background_p=48),r) for r in d.rows])
            push!(rows,(sqrt_s_NN_GeV=energy,T_MeV=bg.T_MeV,muB_MeV=bg.muB_MeV,
                background_p=48,channel=String(channel),density_inv_fm3=d.density,
                passed=d.passed,failed_shells=count(r->!r.passed,d.rows),production_authorized=false))
            println("[gbu-freezeout-reference] s=$(energy) $(channel) n=$(d.density) pass=$(d.passed)")
        end
    end
    CSV.write(joinpath(output,"reference_densities.csv"),rows)
    CSV.write(joinpath(output,"reference_shells.csv"),shells)
    ratios=NamedTuple[]
    for group in values(Dict(e=>filter(r->r.sqrt_s_NN_GeV==e,rows) for e in (200.0,7.7,3.0)))
        m=Dict(r.channel=>r.density_inv_fm3 for r in group)
        e=first(group).sqrt_s_NN_GeV
        push!(ratios,(sqrt_s_NN_GeV=e,Kplus_over_pi_plus=m["K_plus"]/m["pi_plus"],
            Kminus_over_pi_minus=m["K_minus"]/m["pi_minus"],production_authorized=false))
    end
    CSV.write(joinpath(output,"reference_ratios.csv"),ratios)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,Dict(
            "stage"=>"freezeout_q0_reference","status"=>"diagnostic_stage_completed",
            "input_directory"=>abspath(input),"input_hashes"=>input_hashes,
            "settings"=>settings,"background_solves"=>0,"solver_called"=>false,
            "all_conditional_gates"=>length(rows)==12 && all(r.passed for r in rows),
            "production_authorized"=>false,"complete_curve_authorized"=>false))
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
