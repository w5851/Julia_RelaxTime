"""Analysis-only fixed-(T,muB) adapter to the existing charged GBU workflow.

No freezeout map is applied. The BQS solver, channel integration and all gates
are reused without changing their implementation or tolerances.
"""
module ChargedGBUBQSPoint
using CSV, JSON3, TOML, Distributed
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
isdefined(Main, :Models) || Base.include(Main, joinpath(ROOT, "src", "models", "Models.jl"))
isdefined(Main, :ChargedGBUResearchWorkflow) || Base.include(Main,
    joinpath(ROOT, "src", "models", "workflow_apps", "ChargedGBUResearchWorkflow.jl"))
const M = Main.Models
const W = Main.ChargedGBUResearchWorkflow

function parse_args(args)
    opts = Dict{Symbol,Any}(:T_MeV=>50.0, :muB_MeV=>620.0, :workers=>2,
        :output=>joinpath(ROOT,"data","outputs","results","relaxtime","meson_density",
            "charged_gbu_infinite","point_T50_muB620_bqs_20260921"))
    iseven(length(args)) || throw(ArgumentError("expected --option value pairs"))
    for i in 1:2:length(args)
        k,v=args[i:i+1]
        if k=="--T-MeV"; opts[:T_MeV]=parse(Float64,v)
        elseif k=="--muB-MeV"; opts[:muB_MeV]=parse(Float64,v)
        elseif k=="--workers"; opts[:workers]=parse(Int,v)
        elseif k=="--output"; opts[:output]=v
        else; throw(ArgumentError("unknown option: $k"))
        end
    end
    isfinite(opts[:T_MeV]) && opts[:T_MeV]>0 || throw(ArgumentError("T must be finite and positive"))
    isfinite(opts[:muB_MeV]) || throw(ArgumentError("muB must be finite"))
    opts[:workers] in 1:4 || throw(ArgumentError("workers must be 1..4"))
    return opts
end

function diagnostics(model,b,p)
    bg=b.bg; h=Main.Constants_PNJL.ħc_MeV_fm
    mu=M.conserved_mu_from_flavor(bg.mu...)
    rho=M.model_rho(model,b.seed[1:5],b.seed[6:8],bg.T;p_num=p,t_num=8,xi=0.)
    charges=M.conserved_densities_from_flavor(rho)
    return (evaluation_nodes=p,T_MeV=bg.T_MeV,muB_MeV=mu.mu_B*h,
        muQ_MeV=mu.mu_Q*h,muS_MeV=mu.mu_S*h,
        mu_u_MeV=bg.mu.u*h,mu_d_MeV=bg.mu.d*h,mu_s_MeV=bg.mu.s*h,
        mu_pi_plus_MeV=(bg.mu.u-bg.mu.d)*h,mu_K_plus_MeV=(bg.mu.u-bg.mu.s)*h,
        m_u_MeV=bg.m.u*h,m_d_MeV=bg.m.d*h,m_s_MeV=bg.m.s*h,
        Phi=bg.Phi,PhiBar=bg.PhiBar,rho_u_fm3=rho[1],rho_d_fm3=rho[2],rho_s_fm3=rho[3],
        rho_B_fm3=charges.rho_B,rho_Q_fm3=charges.rho_Q,rho_S_fm3=charges.rho_S,
        Q_over_B=charges.rho_Q/charges.rho_B,u_over_d=rho[1]/rho[2],
        charge_constraint_fm3=charges.rho_Q-.4charges.rho_B,
        residual_norm=bg.residual,omega_inv_fm4=b.omega)
end

function summarize(output,pt)
    results=Dict{String,Any}(); rows=NamedTuple[]
    for ch in W.R.CHANNELS
        r=W.readjson(joinpath(output,"$(ch).json"))
        density=r.density===nothing ? NaN : Float64(r.density)
        results[String(ch)]=(passed=Bool(r.passed),density=density)
        push!(rows,(channel=String(ch),density_fm3=density,passed=Bool(r.passed),
            status=String(r.status),reason=hasproperty(r,:reason) ? String(r.reason) : ""))
    end
    # NaN energy explicitly denotes an off-freezeout point, NOT a remapped 3 GeV point.
    r=W.ratio_row(NaN,pt,results)
    ratio=(T_MeV=r.T_MeV,muB_MeV=r.muB_MeV,Kplus_over_pi_plus=r.Kplus_over_pi_plus,
        Kminus_over_pi_minus=r.Kminus_over_pi_minus,plus_passed=r.plus_passed,minus_passed=r.minus_passed)
    CSV.write(joinpath(output,"densities.csv"),rows)
    CSV.write(joinpath(output,"ratios.csv"),[ratio])
    return ratio
end

function run_point(;T_MeV=50.,muB_MeV=620.,output,workers=2)
    ispath(output) && error("refusing to overwrite existing output: $output")
    mkpath(output)
    c=W.validate_config(TOML.parsefile(W.DEFAULT_CONFIG))
    pt=(T_MeV=Float64(T_MeV),muB_MeV=Float64(muB_MeV))
    hashes=W.snapshot(output,W.DEFAULT_CONFIG)
    W.writejson(joinpath(output,"run.json"),(point=pt,config=c,source_hashes=hashes,
        git_head=readchomp(`git -C $ROOT rev-parse HEAD`),julia_version=string(VERSION),
        classification="research",route="existing_charged_gbu_fixed_point_adapter",
        freezeout_mapping_applied=false,meson_feedback=false,xi=0.0))
    model=M.create_model(:PNJL)
    # The saved 3 GeV background is a continuation SEED only, not the new solution.
    reference=joinpath(ROOT,"data","outputs","results","relaxtime","meson_density",
        "charged_gbu_infinite","freezeout_20260907_v2","energy_3.0","background.json")
    old=JSON3.read(read(reference,String))
    b=W.background(model,pt,c;seed=Float64.(old.seed))
    W.writejson(joinpath(output,"background.json"),b)
    cold=W.background(model,pt,c)
    fine_config=deepcopy(c);fine_config["numerics"]["background_nodes"]=96
    fine=W.background(model,pt,fine_config;seed=b.seed)
    W.writejson(joinpath(output,"background_p96.json"),fine)
    W.writejson(joinpath(output,"background_cold.json"),cold)
    rows=[merge((solve_nodes=48,),diagnostics(model,b,p)) for p in (48,96,192)]
    push!(rows,merge((solve_nodes=96,),diagnostics(model,fine,192)))
    CSV.write(joinpath(output,"bqs_diagnostics.csv"),rows)
    seed_difference=maximum(abs,b.seed-cold.seed)
    node_difference=maximum(abs,b.seed-fine.seed)
    bg_ok=seed_difference<1e-6 && node_difference<1e-6 && all(r->
        abs(r.muB_MeV-pt.muB_MeV)<1e-5 && abs(r.Q_over_B-.4)<1e-4 &&
        abs(r.rho_S_fm3)<1e-7 && r.residual_norm<=c["gates"]["background_residual"],rows)
    W.writejson(joinpath(output,"background_checks.json"),(passed=bg_ok,
        maximum_seed_difference=seed_difference,maximum_p48_p96_solution_difference=node_difference,
        scope="cold_vs_continuation_and_background_quadrature_not_global_phase_proof"))
    println("[point-background] ",first(rows)," checks_passed=",bg_ok);flush(stdout)
    jobs=[(bg=b.bg,ch=ch,dir=output,energy="fixed_T$(T_MeV)_muB$(muB_MeV)") for ch in W.R.CHANNELS]
    if workers==1
        foreach(job->W.run_job(job,c),jobs)
    else
        ids=addprocs(workers;exeflags=`--project=$ROOT --threads=1`,dir=ROOT)
        try
            for id in ids
                remotecall_wait(Base.include,id,Main,@__FILE__)
            end
            pmap(job->Main.ChargedGBUResearchWorkflow.run_job(job,c),WorkerPool(ids),jobs)
        finally
            rmprocs(ids)
        end
    end
    ratio=summarize(output,pt)
    all(W.hashfile(joinpath(ROOT,p))==h for (p,h) in hashes) || error("source changed during run")
    passed=bg_ok && ratio.plus_passed && ratio.minus_passed
    W.writejson(joinpath(output,"manifest.json"),(status=passed ? "accepted_fixed_point" : "fixed_point_with_failed_checks",
        point=pt,ratio=ratio,background_checks_passed=bg_ok,source_hashes=hashes,
        observable="fixed_quark_only_gbu_partial_yield",meson_feedback=false,
        charge_to_baryon_ratio=.4,strangeness_density_fm3=0.,thermal_target="infinity",
        freezeout_mapping_applied=false,formal_baseline=false,output_hashes=W.output_hashes(output)))
    println("[point-result] ",ratio," all_checks_passed=",passed)
    return passed
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && (run_point(;parse_args(ARGS)...) || exit(1))
end
