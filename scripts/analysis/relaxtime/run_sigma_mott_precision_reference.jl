#!/usr/bin/env julia
# quote-remote only. Fixed-background reference with bounded inter-s polarization cache.
include(joinpath(@__DIR__, isfile(joinpath(@__DIR__,"run_sigma_mott_spectrum_v2.jl")) ?
    "run_sigma_mott_spectrum_v2.jl" : "run_sigma_mott_spectrum.jl"))

function sigma_reference_rate(row, process, state, runtime, folder)
    qp,tp,kc=state.qp,state.tp,state.kc
    options=runtime.tau_kwargs
    problem=AS._prepare_thermal_cache_problem(process,qp,tp,kc;options...)
    cfg=runtime.profile.sampling
    measure=SIG_TS.thermal_measure(problem.plan,problem.context.bounds;relative_velocity=:cm)
    s0=first(problem.context.bounds)
    feature_state=(quark_params=qp,thermo_params=tp,K_coeffs=kc)
    features=SIG_TS.scattering_feature_seeds(process,feature_state,measure,s0)
    thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
    SIG_PC.reset_cache!()
    calls=0
    function direct(s)
        # Keep all reuse inside each cross-section evaluation, but never allow
        # the cross-s quantized polarization dictionary to grow with the table.
        SIG_PC.reset_cache!()
        calls+=1
        calls%8192==0 && GC.gc()
        SIG_CS.total_cross_section(process,s,qp,tp,kc;n_points=runtime.profile.tau_integral.t_nodes)
    end
    started=time()
    built=SIG_TS.build_thermal_weighted_table(direct,measure,s0;rtol=cfg.rtol,atol=cfg.atol,
        weighted_window=cfg.weighted_window,seed_order=cfg.seed_order,thresholds=thresholds,
        extra_seeds=features.seeds,max_evaluations=cfg.max_evaluations,max_nodes=cfg.max_nodes,
        max_rounds=cfg.max_rounds)
    CSV.write(joinpath(folder,String(row.point_id)*"__"*String(process)*"__history.csv"),built.history)
    table=SIG_TS.validated_table(built)
    SIG_PC.reset_cache!();GC.gc()
    (channel=String(process),seconds=time()-started,rate=built.rate,estimated_error=built.estimated_error,
     evaluations=built.evaluations,nodes=length(table.nodes),events=measure.events,
     unique_events=length(measure.nodes),guard_phase=built.guard_phase,feature_count=length(features.features),reason=built.reason)
end

function run_sigma_precision_reference()
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,ARGS[6]*".toml"))
    runtime=SP_OTHER.build_runtime(profile)
    wanted=get(ENV,"MOTT_POINT_IDS","all")
    backgrounds=collect(CSV.File(joinpath(dirname(OUTPUT_DIR),"spectrum_v2","background.csv")))
    wanted=="all" || filter!(r->String(r.point_id) in split(wanted,','),backgrounds)
    isempty(backgrounds) && error("No precision-reference points")
    previous=collect(CSV.File(joinpath(INPUT_DIR,"previous","observables.csv")))
    oldweights=collect(CSV.File(joinpath(INPUT_DIR,"previous","species_weights.csv")))
    for name in ("left_s6_v1","right_s6_v1")
        append!(previous,collect(CSV.File(joinpath(dirname(OUTPUT_DIR),name,"observables.csv"))))
        append!(oldweights,collect(CSV.File(joinpath(dirname(OUTPUT_DIR),name,"species_weights.csv"))))
    end
    observables,rates,weights,contributions=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    selected_channels=split(get(ENV,"SIGMA_CHANNELS","all"),',')
    process_list=selected_channels==["all"] ? RTM.REQUIRED_RATE_KEYS_FOR_TAU : Symbol.(selected_channels)
    all(p->p in RTM.REQUIRED_RATE_KEYS_FOR_TAU,process_list) || error("Invalid selected precision channel")
    for row in backgrounds
        began=time()
        original=only(r for r in previous if abs(Float64(r.xi)-Float64(row.xi))<1e-12)
        all(abs(Float64(getproperty(original,k))-Float64(getproperty(row,k)))<1e-10
            for k in (:m_u,:m_d,:m_s,:Phi,:Phibar,:xi)) || error("Reference background differs")
        # Reuse the exact observable-background bytes, not a second solve.
        state=frozen_state(original,profile,32,16)
        info=(point_id=String(row.point_id),curve_id="other07",window=String(row.window),xi=Float64(row.xi))
        for process in process_list
            println("REFERENCE_CHANNEL_START ",row.point_id," ",process);flush(stdout)
            r=sigma_reference_rate(row,process,state,runtime,OUTPUT_DIR)
            push!(rates,merge(info,r));CSV.write(joinpath(OUTPUT_DIR,"rates.csv"),rates)
            println("REFERENCE_CHANNEL_DONE ",row.point_id," ",process," seconds=",r.seconds);flush(stdout)
        end
        length(process_list)==17 || continue
        rr=filter(r->r.point_id==info.point_id,rates)
        rate_dict=Dict(Symbol(r.channel)=>r.rate for r in rr)
        density=NamedTuple{SPECIES_OTHER}(Tuple(Float64(getproperty(original,Symbol("n_",sp))) for sp in SPECIES_OTHER))
        inverse=RTM.relaxation_rates(density,rate_dict)
        tau=NamedTuple{SPECIES_OTHER}(Tuple(1/getproperty(inverse,sp) for sp in SPECIES_OTHER))
        for old in filter(r->String(r.point_id)==String(original.point_id),oldweights)
            w=NamedTuple(old)
            push!(weights,merge(w,info,(tau_fm=getproperty(tau,Symbol(old.species)),)))
        end
        terminal=map((:eta,:zeta,:sigma)) do key
            sum(Float64(getproperty(w,Symbol(key,"_coefficient")))*getproperty(tau,Symbol(w.species))
                for w in weights if w.point_id==info.point_id)
        end
        observed=merge(NamedTuple(original),info,(variant=ARGS[6],point_kind="fixed_background_fresh_s_reference",
            delta_xi=Float64(row.delta_xi),eta=terminal[1],zeta=terminal[2],sigma=terminal[3],
            eta_over_s=terminal[1]/Float64(original.entropy),zeta_over_s=terminal[2]/Float64(original.entropy),
            sigma_over_T=terminal[3]/Float64(original.T_fm),seconds=time()-began),
            NamedTuple{Tuple(Symbol("tau_",sp) for sp in SPECIES_OTHER)}(Tuple(tau)),
            NamedTuple{Tuple(Symbol("tauinv_",sp) for sp in SPECIES_OTHER)}(Tuple(inverse)))
        push!(observables,observed)
        for r in rr
            single=Dict(Symbol(v.channel)=>(v.channel==r.channel ? Float64(r.rate) : 0.0) for v in rr)
            partial=RTM.relaxation_rates(density,single)
            for sp in SPECIES_OTHER
                rc=getproperty(partial,sp)
                push!(contributions,merge(info,(species=String(sp),channel=String(r.channel),rate_fm2=r.rate,
                    collision_density_weight_inv_fm3=rc/r.rate,inverse_tau_contribution_inv_fm=rc,
                    inverse_tau_inv_fm=getproperty(inverse,sp))))
            end
        end
        for (name,rows) in (("observables",observables),("species_weights",weights),("channel_contributions",contributions))
            CSV.write(joinpath(OUTPUT_DIR,name*".csv"),rows)
        end
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",points=length(backgrounds),
        channels=length(rates),complete_transport=length(process_list)==17,
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",production_kernel_changed=false,
        effective_profile=SP_OTHER.effective_config(profile),
        interpretation="All direct s samples use unchanged production cross sections after resetting the quantized polarization cache; within-s reuse is retained. Background, densities and terminal coefficients are frozen from completed s6 rows.",
        script_sha256=digest(@__FILE__),library_sha256=SIG_LIBRARY_SHA,
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_sigma_precision_reference()
end
