#!/usr/bin/env julia
# Reuse immutable Fortran channel outputs; compare layers before total eta/s.
using CSV, JSON, SHA
include(joinpath(@__DIR__,"..","..","relaxtime","run_gap_transport_scan.jl"))
include(joinpath(@__DIR__,"ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon: start_run, finish_run, write_json, write_rows

const FORTRAN_W_KEYS=(:ud_to_ud,:uu_to_uu,:us_to_us,:ss_to_ss,:udbar_to_udbar,
    :uubar_to_uubar,:uubar_to_ddbar,:usbar_to_usbar,:uubar_to_ssbar,:ssbar_to_uubar,:ssbar_to_ssbar)
const FORTRAN_WA_KEYS=(:ubardbar_to_ubardbar,:ubarubar_to_ubarubar,:ubarsbar_to_ubarsbar,
    :sbarsbar_to_sbarsbar,:dubar_to_dubar,:subar_to_subar)
const ALIGNMENT_SPECIES=(:u,:d,:s,:ubar,:dbar,:sbar)

function retained_fortran_rates(folder)
    values(name)=parse.(Float64,split(strip(read(joinpath(folder,"results",name),String))))[4:end]
    w,wa=values("w_ij.dat"),values("w_ija.dat")
    wn,wan=values("w_ij_n.dat"),values("w_ij_na.dat")
    length(w)==length(wn)==11 && length(wa)==length(wan)==6 || error("incomplete Fortran channel outputs")
    rates=NamedTuple{(FORTRAN_W_KEYS...,FORTRAN_WA_KEYS...)}((w...,wa...))
    nu,ns,nub,nsb=sqrt.((wn[2]/w[2],wn[4]/w[4],wan[2]/wa[2],wan[4]/wa[4]))
    densities=(u=nu,d=nu,s=ns,ubar=nub,dbar=nub,sbar=nsb)
    tau_inv=RelaxTime.RelaxationTime.relaxation_rates(densities,rates)
    tau=NamedTuple{ALIGNMENT_SPECIES}(Tuple(1/v for v in tau_inv))
    return (rates=rates,densities=densities,tau=tau)
end

function alignment_eta(qp,tp,tau,p,w)
    total=0.0
    for species in ALIGNMENT_SPECIES
        m=RT_ASR.get_mass(species,qp);mu=RT_ASR.get_mu(species,qp)
        for (node,weight) in zip(p,w)
            f=RT_ASR.distribution_with_anisotropy(species,node,m,mu,tp.T,tp.Φ,tp.Φbar,0.0,0.0)
            total+=weight*node^6/(node^2+m^2)*f*(1-f)*getproperty(tau,species)
        end
    end
    return 3total/(15pi^2*tp.T)
end

function fortran_alignment(inputdir,fortran_dir,outdir)
    original=JSON.parsefile(joinpath(fortran_dir,"input_manifest.json"))
    for file in original["source_files"]
        bytes2hex(sha256(read(joinpath(fortran_dir,"source_snapshot",file["path"]))))==file["sha256"] ||
            error("Fortran source snapshot changed: $(file["path"])")
    end
    targets=JSON.parsefile(joinpath(fortran_dir,"summary.json"))["comparison"]
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    sources=filter(p->endswith(p,".jl"),readdir(joinpath(PROJECT_ROOT,"scripts","relaxtime");join=true))
    start_run(outdir,inputdir,@__FILE__,(purpose="four xi=0 Fortran points, layered comparison",
        profile=Models.TransportSamplingProfiles.effective_config(profile),fortran_source=abspath(fortran_dir),
        fortran_manifest_sha256=bytes2hex(sha256(read(joinpath(fortran_dir,"input_manifest.json")))),
        external_digitized_tolerance=0.05,algebraic_replay_tolerance=1e-6,
        new_results="native profile; not a fit to digitized values or a claim that the Fortran grid converged"),
        extra_sources=sources)
    runtime=Models.TransportSamplingProfiles.build_runtime(profile)
    # Match the constants printed in the retained Fortran module, including its
    # rounded hbar*c. The production profile continues using its own constants.
    hc_fortran=197.33
    lambda_fortran=602.3/hc_fortran
    params=runtime.model.params
    values=NamedTuple{fieldnames(typeof(params))}(Tuple(getfield(params,k) for k in fieldnames(typeof(params))))
    matched=merge(values,(hbarc_MeV_fm=hc_fortran,alpha_em=1/137,
        Λ_inv_fm=lambda_fortran,m_ud0_inv_fm=5.5/hc_fortran,m_s0_inv_fm=140.7/hc_fortran,
        G_fm2=1.835/lambda_fortran^2,K_fm5=12.36/lambda_fortran^5,T0_inv_fm=210/hc_fortran,
        thermal_p_max_inv_fm=15.0))
    legacy=(p_num=128,t_num=16,model=Models.PNJLModel(typeof(params)(;matched...)),residual_norm_max=1e-8)
    write_json(joinpath(outdir,"model_parameters.json"),(production=values,matched_fortran=matched))
    current=(p_num=profile.background.p_nodes,t_num=profile.background.cos_nodes,
        model=runtime.model,residual_norm_max=1e-8)
    p15,w15=RelaxTime.GaussLegendre.gauleg(0.0,15.0,128)
    pinf,winf=RT_ASR._cached_semi_infinite_momentum_grid(128,5.0)
    backgrounds=NamedTuple[];channels=NamedTuple[];observables=NamedTuple[];failures=NamedTuple[]
    copied=NamedTuple[]
    for target in targets
        id=target["target_id"];folder=joinpath(fortran_dir,"points",id)
        println("FORTRAN ALIGNMENT START ",id);flush(stdout)
        try
            for name in ("point.csv","status.json","results/w_ij.dat","results/w_ija.dat","results/w_ij_n.dat","results/w_ij_na.dat")
                src=joinpath(folder,name);dst=joinpath(outdir,"inputs","fortran",id,name)
                mkpath(dirname(dst));cp(src,dst)
                push!(copied,(source=src,sha256=bytes2hex(sha256(read(src))),snapshot=relpath(dst,outdir)))
            end
            JSON.parsefile(joinpath(folder,"status.json"))["exitcode"]==0 || error("Fortran point failed")
            row=only(CSV.File(joinpath(folder,"point.csv")))
            old=retained_fortran_rates(folder)
            T,mu=row.T_MeV/ħc_MeV_fm,row.muB_MeV/(3ħc_MeV_fm)
            frozen_T,frozen_mu=row.T_MeV/hc_fortran,row.muB_MeV/(3hc_fortran)
            qp=(m=(u=row.m_u_MeV/hc_fortran,d=row.m_d_MeV/hc_fortran,s=row.m_s_MeV/hc_fortran),
                μ=(u=frozen_mu,d=frozen_mu,s=frozen_mu))
            tp=(T=frozen_T,Φ=row.Phi,Φbar=row.Phibar,ξ=0.0)
            tau=(u=row.tau_u_fm,d=row.tau_u_fm,s=row.tau_s_fm,
                ubar=row.tau_ubar_fm,dbar=row.tau_ubar_fm,sbar=row.tau_sbar_fm)
            eta_replay=alignment_eta(qp,tp,tau,p15,w15)
            eta_inf=alignment_eta(qp,tp,tau,pinf,winf)
            tau_replay=maximum(abs.(Tuple(old.tau)./Tuple(tau).-1))
            density_replay=maximum(ALIGNMENT_SPECIES) do s
                value=RT_ASR.number_density(s,RT_ASR.get_mass(s,qp),RT_ASR.get_mu(s,qp),frozen_T,tp.Φ,tp.Φbar,0.0;
                    p_grid=p15,p_w=w15,angle_nodes=1)
                abs(value/getproperty(old.densities,s)-1)
            end
            for (name,cfg) in (("matched_thermal15_p128",legacy),("production_thermal25_p96",current))
                hc=cfg.model.params.hbarc_MeV_fm
                temperature,chemical=row.T_MeV/hc,row.muB_MeV/(3hc)
                pair=GapTransportScanPhaseEquilibrium.solve_two_branch_candidates(temperature,chemical,0.0,cfg)
                eq=pair.stable
                entropy=Models.model_thermo(cfg.model,eq.x_state,eq.mu_vec,temperature;
                    xi=0.0,p_num=cfg.p_num,t_num=cfg.t_num)[3]
                push!(backgrounds,(target_id=id,configuration=name,
                    m_u_MeV=eq.masses[1]*hc,m_s_MeV=eq.masses[3]*hc,
                    Phi=eq.x_state[4],Phibar=eq.x_state[5],entropy=entropy,
                    m_u_relative=eq.masses[1]*hc/row.m_u_MeV-1,
                    m_s_relative=eq.masses[3]*hc/row.m_s_MeV-1,
                    Phi_relative=eq.x_state[4]/row.Phi-1,Phibar_relative=eq.x_state[5]/row.Phibar-1,
                    entropy_relative=entropy/row.entropy_fm3inv-1,distinct_branches=pair.distinct))
                write_rows(joinpath(outdir,"background.csv"),backgrounds)
                name=="production_thermal25_p96" || continue
                callback=d->begin
                    push!(channels,(target_id=id,channel=d.channel,fortran_rate=getproperty(old.rates,Symbol(d.channel)),
                        julia_rate=d.rate,relative_change=d.rate/getproperty(old.rates,Symbol(d.channel))-1,
                        estimated_error=d.estimated_error,seconds=d.seconds,nodes=d.nodes,evaluations=d.evaluations))
                    write_rows(joinpath(outdir,"channels.csv"),channels)
                end
                result=Models.TransportSamplingProfiles.solve_profile_transport(runtime,eq,T,mu;xi=0.0,on_channel=callback)
                eta_s=result.transport.eta/result.thermo_background.entropy
                push!(observables,(target_id=id,T_MeV=row.T_MeV,muB_MeV=row.muB_MeV,
                    fortran_eta_over_s=row.eta_over_s,digitized_eta_over_s=target["expected_eta_over_s"],
                    julia_eta_over_s=eta_s,relative_to_fortran=eta_s/row.eta_over_s-1,
                    relative_to_digitized=eta_s/target["expected_eta_over_s"]-1,
                    passes_external_5pct=abs(eta_s/target["expected_eta_over_s"]-1)<=0.05,
                    tau_u=result.tau.u,tau_ubar=result.tau.ubar,tau_s=result.tau.s,tau_sbar=result.tau.sbar,
                    tau_u_relative=result.tau.u/row.tau_u_fm-1,tau_ubar_relative=result.tau.ubar/row.tau_ubar_fm-1,
                    tau_s_relative=result.tau.s/row.tau_s_fm-1,tau_sbar_relative=result.tau.sbar/row.tau_sbar_fm-1,
                    frozen_eta_replay_relative=eta_replay/row.eta_fm3inv-1,
                    frozen_eta_tail_relative=eta_inf/eta_replay-1,frozen_tau_replay_max_relative=tau_replay,
                    frozen_density_replay_max_relative=density_replay))
                write_rows(joinpath(outdir,"observables.csv"),observables)
            end
            println("FORTRAN ALIGNMENT DONE ",id);flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(target_id=id,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FORTRAN ALIGNMENT FAILED ",id," ",sprint(showerror,err));flush(stdout)
        end
        write_json(joinpath(outdir,"fortran_inputs.json"),copied)
        GC.gc()
    end
    replay_passed=length(observables)==4 && all(r->max(abs(r.frozen_eta_replay_relative),
        r.frozen_tau_replay_max_relative,r.frozen_density_replay_max_relative)<=1e-6,observables)
    finish_run(outdir;points=length(targets),completed_points=length(observables),failures=length(failures),
        frozen_algebra_replay_passed=replay_passed,external_5pct_passed=count(r->r.passes_external_5pct,observables),
        comparison_complete=isempty(failures) && length(observables)==4)
    isempty(failures) && replay_passed || error("Fortran layer comparison incomplete or algebraic replay failed")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==3 || error("Usage: run_transport_fortran_alignment.jl INPUT_DIR FORTRAN_EVIDENCE_DIR NEW_OUTPUT_DIR")
    fortran_alignment(abspath(ARGS[1]),abspath(ARGS[2]),abspath(ARGS[3]))
end
