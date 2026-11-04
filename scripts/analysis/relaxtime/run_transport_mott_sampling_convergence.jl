#!/usr/bin/env julia
# Fixed-background s-resolution checks for the two muB=900 conductivity bumps.
# Run only on quote-remote. Arguments: SOURCE_ROOT INPUT_DIR OUTPUT_DIR SELECTION [POLE_MESON]
# SELECTION is pilot, remaining, all, or comma-separated point__channel keys.
using CSV, JSON3, SHA, Dates, LinearAlgebra

length(ARGS) in (4,5) || error("expected SOURCE_ROOT INPUT_DIR OUTPUT_DIR SELECTION [POLE_MESON]")
const ROOT = abspath(ARGS[1])
const INPUT = abspath(ARGS[2])
const OUTPUT = abspath(ARGS[3])
isdir(OUTPUT) && !isempty(readdir(OUTPUT)) && error("refusing to overwrite a retained run")
mkpath(OUTPUT)
include(joinpath(ROOT,"src","models","Models.jl"))
const AS = Main.RelaxTime.AverageScatteringRate
const TS = AS.ThermalWeightedSampling
const CS = Main.RelaxTime.TotalCrossSection
const PC = Main.RelaxTime.PolarizationCache
const EC = Main.RelaxTime.EffectiveCouplings
const RTC = Main.RelaxTime.TransportCoefficients
const RT = Main.RelaxTime.RelaxationTime
const HC = Main.Constants_PNJL.ħc_MeV_fm
const POINTS = ("p004","p005","p006","p014","p015","p016")
const CHANNELS = (:udbar_to_udbar,:dubar_to_dubar,:uubar_to_uubar,:uubar_to_ddbar)
const PILOT = "p005__dubar_to_dubar"
const POLE_MESON = length(ARGS)==5 ? ARGS[5] : "pi"
const REFERENCE_TARGET = 2e-5
BLAS.set_num_threads(1)

digest(path) = bytes2hex(sha256(read(path)))
function write_json(path,data)
    open(path,"w") do io
        JSON3.pretty(io,data)
        println(io)
    end
end
function write_rows(path,rows)
    isempty(rows) || CSV.write(path,rows)
end
function frozen_state(row,profile)
    qp0=(m=(u=Float64(row.m_u),d=Float64(row.m_d),s=Float64(row.m_s)),
        μ=(u=Float64(row.muq_fm),d=Float64(row.muq_fm),s=Float64(row.muq_fm)))
    tp=(T=Float64(row.T_fm),Φ=Float64(row.Phi),Φbar=Float64(row.Phibar),ξ=Float64(row.xi))
    abs(tp.T*HC-row.T_MeV)<1e-8 || error("temperature units differ")
    a=profile.a_integral
    p,w=AS._cached_semi_infinite_momentum_grid(a.p_nodes,a.scale_inv_fm)
    av=Main.RelaxTime.AFieldBuilder.build_A_triplet(qp0,tp,p,w;cos_nodes=a.cos_nodes,use_aniso=true)
    qp=merge(qp0,(A=av,))
    gu,gs=EC.calculate_G_from_A(av.u,qp.m.u),EC.calculate_G_from_A(av.s,qp.m.s)
    kc=EC.calculate_effective_couplings(Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    (quark_params=qp,thermo_params=tp,K_coeffs=kc)
end

function rate_options(profile;rtol=1e-3)
    q=profile.tau_integral
    p,w=AS._cached_semi_infinite_momentum_grid(q.p_nodes,q.scale_inv_fm)
    (p_nodes=q.p_nodes,angle_nodes=q.cos_nodes,phi_nodes=q.phi_nodes,p_grid=p,p_w=w,
        n_sigma_points=q.t_nodes,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
        sampling_config=AS.ThermalSamplingConfig(rtol=rtol))
end

# Reference meshes do not use the adaptive cell-refinement rule. They combine
# logarithmic x=s-s0 and uniform sqrt(x) grids; interpolation is in raw sigma,
# unlike the production weighted x*sigma interpolation near threshold.
function reference_grid(measure,s0,order,features,thresholds;phase=0.0)
    a,b=first(measure.nodes),last(measure.nodes)
    xa,xb=a-s0,b-s0
    fractions=phase==0 ? collect(range(0.0,1.0;length=order+1)) :
        vcat(0.0,[(k+phase)/order for k in 0:order-1],1.0)
    fine_fractions=phase==0 ? collect(range(0.0,1.0;length=2order+1)) :
        vcat(0.0,[(k+phase)/(2order) for k in 0:2order-1],1.0)
    grid=vcat([a,b],s0 .+ exp.(log(xa) .+ fractions .* log(xb/xa)),
        s0 .+ (sqrt(xa) .+ fine_fractions .* (sqrt(xb)-sqrt(xa))).^2,
        features)
    for threshold in thresholds
        a<threshold<b || continue
        push!(grid,threshold)
        radius=min(threshold-a,b-threshold,1.0)
        for offset in exp.(range(log(max(128eps(threshold),radius*1e-8)),log(radius);length=24))
            a<threshold-offset<b && push!(grid,threshold-offset)
            a<threshold+offset<b && push!(grid,threshold+offset)
        end
    end
    sort!(unique!(filter(s->a<=s<=b,grid)))
end

function compare_tables(measure,left,right)
    lr=TS.integrate_measure(measure,left)
    rr=TS.integrate_measure(measure,right)
    l1=TS.integrate_measure(measure,s->abs(left(s)-right(s)))
    (left=lr,right=rr,l1=l1,relative_l1=l1/rr,relative_signed=lr/rr-1)
end

function save_table(path,table)
    write_rows(path,[(s_inv_fm2=s,x_inv_fm2=s-table.s0,sigma_fm2=v,
        x_sigma=(s-table.s0)*v) for (s,v) in zip(table.nodes,table.sigma)])
end

function independent_reference(process,state,measure,s0,features,thresholds,base,fine,jobdir,t_nodes)
    values=Dict{Float64,Float64}()
    direct=s->get!(values,s) do
        v=CS.total_cross_section(process,s,state.quark_params,state.thermo_params,state.K_coeffs;n_points=t_nodes)
        isfinite(v) && v>=0 || throw(DomainError((s,v),"invalid direct reference cross section"))
        v
    end
    previous=nothing
    history=NamedTuple[]
    for order in (192,384,768,1536,3072,6144)
        started=time()
        PC.reset_cache!()
        grid=reference_grid(measure,s0,order,features,thresholds)
        table=TS.WeightedSigmaTable(grid,direct.(grid),s0,s0)
        rate=TS.integrate_measure(measure,table)
        delta=previous===nothing ? Inf : compare_tables(measure,previous,table).relative_l1
        push!(history,(order=order,nodes=length(grid),evaluations=length(values),rate=rate,
            previous_l1_relative=delta,seconds=time()-started))
        write_rows(joinpath(jobdir,"reference_history.csv"),history)
        println("  reference order=",order," L1=",delta," seconds=",round(time()-started;digits=2));flush(stdout)
        if previous!==nothing && delta<=REFERENCE_TARGET
            shifted_grid=reference_grid(measure,s0,order,features,thresholds;phase=sqrt(2)-1)
            shifted=TS.WeightedSigmaTable(shifted_grid,direct.(shifted_grid),s0,s0)
            shifted_check=compare_tables(measure,shifted,table)
            if shifted_check.relative_l1<=REFERENCE_TARGET
                save_table(joinpath(jobdir,"reference_table.csv"),table)
                save_table(joinpath(jobdir,"reference_shifted_table.csv"),shifted)
                base_check=compare_tables(measure,base,table)
                fine_check=compare_tables(measure,fine,table)
                return (table=table,rate=rate,order=order,evaluations=length(values),
                    previous_l1=delta,shifted_l1=shifted_check.relative_l1,
                    base_l1=base_check.relative_l1,base_signed=base_check.relative_signed,
                    fine_l1=fine_check.relative_l1,fine_signed=fine_check.relative_signed,
                    direct=direct)
            end
        end
        previous=table
    end
    error("independent s-reference did not converge within its retained budget")
end

function local_audit(measure,s0,pole,base,fine,reference,direct,jobdir)
    count=length(measure.nodes)
    chosen=Set(1:min(count,2048))
    mass=Float64(pole.mass_MeV)/HC
    if mass^2>s0
        center=searchsortedfirst(measure.nodes,mass^2)
        union!(chosen,max(1,center-1024):min(count,center+1023))
    end
    rows=NamedTuple[]
    exact_rate=0.0;base_error=0.0;fine_error=0.0;reference_error=0.0
    for i in sort!(collect(chosen))
        s,w=measure.nodes[i],measure.weights[i]
        actual=direct(s)
        bv,fv,rv=base(s),fine(s),reference(s)
        exact_rate+=w*actual
        base_error+=w*abs(bv-actual);fine_error+=w*abs(fv-actual)
        reference_error+=w*abs(rv-actual)
        push!(rows,(event_index=i,s_inv_fm2=s,x_inv_fm2=s-s0,thermal_weight=w,
            direct_sigma_fm2=actual,base_sigma_fm2=bv,fine_sigma_fm2=fv,reference_sigma_fm2=rv))
    end
    write_rows(joinpath(jobdir,"direct_event_audit.csv"),rows)
    (events_checked=length(rows),exact_rate_subset=exact_rate,
        base_l1_subset=base_error,fine_l1_subset=fine_error,reference_l1_subset=reference_error)
end

function band_contributions(measure,s0,base,fine,reference,jobdir)
    edges=[0.0,1e-6,1e-5,1e-4,1e-3,1e-2,0.05,0.2,0.6,1.5,Inf]
    count=zeros(Int,length(edges)-1)
    a=zeros(length(count));b=similar(a);c=similar(a);err=similar(a)
    fill!(b,0);fill!(c,0);fill!(err,0)
    for (s,w) in zip(measure.nodes,measure.weights)
        band=clamp(searchsortedlast(edges,s-s0),1,length(count))
        bv,fv,rv=base(s),fine(s),reference(s)
        count[band]+=1;a[band]+=w*bv;b[band]+=w*fv;c[band]+=w*rv;err[band]+=w*abs(bv-rv)
    end
    rows=[(x_left_inv_fm2=edges[i],x_right_inv_fm2=edges[i+1],events=count[i],
        base_rate_fm2=a[i],fine_rate_fm2=b[i],reference_rate_fm2=c[i],base_reference_l1_fm2=err[i])
        for i in eachindex(count)]
    write_rows(joinpath(jobdir,"s_band_contributions.csv"),rows)
end

function run_job(row,process,old_rate,pole,profile)
    job=String(row.point_id)*"__"*String(process)
    jobdir=joinpath(OUTPUT,job);mkpath(jobdir)
    state=frozen_state(row,profile)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    options=rate_options(profile)
    PC.reset_cache!();GC.gc()
    t0=time()
    problem=AS._prepare_thermal_cache_problem(process,qp,tp,kc;options...)
    base=AS._build_thermal_weighted_cache(process,problem)
    reproduction=base.rate/Float64(old_rate.rate)-1
    abs(reproduction)<=1e-6 || error("base production replay differs at $job: $reproduction")
    save_table(joinpath(jobdir,"base_table.csv"),base.table)
    println("BASE ",job," replay=",reproduction," seconds=",round(time()-t0;digits=2));flush(stdout)
    measure=TS.thermal_measure(problem.plan,problem.context.bounds;relative_velocity=:cm)
    s0=first(problem.context.bounds)
    features=TS.scattering_feature_seeds(process,state,measure,s0)
    thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
    histories=NamedTuple[(configuration="base",rtol=1e-3,rate=base.rate,
        estimated_error=base.estimated_error,nodes=length(base.table.nodes),evaluations=base.statistics.evaluations)]
    fine_table=base.table
    for rtol in (1e-4,1e-5)
        PC.reset_cache!()
        direct=s->CS.total_cross_section(process,s,qp,tp,kc;n_points=profile.tau_integral.t_nodes)
        built=TS.build_thermal_weighted_table(direct,measure,s0;rtol=rtol,atol=1e-12,
            weighted_window=0.6,seed_order=48,thresholds=thresholds,extra_seeds=features.seeds,
            max_evaluations=160_000,max_nodes=40_000,max_rounds=80)
        table=TS.validated_table(built)
        tag=rtol==1e-4 ? "tight1" : "tight2"
        save_table(joinpath(jobdir,tag*"_table.csv"),table)
        write_rows(joinpath(jobdir,tag*"_history.csv"),built.history)
        push!(histories,(configuration=tag,rtol=rtol,rate=built.rate,
            estimated_error=built.estimated_error,nodes=length(table.nodes),evaluations=built.evaluations))
        fine_table=table
        write_rows(joinpath(jobdir,"adaptive_rates.csv"),histories)
        println("TIGHT ",job," rtol=",rtol," rate=",built.rate);flush(stdout)
    end
    PC.reset_cache!();GC.gc()
    ref=independent_reference(process,state,measure,s0,features.seeds,thresholds,
        base.table,fine_table,jobdir,profile.tau_integral.t_nodes)
    PC.reset_cache!()
    audit=local_audit(measure,s0,pole,base.table,fine_table,ref.table,ref.direct,jobdir)
    band_contributions(measure,s0,base.table,fine_table,ref.table,jobdir)
    stats=(job_id=job,point_id=String(row.point_id),window=String(row.window),channel=String(process),
        T_MeV=Float64(row.T_MeV),muB_MeV=Float64(row.muB_MeV),xi=Float64(row.xi),
        s0_inv_fm2=s0,thermal_events=measure.events,unique_events=length(measure.nodes),audit_meson=String(pole.meson),
        original_rate=Float64(old_rate.rate),base_rate=base.rate,reproduction_relative=reproduction,
        tight1_rate=histories[2].rate,tight2_rate=histories[3].rate,reference_rate=ref.rate,
        reference_order=ref.order,reference_previous_l1=ref.previous_l1,reference_shifted_l1=ref.shifted_l1,
        base_reference_l1=ref.base_l1,base_reference_signed=ref.base_signed,
        tight2_reference_l1=ref.fine_l1,tight2_reference_signed=ref.fine_signed,
        direct_events_checked=audit.events_checked,direct_subset_rate_fraction=audit.exact_rate_subset/ref.rate,
        direct_base_l1_over_total=audit.base_l1_subset/ref.rate,
        direct_tight2_l1_over_total=audit.fine_l1_subset/ref.rate,
        direct_reference_l1_over_total=audit.reference_l1_subset/ref.rate,
        seconds=time()-t0)
    files=Dict(name=>digest(joinpath(jobdir,name)) for name in readdir(jobdir))
    write_json(joinpath(jobdir,"completed.json"),(status="completed",metrics=stats,files=files))
    stats
end

function main()
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    rows=collect(CSV.File(joinpath(INPUT,"source_rows.csv")))
    rates=collect(CSV.File(joinpath(INPUT,"source_rates.csv")))
    poles=collect(CSV.File(joinpath(INPUT,"poles.csv")))
    selection=ARGS[4]
    legacy_selection=selection in ("all","pilot","remaining")
    jobs=[(row,process) for row in rows if !legacy_selection || String(row.point_id) in POINTS for process in CHANNELS]
    jobs=filter(jobs) do (r,p)
        key=String(r.point_id)*"__"*String(p)
        selection=="all" || (selection=="pilot" && key==PILOT) ||
            (selection=="remaining" && key!=PILOT) || key in split(selection,',')
    end
    isempty(jobs) && error("empty study selection")
    if !legacy_selection
        requested=Set(split(selection,','))
        actual=Set(String(r.point_id)*"__"*String(p) for (r,p) in jobs)
        actual==requested || error("one or more requested point/channel keys are absent or unsupported")
    end
    write_json(joinpath(OUTPUT,"run_config.json"),(purpose="s-resolution at fixed formal backgrounds and outer measure",
        selection=selection,jobs=[String(r.point_id)*"__"*String(p) for (r,p) in jobs],
        source_root=ROOT,source_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        script_sha256=digest(@__FILE__),julia_version=string(VERSION),threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(),profile=Models.TransportSamplingProfiles.effective_config(profile),
        adaptive_rtols=(1e-3,1e-4,1e-5),reference_target=REFERENCE_TARGET,audit_meson=POLE_MESON,
        inputs=Dict(name=>digest(joinpath(INPUT,name)) for name in ("source_rows.csv","source_rates.csv","poles.csv"))))
    results=NamedTuple[]
    for (row,process) in jobs
        old=only(r for r in rates if String(r.point_id)==String(row.point_id) && Symbol(r.channel)==process)
        pole=only(r for r in poles if String(r.point_id)==String(row.point_id) && r.meson==POLE_MESON && r.direction=="forward")
        println("START ",row.point_id," ",process);flush(stdout)
        result=run_job(row,process,old,pole,profile)
        push!(results,result)
        write_rows(joinpath(OUTPUT,"channels.csv"),results)
        write_json(joinpath(OUTPUT,"progress.json"),(completed=length(results),expected=length(jobs),last_job=result.job_id))
        println("DONE ",result.job_id," base_L1=",result.base_reference_l1," seconds=",round(result.seconds;digits=2));flush(stdout)
        PC.reset_cache!();GC.gc()
    end
    write_json(joinpath(OUTPUT,"summary.json"),(status="completed_analysis",jobs=length(results),
        finished_utc=string(now(UTC)),channels_sha256=digest(joinpath(OUTPUT,"channels.csv")),
        max_base_reference_l1=maximum(r.base_reference_l1 for r in results),
        max_base_reference_signed=maximum(abs(r.base_reference_signed) for r in results),
        max_tight2_reference_l1=maximum(r.tight2_reference_l1 for r in results)))
end

main()
