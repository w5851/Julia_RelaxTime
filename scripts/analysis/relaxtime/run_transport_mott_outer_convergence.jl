#!/usr/bin/env julia
# Stream actual thermal events on refined outer grids; no event histogram.
# quote-remote only: SOURCE_ROOT INPUT_DIR OUTPUT_DIR POINT_IDS [CHANNELS S_RUN_DIRS]
using CSV, JSON3, SHA, Dates, LinearAlgebra
length(ARGS) in (4,6) || error("expected SOURCE_ROOT INPUT_DIR OUTPUT_DIR POINT_IDS [CHANNELS S_RUN_DIRS]")
const ROOT=abspath(ARGS[1])
const INPUT=abspath(ARGS[2])
const OUTPUT=abspath(ARGS[3])
const CASE=dirname(OUTPUT)
isdir(OUTPUT) && !isempty(readdir(OUTPUT)) && error("refusing to overwrite retained analysis")
mkpath(OUTPUT)
include(joinpath(ROOT,"src","models","Models.jl"))
const AS=Main.RelaxTime.AverageScatteringRate
const TS=AS.ThermalWeightedSampling
const CS=Main.RelaxTime.TotalCrossSection
const PC=Main.RelaxTime.PolarizationCache
const EC=Main.RelaxTime.EffectiveCouplings
const CHANNELS=length(ARGS)==6 ? Tuple(Symbol.(split(ARGS[5],','))) : (:udbar_to_udbar,:dubar_to_dubar)
const S_RUN_DIRS=length(ARGS)==6 ? Tuple(split(ARGS[6],',')) : ("pilot_v1","remaining_v1")
length(unique(CHANNELS))==length(CHANNELS)>0 || error("channels must be nonempty and unique")
all(c->c in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU,CHANNELS) || error("unsupported channel")
all(r->occursin(r"^[a-zA-Z0-9_-]+$",r),S_RUN_DIRS) || error("invalid retained s-run names")
const CONFIGS=((name="base",p=48,c=36,phi=36),
    (name="joint1",p=72,c=48,phi=54),(name="joint2",p=96,c=60,phi=72))
const S_TARGET=2e-5
const OUTER_TARGET=1e-3
BLAS.set_num_threads(1)
digest(p)=bytes2hex(sha256(read(p)))
function write_json(path,data)
    open(path,"w") do io
        JSON3.pretty(io,data);println(io)
    end
end

function frozen_state(row,profile)
    qp0=(m=(u=Float64(row.m_u),d=Float64(row.m_d),s=Float64(row.m_s)),
        μ=(u=Float64(row.muq_fm),d=Float64(row.muq_fm),s=Float64(row.muq_fm)))
    tp=(T=Float64(row.T_fm),Φ=Float64(row.Phi),Φbar=Float64(row.Phibar),ξ=Float64(row.xi))
    p,w=AS._cached_semi_infinite_momentum_grid(profile.a_integral.p_nodes,profile.a_integral.scale_inv_fm)
    av=Main.RelaxTime.AFieldBuilder.build_A_triplet(qp0,tp,p,w;
        cos_nodes=profile.a_integral.cos_nodes,use_aniso=true)
    qp=merge(qp0,(A=av,))
    gu,gs=EC.calculate_G_from_A(av.u,qp.m.u),EC.calculate_G_from_A(av.s,qp.m.s)
    kc=EC.calculate_effective_couplings(Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    (quark_params=qp,thermo_params=tp,K_coeffs=kc)
end

mutable struct Support
    lower::Float64
    upper::Float64
end
@inline visit!(a::Support,s,w)=(a.lower=min(a.lower,s);a.upper=max(a.upper,s);nothing)
mutable struct Accumulator{F,N}
    query::F
    totals::NTuple{N,Float64}
    corrections::NTuple{N,Float64}
end
@inline function visit!(a::Accumulator{F,N},s,w) where {F,N}
    v=a.query(s)
    terms=ntuple(k->w*v[k]-a.corrections[k],Val(N))
    next=ntuple(k->a.totals[k]+terms[k],Val(N))
    a.corrections=ntuple(k->(next[k]-a.totals[k])-terms[k],Val(N))
    a.totals=next
    nothing
end

# Same event construction and CM velocity as the frozen production measure.
# Only accumulation order differs: Kahan streaming replaces sorting/storage.
function stream!(acc,plan,bounds)
    groups=AS.angular_groups(plan.cos_relative,plan.phiw)
    lo,hi=bounds;p,pw,ei,ej=plan.p,plan.pw,plan.ei,plan.ej
    folded=plan.mi==plan.mj && ei==ej
    events=0
    @inbounds for i in eachindex(p),j in (folded ? (i:length(p)) : eachindex(p))
        momentum=pw[i]*pw[j]*p[i]^2*p[j]^2*plan.prefactor
        reverse_momentum=folded && i!=j ? pw[j]*pw[i]*p[j]^2*p[i]^2*plan.prefactor : 0.0
        for group in groups
            angular=0.0;reverse_angular=0.0
            for (a,b) in group.pairs
                angular+=plan.cw[a]*plan.cw[b]*plan.fi[i,a]*plan.fj[j,b]
                if folded && i!=j
                    reverse_angular+=plan.cw[a]*plan.cw[b]*plan.fi[j,a]*plan.fj[i,b]
                end
            end
            combined=momentum*angular+reverse_momentum*reverse_angular
            combined==0 && continue
            for k in eachindex(group.cosines)
                s=plan.mi^2+plan.mj^2+2.0*(ei[i]*ej[j]-p[i]*p[j]*group.cosines[k])
                lo<s<hi || continue
                v=AS._relative_velocity_for_mode(:cm,s,plan.mi,plan.mj,ei[i],ej[j])
                0<v<=2 || continue
                w=combined*group.weights[k]*v
                w==0 && continue
                isfinite(w) && w>0 || error("invalid thermal event weight")
                visit!(acc,s,w);events+=1
            end
        end
    end
    events>0 || error("empty thermal support")
    events
end
function support(plan,bounds)
    acc=Support(Inf,-Inf)
    n=stream!(acc,plan,bounds)
    (lower=acc.lower,upper=acc.upper,events=n)
end
function integrate(plan,bounds,query,::Val{N}) where {N}
    z=ntuple(_->0.0,Val(N))
    acc=Accumulator(query,z,z)
    n=stream!(acc,plan,bounds)
    all(isfinite,acc.totals) || error("nonfinite thermal integral")
    (values=acc.totals,events=n)
end

# All three inputs interpolate raw sigma linearly. A union of their knots
# permits one exact segment lookup, preserving each piecewise linear curve.
struct RawComparison
    nodes::Vector{Float64}
    values::Vector{NTuple{3,Float64}}
end
function RawComparison(current,previous,shifted)
    tables=(current,previous,shifted)
    nodes=sort!(unique(vcat((t.nodes for t in tables)...)))
    RawComparison(nodes,[map(t->t(s),tables) for s in nodes])
end
@inline function (q::RawComparison)(s)
    first(q.nodes)<=s<=last(q.nodes) || error("outer event outside checked sigma support")
    i=searchsortedfirst(q.nodes,s)
    v=if q.nodes[i]==s
        @inbounds q.values[i]
    else
        @inbounds f=(s-q.nodes[i-1])/(q.nodes[i]-q.nodes[i-1])
        @inbounds ntuple(k->(1-f)*q.values[i-1][k]+f*q.values[i][k],Val(3))
    end
    (v[1],abs(v[1]-v[2]),abs(v[1]-v[3]))
end

function refined_grid(original,level;shifted=false)
    grid=copy(original);pieces=2^level
    fractions=shifted ? [(k+sqrt(2)-1)/pieces for k in 0:pieces-1] :
        [k/pieces for k in 1:pieces-1]
    for (a,b) in zip(original[1:end-1],original[2:end]),f in fractions
        a<a+(b-a)*f<b && push!(grid,a+(b-a)*f)
    end
    sort!(unique!(grid))
end

function run_job(row,channel,profile)
    started=time();job=String(row.point_id)*"__"*String(channel)
    jobdir=joinpath(OUTPUT,job);mkpath(jobdir)
    origin=only(filter(isdir,[joinpath(CASE,r,job) for r in S_RUN_DIRS]))
    input_table=joinpath(origin,"reference_table.csv")
    oldstats=JSON3.read(read(joinpath(origin,"completed.json"),String)).metrics
    state=frozen_state(row,profile);qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    problems=map(CONFIGS) do cfg
        p,w=AS._cached_semi_infinite_momentum_grid(cfg.p,profile.tau_integral.scale_inv_fm)
        AS._prepare_thermal_cache_problem(channel,qp,tp,kc;p_nodes=cfg.p,
            angle_nodes=cfg.c,phi_nodes=cfg.phi,p_grid=p,p_w=w,
            n_sigma_points=profile.tau_integral.t_nodes,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm)
    end
    normalization=problems[1].plan.prefactor
    plans=map(q->merge(q.plan,(prefactor=normalization,)),problems)
    bounds=problems[1].context.bounds
    all(p->p.context.bounds==bounds,problems) || error("outer physical domains differ")
    supports=map(p->support(p,bounds),plans)
    a,b=minimum(s.lower for s in supports),maximum(s.upper for s in supports)
    s0=first(bounds)
    old=collect(CSV.File(input_table))
    values=Dict(Float64(r.s_inv_fm2)=>Float64(r.sigma_fm2) for r in old)
    raw=s->CS.total_cross_section(channel,s,qp,tp,kc;n_points=profile.tau_integral.t_nodes)
    PC.reset_cache!()
    anchors=unique(round.(Int,range(1,length(old);length=9)))
    anchor_error=maximum(abs(raw(Float64(old[i].s_inv_fm2))/Float64(old[i].sigma_fm2)-1) for i in anchors)
    anchor_error<=1e-9 || error("frozen sigma replay differs")
    direct=s->get!(values,s) do
        v=raw(s)
        isfinite(v) && v>=0 || error("invalid direct outer reference sigma")
        v
    end
    original=sort!(unique(vcat(a,b,[Float64(r.s_inv_fm2) for r in old])))
    filter!(s->a<=s<=b,original)
    oldmin=Float64(first(old).s_inv_fm2)
    if a<oldmin
        append!(original,s0 .+ exp.(range(log(a-s0),log(oldmin-s0);length=65)))
        sort!(unique!(original))
    end
    previous=TS.WeightedSigmaTable(original,direct.(original),s0,s0)
    # Verify the streaming accumulator against the original fixed-grid result.
    replay=integrate(plans[1],bounds,s->(previous(s),),Val(1)).values[1]
    replay_error=replay/Float64(oldstats.reference_rate)-1
    abs(replay_error)<=1e-10 || error("streaming/reference replay differs: $replay_error")
    history=NamedTuple[]
    for level in 1:3
        PC.reset_cache!()
        grid=refined_grid(original,level)
        shifted_grid=refined_grid(original,level;shifted=true)
        table=TS.WeightedSigmaTable(grid,direct.(grid),s0,s0)
        shifted=TS.WeightedSigmaTable(shifted_grid,direct.(shifted_grid),s0,s0)
        query=RawComparison(table,previous,shifted)
        checks=map(zip(CONFIGS,plans)) do (cfg,plan)
            t0=time();out=integrate(plan,bounds,query,Val(3));rate,delta,shift=out.values
            metric=(level=level,configuration=cfg.name,p_nodes=cfg.p,cos_nodes=cfg.c,phi_nodes=cfg.phi,
                rate=rate,previous_l1_relative=delta/rate,shifted_l1_relative=shift/rate,
                events=out.events,stream_seconds=time()-t0)
            push!(history,metric)
            CSV.write(joinpath(jobdir,"history.csv"),history)
            println("OUTER ",job," ",cfg.name," level=",level," rate=",rate,
                " s_L1=",delta/rate," seconds=",round(time()-t0;digits=2));flush(stdout)
            metric
        end
        if all(c->max(c.previous_l1_relative,c.shifted_l1_relative)<=S_TARGET,checks)
            CSV.write(joinpath(jobdir,"reference_table.csv"),[(s_inv_fm2=s,sigma_fm2=v) for (s,v) in zip(table.nodes,table.sigma)])
            CSV.write(joinpath(jobdir,"shifted_table.csv"),[(s_inv_fm2=s,sigma_fm2=v) for (s,v) in zip(shifted.nodes,shifted.sigma)])
            base,one,two=checks
            change=two.rate/one.rate-1
            stats=(job_id=job,point_id=String(row.point_id),window=String(row.window),xi=Float64(row.xi),
                channel=String(channel),base_rate=base.rate,joint1_rate=one.rate,joint2_rate=two.rate,
                joint1_relative_to_base=one.rate/base.rate-1,joint2_relative_to_base=two.rate/base.rate-1,
                joint2_relative_to_joint1=change,outer_target=OUTER_TARGET,outer_converged=abs(change)<=OUTER_TARGET,
                max_s_l1=maximum(max(c.previous_l1_relative,c.shifted_l1_relative) for c in checks),
                stream_replay_relative=replay_error,sigma_anchor_relative=anchor_error,seconds=time()-started)
            write_json(joinpath(jobdir,"completed.json"),(status="completed_analysis",metrics=stats,
                normalization="Incoming density normalization fixed to production p32/c36/scale10",
                source_reference_table=input_table,source_reference_sha256=digest(input_table),
                files=Dict(n=>digest(joinpath(jobdir,n)) for n in readdir(jobdir))))
            return stats
        end
        previous=table
    end
    error("outer-grid sigma reference failed its unchanged L1 target")
end

function main()
    selected=Set(split(ARGS[4],','))
    rows=filter(r->String(r.point_id) in selected,collect(CSV.File(joinpath(INPUT,"source_rows.csv"))))
    length(rows)==length(selected) || error("requested outer point is absent")
    profile=Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    results=NamedTuple[]
    write_json(joinpath(OUTPUT,"run_config.json"),(purpose="outer quadrature at fixed formal backgrounds and density normalization",
        source_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",script_sha256=digest(@__FILE__),
        configurations=CONFIGS,points=sort!(collect(selected)),channels=CHANNELS,s_run_dirs=S_RUN_DIRS,s_target=S_TARGET,outer_target=OUTER_TARGET,
        threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads(),julia_version=string(VERSION)))
    for row in rows,channel in CHANNELS
        println("START OUTER ",row.point_id," ",channel);flush(stdout)
        push!(results,run_job(row,channel,profile))
        CSV.write(joinpath(OUTPUT,"channels.csv"),results)
        write_json(joinpath(OUTPUT,"progress.json"),(completed=length(results),expected=length(rows)*length(CHANNELS)))
        PC.reset_cache!();GC.gc()
    end
    write_json(joinpath(OUTPUT,"summary.json"),(status="completed_analysis",jobs=length(results),
        all_outer_converged=all(r.outer_converged for r in results),
        max_joint_change=maximum(abs(r.joint2_relative_to_joint1) for r in results),
        max_base_change=maximum(abs(r.joint2_relative_to_base) for r in results),
        finished_utc=string(now(UTC)),channels_sha256=digest(joinpath(OUTPUT,"channels.csv"))))
end
main()
