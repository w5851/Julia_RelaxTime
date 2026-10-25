#!/usr/bin/env julia
using CSV, SHA, JSON
include(joinpath(@__DIR__, "ThresholdCacheCandidates.jl"))
include(joinpath(@__DIR__, "ThermalAverageCandidates.jl"))
if !isdefined(Main, :GroupedThermalProjection)
    include(joinpath(@__DIR__, "GroupedThermalProjection.jl"))
end
using .ThresholdCacheCandidates, .ThresholdUpgradeCommon, .ThermalAverageCandidates
using .GroupedThermalProjection

struct DirectSampleTable{F}
    nodes::Vector{Float64}
    sigma::Vector{Float64}
    direct::F
end

@inline function (table::DirectSampleTable)(s::Float64)
    xs, ys = table.nodes, table.sigma
    s < xs[1] && return table.direct(s)
    s > xs[end] && return 0.0
    i = searchsortedfirst(xs, s)
    i == 1 && return ys[1]
    @inbounds begin
        x0, x1 = xs[i-1], xs[i]
        a = (s-x0)/(x1-x0)
        return ys[i-1] + a*(ys[i]-ys[i-1])
    end
end

function local_table_integral(table, s0, xmax)
    edges = sort!(unique(vcat([s0, s0+xmax],filter(s->s0<s<s0+xmax, table.nodes))))
    z, w = Main.RelaxTime.GaussLegendre.gauleg(-1.0,1.0,12)
    value=0.0
    for (a,b) in zip(edges[1:end-1],edges[2:end]), (zi,wi) in zip(z,w)
        s=(a+b)/2+(b-a)*zi/2
        value+=(b-a)/2*wi*(s-s0)*table(s)
    end
    value
end

function build_anomaly_tables(process, state; rtol=1e-4,
    known_samples=Dict{Float64,Float64}(),reference_levels=(192,384))
    length(reference_levels)==2 && first(reference_levels)<last(reference_levels) ||
        throw(ArgumentError("two increasing reference levels are required"))
    all(n->n>=96 && ispow2(n÷96) && n%96==0,reference_levels) ||
        throw(ArgumentError("reference levels must be 96 times a power of two"))
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    i,j,c,d=CS.parse_particles_from_process(process)
    mi,mj,mc,md=(AS.get_mass(q,qp) for q in (i,j,c,d))
    cutoff=Main.Constants_PNJL.Λ_inv_fm
    _,s0,sup=AS._resolve_sigma_support_bounds(mi,mj,mc,md,cutoff,nothing,nothing)
    stock=AS.build_w0cdf_pchip_cache(process,qp,tp,kc;N=560,p_cutoff=cutoff,
        n_sigma_points=24,threshold_subtraction=true,sigma_cache_policy=:validated_anchored)
    values=Dict(s=>AS._exact_cached_sigma_value(stock,s,stock.sigma_vals[k])
        for (k,s) in enumerate(stock.s_vals) if s0<s<sup)
    for (s,value) in known_samples
        s0<s<sup || continue
        isfinite(value) && value>=0 || throw(DomainError((s,value),"invalid retained cross section"))
        # Recomputed stock knots take precedence; retained new samples belong
        # to the identical background and kernel, checked by the caller.
        haskey(values,s) || (values[s]=value)
    end
    new_calls=Ref(0)
    function direct(s)
        s<=s0 && return 0.0
        get!(values,s) do
            new_calls[]+=1
            CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
        end
    end
    xmin=max(s0*1e-10,64eps(s0)); xend=(sup-s0)*(1-1e-10)
    xmax=min(0.6,xend)
    # The old grid and logarithmic threshold seeds supply the adaptive design.
    seeds=sort!(unique(vcat([xmin,xmax],10.0.^range(log10(xmin),log10(xmax);length=32),
        [s-s0 for s in stock.s_vals if xmin<s-s0<xmax])))
    known=Dict(s-s0=>v for (s,v) in values if xmin<=s-s0<=xmax)
    local_table=build_adaptive_sigma(x->direct(s0+x),seeds;representation=:x_sigma,
        rtol=rtol,max_nodes=4096,known_values=known)
    local_table.converged || error("local adaptive budget exhausted")
    tail_seeds=sort!(unique(vcat([xmax,xend],collect(range(xmax,xend;length=65)),
        [s-s0 for s in stock.s_vals if xmax<s-s0<xend])))
    # Internal qqbar openings may lie inside the full support, not at its start.
    for threshold in (4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2)
        x=threshold-s0
        if xmax<x<xend
            for offset in (-1e-3,-1e-5,0.0,1e-5,1e-3)
                xmax<x+offset<xend && push!(tail_seeds,x+offset)
            end
        end
    end
    sort!(unique!(tail_seeds))
    tail_table=build_adaptive_sigma(x->direct(s0+x),tail_seeds;representation=:raw,
        rtol=rtol,max_nodes=8192,known_values=Dict(s-s0=>v for (s,v) in values))
    tail_table.converged || error("tail adaptive budget exhausted")
    candidate = s -> s<=s0 ? 0.0 : s<s0+xmin ? direct(s) :
        s<=s0+xmax ? local_table(clamp(s-s0,xmin,xmax)) :
        s<=s0+xend ? tail_table(clamp(s-s0,xmax,xend)) : 0.0
    references=DirectSampleTable[]
    reference_grids=Vector{Float64}[]
    for n in reference_levels
        grid=Float64[s0+xmin,s0+xmax,s0+xend]
        append!(grid,[s for s in stock.s_vals if s0+xmin<s<s0+xend])
        # Independent fixed design, nested between levels, covering the full support.
        for exponent in 0:trailing_zeros(n÷96)
            level=96*2^exponent
            append!(grid,s0 .+ 10.0.^range(log10(xmin),log10(xmax);length=level))
            append!(grid,[s0+xmax*(k/level)^2 for k in 1:level])
            append!(grid,range(s0+xmax,s0+xend;length=2level+1))
        end
        grid=sort!(unique(grid))
        push!(reference_grids,grid)
        push!(references,DirectSampleTable(grid,direct.(grid),direct))
    end
    return (;stock,local_table,tail_table,candidate,references,direct,s0,xmin,xmax,xend,
        new_calls=new_calls[],values)
end

function anomaly_main(inputdir,outdir; shard=0,shards=1,limit=typemax(Int),
    group_angles=false,completed_runs=String[])
    points=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    sort!(points;by=r->(r.muB_MeV==0.0 && r.alpha_T==1.0 && r.xi==0.37 ? 0 : 1,
        Float64(r.muB_MeV),Float64(r.alpha_T),Float64(r.xi)))
    plan=collect(CSV.File(joinpath(inputdir,"channel_plan.csv")))
    previous=Set{Tuple{String,Symbol}}()
    resumed_sources=NamedTuple[]
    for run in completed_runs
        for file in ("review_points.csv","channel_plan.csv")
            read(joinpath(run,"inputs",file))==read(joinpath(inputdir,file)) ||
                error("Continuation input differs: $run/$file")
        end
        rates_file=joinpath(run,"rates.csv")
        for r in CSV.File(rates_file)
            push!(previous,(String(r.point_id),Symbol(r.channel)))
        end
        push!(resumed_sources,(directory=abspath(run),rates_sha256=bytes2hex(sha256(read(rates_file)))))
    end
    projector=group_angles ? grouped_project_linear_rates : project_linear_rates
    start_run(outdir,inputdir,@__FILE__,(
        purpose="selected mode-A non-first-order anomalies: original backgrounds, adaptive full support and independent fixed designs",
        selection_file=joinpath(inputdir,"selection.json"),selection_sha256=bytes2hex(sha256(read(joinpath(inputdir,"selection.json")))),
        shard=shard,shards=shards,limit=limit,p_nodes=128,angle_nodes=20,phi_nodes=36,
        cache_rtol=1e-4,reference_levels=[192,384],reference_rate_target=2e-4,
        original_replay_target=1e-8,background_solve=false,production_write=false),
        extra_sources=[joinpath(@__DIR__,"ThresholdCacheCandidates.jl"),
            joinpath(@__DIR__,"ThermalAverageCandidates.jl"),
            joinpath(@__DIR__,"prepare_threshold_anomaly_validation.py"),
            joinpath(@__DIR__,"GroupedThermalProjection.jl")])
    write_json(joinpath(outdir,"continuation.json"),(
        angular_projection=group_angles ? "exact_equal_geometry_grouping" : "full_discrete_sum",
        resumed_sources=resumed_sources,previous_cases=length(previous)))
    results=NamedTuple[];failures=NamedTuple[];completed=0
    skipped=0
    projection_checked=Set{Symbol}()
    for (index,row) in enumerate(points)
        mod(index-1,shards)==shard || continue
        Main.RelaxTime.PolarizationCache.reset_cache!()
        state=retained_state(row)
        records=filter(r->r.point_id==row.point_id && Bool(r.recompute),plan)
        channels=sort!(unique(Symbol(r.channel) for r in records))
        for process in channels
            if (String(row.point_id),process) in previous
                skipped+=1
                continue
            end
            completed>=limit && break
            id=String(row.point_id)*"_"*String(process)
            println("CASE ",id);flush(stdout)
            try
                measured=@timed build_anomaly_tables(process,state)
                built=measured.value
                println("BUILT seconds=",measured.time," calls=",built.new_calls,
                    " nodes=",length(built.local_table.nodes),"/",length(built.tail_table.nodes));flush(stdout)
                stock_query=let cache=built.stock, fixed_state=state
                    s->evaluate_cache(cache,fixed_state,s)
                end
                grid=sort!(unique(vcat(built.s0 .+ built.local_table.nodes,built.s0 .+ built.tail_table.nodes,
                    built.references[1].nodes,built.references[2].nodes)))
                thermal_plan=prepare_thermal_plan(process,state)
                bounds=ThermalAverageCandidates.process_bounds(process,state)
                timed=@timed projector(thermal_plan,grid,built.s0,stock_query,
                    (built.candidate,built.references[1],built.references[2]),(built.s0+built.xmax,0.0,0.0),bounds)
                println("PROJECTED seconds=",timed.time," bytes=",timed.bytes);flush(stdout)
                old,candidate,coarse,fine=timed.value
                write_json(joinpath(outdir,"partial",id*".json"),(
                    rates=timed.value,build_seconds=measured.time,thermal_seconds=timed.time,
                    return_types=[string(Core.Compiler.return_type(q,Tuple{Float64})) for q in
                        (stock_query,built.candidate,built.references[1],built.references[2])]))
                projection_error=NaN
                if !(process in projection_checked)
                    small_plan=prepare_thermal_plan(process,state;p_nodes=12,angle_nodes=4,phi_nodes=8)
                    small_projected=projector(small_plan,grid,built.s0,stock_query,
                        (built.candidate,built.references[1],built.references[2]),(built.s0+built.xmax,0.0,0.0),bounds)
                    direct_query_rates=shared_rates((stock_query,built.candidate,built.references[1],built.references[2]),
                        (process,process,process,process),state;p_nodes=12,angle_nodes=4,phi_nodes=8)
                    projection_error=maximum(abs(a/b-1) for (a,b) in zip(small_projected,direct_query_rates))
                    projection_error<=1e-10 || error("union-grid projection parity failed: $projection_error")
                    push!(projection_checked,process)
                end
                original_values=unique(Float64(r.rate) for r in records if Symbol(r.channel)==process)
                maximum(original_values)-minimum(original_values)<1e-10*max(maximum(original_values),1e-12) ||
                    error("different source rates for one channel")
                source=first(original_values)
                parity=source==0 ? old : old/source-1
                reference_change=coarse/fine-1
                adaptive_change=candidate/fine-1
                jcoarse=local_table_integral(built.references[1],built.s0,built.xmax)
                jfine=local_table_integral(built.references[2],built.s0,built.xmax)
                push!(results,(point_id=String(row.point_id),channel=String(process),window_ids=String(row.window_ids),
                    muB_MeV=Float64(row.muB_MeV),alpha_T=Float64(row.alpha_T),xi=Float64(row.xi),
                    source_rate=source,replayed_rate=old,adaptive_rate=candidate,dense192_rate=coarse,dense384_rate=fine,
                    original_relative_error=parity,reference_relative_change=reference_change,adaptive_relative_error=adaptive_change,
                    projection_relative_error=projection_error,
                    J192=jcoarse,J384=jfine,J_reference_relative_change=jcoarse/jfine-1,
                    s0=built.s0,local_nodes=length(built.local_table.nodes),tail_nodes=length(built.tail_table.nodes),
                    direct_evaluations=built.new_calls,build_seconds=measured.time,thermal_seconds=timed.time,
                    thermal_allocated_bytes=timed.bytes,
                    original_replayed=abs(parity)<=1e-8,reference_converged=abs(reference_change)<=2e-4,
                    adaptive_validated=abs(adaptive_change)<=3e-4))
                write_rows(joinpath(outdir,"rates.csv"),results)
                write_rows(joinpath(outdir,"tables",id*".csv"),
                    [(s=s,sigma=v) for (s,v) in sort!(collect(built.values);by=first)])
                println("RESULT old=",parity," ref=",reference_change," adaptive=",adaptive_change,
                    " rate_change=",fine/source-1," seconds=",measured.time+timed.time);flush(stdout)
            catch err
                message=sprint(showerror,err,catch_backtrace())
                push!(failures,(point_id=String(row.point_id),channel=String(process),error=message))
                write_rows(joinpath(outdir,"failures.csv"),failures)
                println("FAILED ",id," ",message);flush(stdout)
            end
            completed+=1
        end
        completed>=limit && break
    end
    finish_run(outdir;attempted=completed,successful=length(results),failures=length(failures),
        skipped_previous=skipped,
        reference_converged=count(r->r.reference_converged,results),
        adaptive_validated=count(r->r.adaptive_validated,results),
        original_replayed=count(r->r.original_replayed,results))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    positional=filter(!=("--group-angles"),ARGS)
    length(positional)>=2 || error("Usage: run_threshold_anomaly_validation.jl INPUT_DIR NEW_OUTPUT_DIR [SHARD SHARDS LIMIT COMPLETED_RUNS...] [--group-angles]")
    anomaly_main(abspath(positional[1]),abspath(positional[2]);shard=length(positional)>=3 ? parse(Int,positional[3]) : 0,
        shards=length(positional)>=4 ? parse(Int,positional[4]) : 1,
        limit=length(positional)>=5 ? parse(Int,positional[5]) : typemax(Int),
        group_angles="--group-angles" in ARGS,
        completed_runs=length(positional)>=6 ? abspath.(positional[6:end]) : String[])
end
