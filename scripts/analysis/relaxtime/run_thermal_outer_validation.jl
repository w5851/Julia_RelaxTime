#!/usr/bin/env julia
using CSV, SHA
include(joinpath(@__DIR__, "ThermalWeightedSampling.jl"))
include(joinpath(@__DIR__, "StreamedThermalValidation.jl"))
using .ThermalWeightedSampling, .ThermalAverageCandidates, .StreamedThermalValidation
using .ThresholdUpgradeCommon

function outer_grid_configs(scope)
    base = (p_nodes=128, angle_nodes=20, phi_nodes=36, scale=AS.DEFAULT_SEMI_INF_SCALE,
        density_nodes=AS.DEFAULT_SEMI_INF_NODES)
    configs = [(name="base", options=base)]
    for (name, options) in (("p96", (p_nodes=96,)), ("p192", (p_nodes=192,)), ("p256", (p_nodes=256,)),
        ("polar28", (angle_nodes=28,)), ("polar36", (angle_nodes=36,)),
        ("phi54", (phi_nodes=54,)), ("phi72", (phi_nodes=72,)),
        ("scale5", (scale=5.0,)), ("scale20", (scale=20.0,)),
        ("density256", (density_nodes=256,)),
        ("cost_candidate", (p_nodes=96, angle_nodes=16, phi_nodes=24,)),
        ("joint", (p_nodes=192, angle_nodes=28, phi_nodes=54, density_nodes=256,)))
        push!(configs, (name=name, options=merge(base, options)))
    end
    if scope=="outer_sensitive"
        for (name,options) in (("p96_polar28",(p_nodes=96,angle_nodes=28,)),
            ("joint_fine",(p_nodes=256,angle_nodes=36,phi_nodes=72,density_nodes=256,)),
            ("joint_ultra",(p_nodes=384,angle_nodes=48,phi_nodes=96,density_nodes=256,)))
            push!(configs,(name=name,options=merge(base,options)))
        end
    end
    return configs
end

function outer_plan(process, state, options)
    p, pw = AS._cached_semi_infinite_momentum_grid(options.p_nodes, options.scale)
    dp, dw = AS._cached_semi_infinite_momentum_grid(options.density_nodes, AS.DEFAULT_SEMI_INF_SCALE)
    return prepare_thermal_plan(process, state; p_nodes=options.p_nodes,
        angle_nodes=options.angle_nodes, phi_nodes=options.phi_nodes,
        p_grid=p, p_w=pw, density_p_grid=dp, density_p_w=dw)
end

function load_qualified_table(basedir, point, process)
    for shard in 0:1
        dir = joinpath(basedir, "thermal_weighted_qualification_s$(shard)_v1")
        path = joinpath(dir, "rates.csv")
        isfile(path) || continue
        rows = [r for r in CSV.File(path) if String(r.point_id)==point && Symbol(r.channel)==process]
        isempty(rows) && continue
        row = only(rows)
        row.validated && row.adaptive_converged || error("source table has not passed qualification")
        tablepath = joinpath(dir, "tables", point*"_"*String(process)*".csv")
        table = CSV.File(tablepath)
        return (nodes=Float64.(table.s), sigma=Float64.(table.sigma), rate=Float64(row.rate),
            path=abspath(tablepath), sha256=bytes2hex(sha256(read(tablepath))))
    end
    error("no qualified table for $(point) $(process)")
end

function extend_qualified_table(source, lower, upper, s0, direct)
    nodes = copy(source.nodes)
    if lower < first(nodes)
        append!(nodes, s0 .+ 10.0.^range(log10(lower-s0), log10(first(nodes)-s0); length=33))
        push!(nodes, lower)
    end
    upper > last(nodes) && append!(nodes, range(last(nodes), upper; length=9))
    sort!(unique!(nodes))
    existing = Dict(zip(source.nodes, source.sigma))
    values = [haskey(existing,s) ? existing[s] : direct(s) for s in nodes]
    return WeightedSigmaTable(nodes, values, s0, s0+0.6)
end

function outer_validation(inputdir, outdir; scope="pilot",shard=0,shards=1,retry_from=nothing)
    scope in ("pilot", "representatives", "outer_sensitive") || error("unknown scope")
    points = collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    selected = scope=="pilot" ? Set([(0.0,1.0,0.37)]) :
        Set([(0.0,1.0,-0.22),(0.0,1.0,0.37),(0.0,1.2,-0.47),
            (450.0,1.0,0.26),(450.0,1.1,0.35),(900.0,1.0,-0.04)])
    jobs = [(r,p) for r in points if (r.muB_MeV,r.alpha_T,r.xi) in selected
        for p in (scope=="pilot" ? (:uubar_to_uubar,:us_to_us) : Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU)]
    if scope=="outer_sensitive"
        filter!(job->last(job)===:udbar_to_udbar &&
            (first(job).muB_MeV,first(job).alpha_T,first(job).xi) in
            ((0.0,1.0,-0.22),(450.0,1.1,0.35),(900.0,1.0,-0.04)),jobs)
    end
    shards>0 && 0<=shard<shards || error("invalid shard parameters")
    jobs=[job for (i,job) in enumerate(jobs) if mod(i-1,shards)==shard]
    retry_evidence=nothing
    if retry_from!==nothing
        path=joinpath(retry_from,"failures.csv")
        sha256(read(joinpath(retry_from,"inputs","review_points.csv")))==
            sha256(read(joinpath(inputdir,"review_points.csv"))) || error("retry input differs")
        failed=Set((String(r.point_id),Symbol(r.channel)) for r in CSV.File(path))
        filter!(job->(String(first(job).point_id),last(job)) in failed,jobs)
        isempty(jobs) && error("no failed jobs match the requested outer retry")
        retry_evidence=(path=abspath(path),sha256=bytes2hex(sha256(read(path))))
    end
    configs = outer_grid_configs(scope)
    start_run(outdir,inputdir,@__FILE__,(purpose="independent outer-grid and joint thermal quadrature validation",
        scope=scope,shard=shard,shards=shards,jobs=length(jobs),grid_configs=configs,retry_from=retry_evidence,
        independent_axes_processes=["uubar_to_uubar","uubar_to_ddbar","ssbar_to_ssbar","ssbar_to_uubar"],
        all_process_grids=["base","cost_candidate","joint"],
        t_nodes=24,cache_target=1e-3,reference_target=2e-4,
        final_observable_target=0.01,reference_orders=[192,384,768,1536,3072],
        memory_strategy="stream every actual quadrature event without storing an event array"),
        extra_sources=[joinpath(@__DIR__,f) for f in ("ThermalWeightedSampling.jl",
            "ThermalAverageCandidates.jl","ThresholdCacheCandidates.jl","StreamedThermalValidation.jl")])
    results=NamedTuple[];failures=NamedTuple[];sources=NamedTuple[]
    basedir=dirname(inputdir)
    for (row,process) in jobs
        point=String(row.point_id)
        println("START ",point," ",process);flush(stdout)
        try
            source=load_qualified_table(basedir,point,process)
            push!(sources,(point_id=point,channel=String(process),path=source.path,sha256=source.sha256))
            write_rows(joinpath(outdir,"table_inputs.csv"),sources)
            state=retained_state(row);bounds=ThermalAverageCandidates.process_bounds(process,state);s0=bounds[1]
            used_configs=scope in ("pilot","outer_sensitive") || process in (:uubar_to_uubar,:uubar_to_ddbar,:ssbar_to_ssbar,:ssbar_to_uubar) ?
                configs : filter(c->c.name in ("base","cost_candidate","joint"),configs)
            plans=[outer_plan(process,state,c.options) for c in used_configs]
            supports=[thermal_support(p,bounds) for p in plans]
            lower=minimum(s.lower for s in supports);upper=maximum(s.upper for s in supports)
            seed_support=ThermalMeasure([lower,upper],[1.0,1.0]) # support only, never used as thermal weights
            feature=scattering_feature_seeds(process,state,seed_support,s0)
            qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
            thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
            Main.RelaxTime.PolarizationCache.reset_cache!();GC.gc()
            samples=Dict{Float64,Float64}()
            direct=s->get!(samples,s) do
                value=CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
                isfinite(value) && value>=0 || error("nonfinite direct cross section at s=$s")
                value
            end
            candidate=extend_qualified_table(source,lower,upper,s0,direct)
            previous=nothing;accepted=false
            for order in (192,384,768,1536,3072)
                grid=thermal_seed_grid(seed_support,s0;thresholds=thresholds,order=order)
                append!(grid,feature.seeds);sort!(unique!(grid))
                filter!(s->lower<=s<=upper,grid)
                reference=WeightedSigmaTable(grid,direct.(grid),s0,s0)
                if previous===nothing
                    previous=reference
                    continue
                end
                current=NamedTuple[]
                for (config,plan,support) in zip(used_configs,plans,supports)
                    query=TableComparison(candidate,reference,previous)
                    elapsed=@timed streamed_integrals(plan,query,bounds,Val(4))
                    c,r,cr,pr=elapsed.value.values
                    push!(current,(point_id=point,channel=String(process),muB_MeV=Float64(row.muB_MeV),
                        alpha_T=Float64(row.alpha_T),xi=Float64(row.xi),grid=String(config.name),
                        p_nodes=config.options.p_nodes,angle_nodes=config.options.angle_nodes,
                        phi_nodes=config.options.phi_nodes,scale=config.options.scale,
                        density_nodes=config.options.density_nodes,reference_order=order,
                        candidate_rate=c,reference_rate=r,cache_l1_relative=cr/max(r,1e-30),
                        reference_l1_relative=pr/max(r,1e-30),cache_validated=cr<=1e-12+1e-3r,
                        reference_converged=pr<=1e-12+2e-4r,events=support.events,
                        support_lower=support.lower,support_upper=support.upper,seconds=elapsed.time))
                    push!(results,last(current));write_rows(joinpath(outdir,"rates.csv"),results)
                    println("GRID ",config.name," order=",order," rate=",r,
                        " cache_l1=",last(current).cache_l1_relative,
                        " ref_l1=",last(current).reference_l1_relative," seconds=",elapsed.time);flush(stdout)
                    # Reject an under-resolved sigma reference on the cheapest
                    # common grid before repeating the expensive outer grids.
                    # Final acceptance still requires every requested grid.
                    if config.name=="base" && !last(current).reference_converged
                        break
                    end
                end
                # A failed reused candidate is retained as such. A converged
                # reference still permits an independent outer-grid diagnosis.
                if length(current)==length(used_configs) && all(r->r.reference_converged,current)
                    accepted=true
                    write_rows(joinpath(outdir,"references",point*"_"*String(process)*".csv"),
                        [(s=s,sigma=v) for (s,v) in zip(reference.nodes,reference.sigma)])
                    break
                end
                previous=reference
            end
            accepted || push!(failures,(point_id=point,channel=String(process),error="reference budget exhausted"))
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,channel=String(process),error=sprint(showerror,err,catch_backtrace())))
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
        write_rows(joinpath(outdir,"failures.csv"),failures)
        GC.gc()
    end
    finish_run(outdir;attempted=length(jobs),failures=length(failures),references_converged=isempty(failures))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=6 || error("Usage: run_thermal_outer_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|representatives|outer_sensitive SHARD SHARDS RETRY_DIR]")
    outer_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)>=3 ? ARGS[3] : "pilot",
        shard=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 0,shards=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 1,
        retry_from=length(ARGS)>=6 ? abspath(ARGS[6]) : nothing)
end
