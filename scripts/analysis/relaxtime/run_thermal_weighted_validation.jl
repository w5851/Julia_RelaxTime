#!/usr/bin/env julia
using CSV, SHA
include(joinpath(@__DIR__,"ThermalWeightedSampling.jl"))
using .ThermalWeightedSampling, .ThermalAverageCandidates, .ThresholdUpgradeCommon

function fixed_reference(measure,s0,direct,features,thresholds,candidate;target=2e-4,candidate_target=1e-3)
    previous=nothing;history=NamedTuple[]
    for order in (96,192,384,768,1536,3072,6144)
        grid=thermal_seed_grid(measure,s0;thresholds=thresholds,order=order)
        append!(grid,features);sort!(unique!(grid))
        filter!(s->first(measure.nodes)<=s<=last(measure.nodes),grid)
        table=WeightedSigmaTable(grid,direct.(grid),s0,s0)
        rate=integrate_measure(measure,table)
        difference=previous===nothing ? Inf : let old_table=previous,new_table=table
            integrate_measure(measure,s->abs(new_table(s)-old_table(s)))
        end
        candidate_difference=integrate_measure(measure,s->abs(table(s)-candidate(s)))
        push!(history,(order=order,nodes=length(grid),rate=rate,
            previous_l1_relative=difference/max(rate,1e-30),candidate_l1_relative=candidate_difference/max(rate,1e-30)))
        if previous!==nothing && difference<=1e-12+target*rate && candidate_difference<=1e-12+candidate_target*rate
            return (converged=true,rate=rate,order=order,history=history)
        end
        previous=table
    end
    return (converged=false,rate=last(history).rate,order=last(history).order,history=history)
end

function weighted_validation(inputdir,outdir;scope="pilot",rtol=1e-3,shard=0,shards=1,limit=typemax(Int))
    scope in ("pilot","selected","qualification") || error("unknown validation scope")
    points=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    records=collect(CSV.File(joinpath(inputdir,"channel_plan.csv")))
    oraclefile=joinpath(dirname(inputdir),"anomaly_analysis_v1","rates.csv")
    oracle=Dict((String(r.point_id),Symbol(r.channel))=>r for r in CSV.File(oraclefile))
    jobs=Set((String(r.point_id),Symbol(r.channel)) for r in records if Bool(r.recompute))
    if scope=="pilot"
        pilot=Set([(0.0,1.0,0.37,:uubar_to_uubar),(0.0,1.0,0.37,:us_to_us),
            (450.0,1.0,0.26,:uubar_to_uubar),(900.0,1.0,-0.04,:uubar_to_uubar),
            (900.0,1.0,-0.04,:uubar_to_ddbar),(900.0,1.0,-0.02,:uubar_to_uubar),
            (900.0,1.0,-0.02,:uubar_to_ddbar)])
        jobs=Set((String(r.point_id),p) for r in points for p in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
            if (r.muB_MeV,r.alpha_T,r.xi,p) in pilot)
    elseif scope=="qualification"
        full=Set([(0.0,1.0,-0.22),(0.0,1.0,0.37),(0.0,1.2,-0.47),
            (450.0,1.0,0.26),(450.0,1.1,0.35),(900.0,1.0,-0.04)])
        for r in points
            (r.muB_MeV,r.alpha_T,r.xi) in full || continue
            union!(jobs,((String(r.point_id),p) for p in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
        end
    end
    shards>0 && 0<=shard<shards || error("invalid shard parameters")
    ordered=[pair for (index,pair) in enumerate(sort!(collect(jobs);by=identity))
        if mod(index-1,shards)==shard]
    resize!(ordered,min(length(ordered),limit))
    start_run(outdir,inputdir,@__FILE__,(purpose="uniform thermal-weighted cache qualification without per-background exceptions",
        scope=scope,rtol=rtol,p_nodes=128,angle_nodes=20,phi_nodes=36,t_nodes=24,
        shard=shard,shards=shards,expected_jobs=length(ordered),
        independent_reference_target=max(rtol,3e-4),historical_reference_used_for_sampling=false,
        final_observable_target=0.01,author_precision_choice="20261022 final observables 1 percent, cost priority",
        refinement="nested midpoint with two single shifted guards and a thermal supported probe",
        feature_backend="endpoint-subtracted equal-mass k0 PV; unchanged unequal-mass kernel",
        sigma_backend="historical production kernel, including original B0",
        oracle_path=abspath(oraclefile),oracle_sha256=bytes2hex(sha256(read(oraclefile)))),
        extra_sources=[joinpath(@__DIR__,"ThermalWeightedSampling.jl"),joinpath(@__DIR__,"ThermalAverageCandidates.jl"),
            joinpath(@__DIR__,"ThresholdCacheCandidates.jl")])
    results=NamedTuple[];failures=NamedTuple[];current_point="";state=nothing
    measures=Dict{Tuple,Any}()
    for (point,process) in ordered
        row=only(r for r in points if String(r.point_id)==point)
        if point!=current_point
            empty!(measures);GC.gc()
            state=retained_state(row);current_point=point
        end
        Main.RelaxTime.PolarizationCache.reset_cache!();GC.gc()
        println("START ",point," ",process);flush(stdout)
        started=time_ns()
        try
            bounds=ThermalAverageCandidates.process_bounds(process,state)
            incoming=CS.parse_particles_from_process(process)[1:2]
            key=(incoming...,bounds...)
            measured_plan=@timed get!(measures,key) do
                plan=prepare_thermal_plan(process,state)
                thermal_measure(plan,bounds)
            end
            measure=measured_plan.value;s0=bounds[1]
            println("MEASURE seconds=",measured_plan.time," events=",measure.events,
                " unique=",length(measure.nodes));flush(stdout)
            qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
            feature_t=@timed scattering_feature_seeds(process,state,measure,s0)
            feature=feature_t.value
            println("FEATURES seconds=",feature_t.time," count=",length(feature.features));flush(stdout)
            thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
            samples=Dict{Float64,Float64}()
            direct=s->get!(samples,s) do
                CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
            end
            built_t=@timed build_thermal_weighted_table(direct,measure,s0;rtol=rtol,
                thresholds=thresholds,extra_seeds=feature.seeds)
            built=built_t.value
            println("BUILT seconds=",built_t.time," calls=",built.evaluations,
                " converged=",built.converged);flush(stdout)
            id=point*"_"*String(process)
            write_rows(joinpath(outdir,"histories",id*".csv"),built.history)
            write_rows(joinpath(outdir,"features",id*".csv"),feature.features)
            write_rows(joinpath(outdir,"samples",id*".csv"),[(s=s,sigma=v) for (s,v) in sort!(collect(samples);by=first)])
            write_rows(joinpath(outdir,"tables",id*".csv"),[(s=s,sigma=v) for (s,v) in zip(built.table.nodes,built.table.sigma)])
            if !built.converged
                push!(failures,(point_id=point,channel=String(process),reason=built.reason,
                    error="adaptive thermal error target not reached"))
                write_rows(joinpath(outdir,"failures.csv"),failures)
                println("UNRESOLVED ",built.reason);flush(stdout)
                continue
            end
            table=validated_table(built)
            independent_table_rate=integrate_measure(measure,table)
            abs(independent_table_rate-built.rate)<=1e-11*max(built.rate,1e-30) ||
                error("cell integral and complete discrete measure differ")
            reference_source="retained_independent_reference";reference_order=0
            reference_t=@timed if haskey(oracle,(point,process))
                old=oracle[(point,process)]
                reference_order=Int(old.reference_fine_level)
                (converged=Bool(old.reference_converged),rate=Float64(old.reference_fine_rate))
            else
                reference_source="new_independent_fixed_design"
                ref=fixed_reference(measure,s0,direct,feature.seeds,thresholds,table;
                    candidate_target=max(rtol,3e-4))
                reference_order=ref.order
                write_rows(joinpath(outdir,"reference_histories",id*".csv"),ref.history)
                ref
            end
            reference=reference_t.value
            difference=abs(built.rate/reference.rate-1)
            passed=reference.converged && difference<=max(rtol,3e-4)
            push!(results,(point_id=point,channel=String(process),muB_MeV=Float64(row.muB_MeV),
                alpha_T=Float64(row.alpha_T),xi=Float64(row.xi),rate=built.rate,
                reference_rate=reference.rate,reference_order=reference_order,reference_source=reference_source,
                reference_relative_difference=difference,reference_converged=reference.converged,
                adaptive_converged=built.converged,validated=passed,
                estimated_relative_error=built.estimated_error/max(built.rate,1e-30),
                guard_phase=built.guard_phase,nodes=length(table.nodes),evaluations=built.evaluations,
                thermal_events=measure.events,unique_thermal_nodes=length(measure.nodes),
                feature_count=length(feature.features),measure_seconds=measured_plan.time,
                feature_seconds=feature_t.time,build_seconds=built_t.time,reference_seconds=reference_t.time,
                total_seconds=(time_ns()-started)*1e-9))
            write_rows(joinpath(outdir,"rates.csv"),results)
            println("RESULT pass=",passed," rel=",difference," rate=",built.rate,
                " nodes=",length(table.nodes)," calls=",built.evaluations," seconds=",last(results).total_seconds);flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,channel=String(process),reason="exception",
                error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
    end
    finish_run(outdir;attempted=length(ordered),successful=length(results),failures=length(failures),
        validated=count(r->r.validated,results),all_validated=isempty(failures) && all(r->r.validated,results))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)>=2 || error("Usage: run_thermal_weighted_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|selected|qualification] [RTOL SHARD SHARDS LIMIT]")
    weighted_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)>=3 ? ARGS[3] : "pilot",
        rtol=length(ARGS)>=4 ? parse(Float64,ARGS[4]) : 1e-3,
        shard=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 0,
        shards=length(ARGS)>=6 ? parse(Int,ARGS[6]) : 1,
        limit=length(ARGS)>=7 ? parse(Int,ARGS[7]) : typemax(Int))
end
