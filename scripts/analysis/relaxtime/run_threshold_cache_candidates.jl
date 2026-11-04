#!/usr/bin/env julia
using CSV, SHA
if !isdefined(Main,:ThresholdCacheCandidates)
    include(joinpath(@__DIR__,"ThresholdCacheCandidates.jl"))
end
using .ThresholdCacheCandidates, .ThresholdUpgradeCommon

function load_upgrade_case(ablationdir,xi,process)
    inputdir=joinpath(ablationdir,"inputs")
    row=retained_row(inputdir,xi); state=retained_state(row)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    s0=4qp.m.u^2
    file=joinpath(ablationdir,"knots","xi$(xi)_$(process).csv")
    source_rate=NaN
    if isfile(file)
        data=collect(CSV.File(file))
        raw=Dict(Float64(r.s)=>Float64(r.sigma) for r in data if Bool(r.stock_node))
        record=only(r for r in CSV.File(joinpath(ablationdir,"ablation.csv"))
                    if r.xi==xi && String(r.channel)==String(process) && String(r.method)=="stock_on")
        template=AS.CrossSectionCache(process,sort!(collect(keys(raw))),[raw[s] for s in sort!(collect(keys(raw)))])
        template.asym_enabled=true; template.asym_requested=true
        template.asym_s0=s0; template.asym_A=Float64(record.template_A)
        template.asym_taper_start=0.6; template.asym_taper_end=1.2
        stock=raw_cache(process,template.s_vals,[raw[s] for s in template.s_vals],template;subtract=true)
        source_rate=Float64(record.rate)
        ref=only(r for r in CSV.File(joinpath(ablationdir,"direct_reference.csv"))
                 if r.xi==xi && String(r.channel)==String(process))
        reference=Float64(ref.J256)
    else
        stock=AS.build_w0cdf_pchip_cache(process,qp,tp,kc;N=560,
            p_cutoff=Main.Constants_PNJL.Λ_inv_fm,n_sigma_points=24,
            threshold_subtraction=true,asym_window=0.6,asym_fit_min_points=8,
            asym_extra_points=10,sigma_cache_policy=:validated_anchored)
        raw=Dict(s=>AS._exact_cached_sigma_value(stock,s,stock.sigma_vals[i])
                 for (i,s) in enumerate(stock.s_vals))
        reference=direct_integral(process,state;order=256)
    end
    # Include the interface point explicitly so a composite table has no seam.
    direct(s)=CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
    xmin=s0*1e-10; xmax=0.6
    for x in (xmin,xmax)
        s=s0+x
        haskey(raw,s) || (raw[s]=direct(s))
    end
    grid=sort!(collect(keys(raw)))
    tail=AS.CrossSectionCache(process,grid,[raw[s] for s in grid])
    return (state=state,process=process,s0=s0,xmin=xmin,xmax=xmax,raw=raw,
            stock=stock,tail=tail,reference=reference,source_rate=source_rate)
end

function physical_seeds(case,estimates)
    seeds=[case.xmin,case.xmax]
    append!(seeds,[s-case.s0 for s in case.stock.s_vals if case.xmin<s-case.s0<case.xmax])
    append!(seeds,10.0 .^ range(log10(case.xmin),log10(case.xmax);length=12))
    for estimate in estimates, factor in (0.1,0.3,1.0,3.0,10.0)
        x=estimate.x_delta*factor
        isfinite(x) && case.xmin<x<case.xmax && push!(seeds,x)
    end
    return sort!(unique(seeds))
end

function table_weighted_integral(table::AdaptiveSigmaTable{R}) where {R}
    total=0.0
    for i in 1:length(table.nodes)-1
        a,b=table.nodes[i],table.nodes[i+1]
        fa,fb=table.sigma[i],table.sigma[i+1]
        total+=R===:x_sigma ? (b-a)*(a*fa+b*fb)/2 :
                                (b-a)*((2a+b)*fa+(a+2b)*fb)/6
    end
    return total
end

function run_candidate_study(ablationdir,outdir;tables_only=false)
    inputdir=joinpath(ablationdir,"inputs")
    start_run(outdir,inputdir,@__FILE__,(
        purpose="adaptive raw/x-sigma interpolation and local paired-PV sensitivity",
        ablation_dir=abspath(ablationdir),
        ablation_manifest_sha256=bytes2hex(sha256(read(joinpath(ablationdir,"manifest.json")))),
        tolerances=[1e-3,1e-4],weight="x=s-sth; independent direct J and full 5D rate checks",
        tables_only=tables_only,refinement_probe="nested Simpson quarter/midpoint values",
        local_x_max=0.6,physical_seeds="pi/mixed delta^2/c^2 plus original nodes and logarithmic probes",
        extra_checks="xi=-0.25,-0.2,0 at muB=0 alphaT=1"),
        extra_sources=[joinpath(@__DIR__,"ThresholdCacheCandidates.jl"),
                       joinpath(ROOT,"tests","unit","relaxtime","test_threshold_cache_candidates.jl")])
    results=NamedTuple[]; rate_rows=NamedTuple[]; scale_rows=NamedTuple[]
    cases=vcat([(xi,p) for xi in (0.37,0.36,0.38)
                for p in (:uubar_to_uubar,:uubar_to_ddbar)],
               [(xi,:uubar_to_uubar) for xi in (-0.25,-0.2,0.0)])
    parity_done=false
    for (xi,process) in cases
        println("ADAPTIVE CASE ",xi," ",process);flush(stdout)
        case=load_upgrade_case(ablationdir,xi,process)
        state=case.state; qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
        estimates=threshold_estimates(state)
        append!(scale_rows,[merge((xi=xi,channel=String(process)),e) for e in estimates])
        write_rows(joinpath(outdir,"threshold_scales.csv"),scale_rows)
        seeds=physical_seeds(case,estimates)
        known=Dict(s-case.s0=>v for (s,v) in case.raw)
        direct(s)=CS.total_cross_section(process,s,qp,tp,kc;n_points=24)
        tail(s)=evaluate_cache(case.tail,state,s)
        if !parity_done && !tables_only
            stock_query(s)=evaluate_cache(case.stock,state,s)
            parity=functional_rate(stock_query,process,state)
            isapprox(parity,case.source_rate;rtol=1e-10,atol=0) ||
                error("Callable consumer parity failed: $parity vs $(case.source_rate)")
            write_json(joinpath(outdir,"consumer_parity.json"),
                (xi=xi,channel=String(process),source=case.source_rate,
                 callable=parity,relative_error=parity/case.source_rate-1))
            parity_done=true
        end
        for representation in (:raw,:x_sigma), rtol in (1e-3,1e-4)
            timed=@timed build_adaptive_sigma(x->direct(case.s0+x),seeds;
                representation=representation,rtol=rtol,known_values=known)
            table=timed.value
            table.converged || error("Adaptive node budget exhausted")
            J=table_weighted_integral(table)
            push!(results,(xi=xi,channel=String(process),loop="legacy",
                representation=String(representation),rtol=rtol,nodes=length(table.nodes),
                direct_evaluations=table.evaluations,estimated_error=table.estimated_error,
                J=J,direct_reference=case.reference,J_relative_error=J/case.reference-1,
                build_seconds=timed.time))
            method="$(representation)_$(rtol)"
            write_rows(joinpath(outdir,"tables","xi$(xi)_$(process)_$(method).csv"),
                [(x=x,sigma=y) for (x,y) in zip(table.nodes,table.sigma)])
            write_rows(joinpath(outdir,"adaptive_integrals.csv"),results)
            println("TABLE ",method," n=",length(table.nodes)," new=",table.evaluations,
                    " Jerr=",J/case.reference-1);flush(stdout)
            # Compare both representations at the tighter target in the mixed
            # window, and the weighted candidate in the extra pi/iso cases.
            if !tables_only && rtol==1e-4 && (xi>0.3 || representation===:x_sigma)
                query=composite_sigma(table,case.s0,tail,direct)
                rate=functional_rate(query,process,state)
                push!(rate_rows,(xi=xi,channel=String(process),loop="legacy",
                    method=method,rate=rate))
                write_rows(joinpath(outdir,"candidate_rates.csv"),rate_rows)
            end
        end
        if xi>0.3
            plans=(u=PairLoopPlan(qp.m.u,qp.μ.u,tp;n_nodes=64),
                   s=PairLoopPlan(qp.m.s,qp.μ.s,tp;n_nodes=64))
            paired(s)=paired_cross_section(process,s,state,plans)
            # Local loop upgrade only. The tail remains the same raw legacy
            # cache; this is a controlled diagnostic, not a production hybrid.
            table=build_adaptive_sigma(x->paired(case.s0+x),seeds;rtol=1e-4)
            table.converged || error("Paired-loop adaptive table exhausted its budget")
            q,w=Main.RelaxTime.GaussLegendre.gauleg(0.0,sqrt(case.xmax),256)
            reference=sum(wi*2qi^3*paired(case.s0+qi^2) for (qi,wi) in zip(q,w))
            J=table_weighted_integral(table)
            push!(results,(xi=xi,channel=String(process),loop="paired_pv64",
                representation="x_sigma",rtol=1e-4,nodes=length(table.nodes),
                direct_evaluations=table.evaluations,estimated_error=table.estimated_error,
                J=J,direct_reference=reference,J_relative_error=J/reference-1,build_seconds=NaN))
            if !tables_only
                query=composite_sigma(table,case.s0,tail,paired)
                rate=functional_rate(query,process,state)
                push!(rate_rows,(xi=xi,channel=String(process),loop="paired_pv64",method="x_sigma_0.0001",rate=rate))
                println("PAIRED rate=",rate," Jerr=",J/reference-1);flush(stdout)
            end
            write_rows(joinpath(outdir,"adaptive_integrals.csv"),results)
            write_rows(joinpath(outdir,"candidate_rates.csv"),rate_rows)
            write_rows(joinpath(outdir,"tables","xi$(xi)_$(process)_paired_pv64.csv"),
                [(x=x,sigma=y) for (x,y) in zip(table.nodes,table.sigma)])
            println("PAIRED Jerr=",J/reference-1);flush(stdout)
        end
    end
    finish_run(outdir; integral_comparisons=length(results),rates=length(rate_rows),
        consumer_parity=parity_done,tables_only=tables_only)
    println("COMPLETE ",outdir)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (2,3) || error("Usage: run_threshold_cache_candidates.jl ABLATION_DIR NEW_OUTPUT_DIR [--tables-only]")
    length(ARGS)==3 && ARGS[3]!="--tables-only" && error("Unknown option $(ARGS[3])")
    run_candidate_study(abspath(ARGS[1]),abspath(ARGS[2]);tables_only=length(ARGS)==3)
end
