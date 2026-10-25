#!/usr/bin/env julia
# Independent t and loop checks through the real scattering amplitude.
include(joinpath(@__DIR__,"LoopQuadratureReference.jl"))
include(joinpath(LOOP_REFERENCE_ROOT,"scripts","analysis","relaxtime","ThermalWeightedSampling.jl"))
include(joinpath(LOOP_REFERENCE_ROOT,"scripts","analysis","relaxtime","CrossSectionLoopProbe.jl"))
using CSV, SHA
using .LoopQuadratureReference, .ThermalWeightedSampling, .CrossSectionLoopProbe
using .ThermalAverageCandidates, .ThresholdUpgradeCommon

mutable struct ReferencePolarizationBackend
    rtol::Float64
    loops::Dict{Tuple,ComplexF64}
    max_loop_error::Float64
    evaluations::Int
end
ReferencePolarizationBackend(;rtol=2e-9)=ReferencePolarizationBackend(rtol,Dict{Tuple,ComplexF64}(),0.0,0)

function (backend::ReferencePolarizationBackend)(channel,k0,k,m1,m2,mu1,mu2,T,Phi,Phibar,xi,a1,a2,ns)
    lambda=k0+mu1-mu2
    tp=(T=T,Φ=Phi,Φbar=Phibar,ξ=xi)
    function loop_at(w)
        get!(backend.loops,(w,k,m1,m2,mu1,mu2,T,Phi,Phibar,xi)) do
            ref=loop_reference(w,k,m1,m2,mu1,mu2,tp;rtol=backend.rtol)
            backend.max_loop_error=max(backend.max_loop_error,ref.error)
            backend.evaluations+=1
            ref.value
        end
    end
    b0=loop_at(lambda)
    ns==1 && (b0=(b0+loop_at(-k0+mu1-mu2))/2)
    coefficient=k*k-lambda*lambda+(channel===:P ? (m1-m2)^2 : (m1+m2)^2)
    value=-Main.Constants_PNJL.N_color/(8pi^2)*(a1+a2+coefficient*b0)
    return (real(value),imag(value))
end

function qualified_kernel_source(basedir,point,process)
    for shard in 0:1
        dir=joinpath(basedir,"thermal_weighted_qualification_s$(shard)_v1")
        rows=[r for r in CSV.File(joinpath(dir,"rates.csv")) if String(r.point_id)==point && Symbol(r.channel)==process]
        isempty(rows) && continue
        row=only(rows)
        row.validated || error("baseline cache failed qualification")
        path=joinpath(dir,"tables",point*"_"*String(process)*".csv")
        data=CSV.File(path)
        return (nodes=Float64.(data.s),sigma=Float64.(data.sigma),rate=Float64(row.rate),
            path=abspath(path),sha256=bytes2hex(sha256(read(path))))
    end
    error("missing qualified baseline")
end

function variant_table(direct,measure,s0,seeds;target=2e-4)
    nodes=sort!(unique(seeds));values=Dict{Float64,Float64}();history=NamedTuple[]
    sample(s)=get!(values,s) do
        value=direct(s)
        isfinite(value) && value>=0 || throw(DomainError((s,value),"invalid kernel variant cross section"))
        value
    end
    previous=WeightedSigmaTable(nodes,sample.(nodes),s0,s0+0.6)
    for level in 1:4
        # Nested midpoints, then a non-midpoint audit at the next level.
        fraction=level==2 ? sqrt(2.0)-1 : 0.5
        added=[a+(b-a)*fraction for (a,b) in zip(nodes[1:end-1],nodes[2:end])]
        nodes=sort!(unique(vcat(nodes,added)))
        current=WeightedSigmaTable(nodes,sample.(nodes),s0,s0+0.6)
        rate=integrate_measure(measure,current)
        difference=let coarse_table=previous,fine_table=current
            integrate_measure(measure,s->abs(coarse_table(s)-fine_table(s)))
        end
        push!(history,(level=level,nodes=length(nodes),evaluations=length(values),rate=rate,
            l1_relative=difference/max(rate,1e-30)))
        if level>=2 && difference<=1e-12+target*rate
            return (converged=true,table=current,rate=rate,history=history,evaluations=length(values))
        end
        previous=current
    end
    return (converged=false,table=previous,rate=last(history).rate,history=history,evaluations=length(values))
end

function refined_a_state(state;couplings=false,nodes=32,cos_nodes=16)
    rt=Main.RelaxTime;qp,tp=state.quark_params,state.thermo_params
    p,w=AS._cached_semi_infinite_momentum_grid(nodes,5.0)
    av=rt.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=cos_nodes,use_aniso=true)
    kc=state.K_coeffs
    if couplings
        iso=rt.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=cos_nodes,use_aniso=false)
        gu=rt.EffectiveCouplings.calculate_G_from_A(iso.u,qp.m.u)
        gs=rt.EffectiveCouplings.calculate_G_from_A(iso.s,qp.m.s)
        kc=rt.EffectiveCouplings.calculate_effective_couplings(
            Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    end
    return (quark_params=merge(qp,(A=av,)),thermo_params=tp,K_coeffs=kc)
end

function kernel_validation(inputdir,outdir;scope="pilot",shard=0,shards=1,reuse_dir=nothing,variant_scope="all")
    scope in ("pilot","representatives") || error("unknown scope")
    variant_scope in ("all","a_convergence") || error("unknown kernel variant scope")
    points=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    cases=scope=="pilot" ? Set([(0.0,1.0,0.37,:uubar_to_uubar),(0.0,1.0,0.37,:us_to_us)]) :
        Set([(0.0,1.0,0.37,:uubar_to_uubar),(0.0,1.0,-0.22,:uubar_to_uubar),
            (0.0,1.2,-0.47,:ssbar_to_ssbar),(450.0,1.0,0.26,:uubar_to_uubar),
            (450.0,1.1,0.35,:usbar_to_usbar),(900.0,1.0,-0.04,:uubar_to_uubar),
            (900.0,1.0,-0.04,:usbar_to_usbar),(900.0,1.0,-0.04,:ssbar_to_ssbar)])
    jobs=[(r,p) for r in points for p in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
        if (r.muB_MeV,r.alpha_T,r.xi,p) in cases]
    shards>0 && 0<=shard<shards || throw(ArgumentError("invalid kernel shards"))
    base_variants=((name="t48",t_nodes=48,loop=:production,inputs=:retained),
        (name="t96",t_nodes=96,loop=:production,inputs=:retained),
        (name="loop_s",t_nodes=24,loop=:s,inputs=:retained),
        (name="loop_all",t_nodes=24,loop=:all,inputs=:retained),
        (name="joint_t48_loop_all",t_nodes=48,loop=:all,inputs=:retained),
        (name="A_refined",t_nodes=24,loop=:production,inputs=:a),
        (name="A_K_refined",t_nodes=24,loop=:production,inputs=:a_k))
    variants=variant_scope=="all" ?
        map(v->merge(v,(a_nodes=32,a_cos=16)),base_variants) :
        ((name="A_K_64",t_nodes=24,loop=:production,inputs=:a_k,a_nodes=64,a_cos=32),
         (name="A_K_128",t_nodes=24,loop=:production,inputs=:a_k,a_nodes=128,a_cos=64))
    reused=Tuple{String,Symbol}[]
    reused_variants=Set{Tuple{String,Symbol,String}}()
    reuse_evidence=nothing
    if reuse_dir!==nothing
        old_input=joinpath(reuse_dir,"inputs","review_points.csv")
        sha256(read(old_input))==sha256(read(joinpath(inputdir,"review_points.csv"))) ||
            error("reused kernel results have different background inputs")
        old_path=joinpath(reuse_dir,"rates.csv")
        old_rows=collect(CSV.File(old_path))
        for (r,p) in jobs
            found=filter(x->String(x.point_id)==String(r.point_id) && Symbol(x.channel)==p,old_rows)
            for x in found
                x.interpolation_converged && String(x.variant) in Set(v.name for v in variants) &&
                    push!(reused_variants,(String(r.point_id),p,String(x.variant)))
            end
            if all(v->(String(r.point_id),p,v.name) in reused_variants,variants)
                push!(reused,(String(r.point_id),p))
            end
        end
        reuse_evidence=(path=abspath(old_path),sha256=bytes2hex(sha256(read(old_path))),
            cases=reused,comparisons=collect(reused_variants))
        filter!(job->!((String(first(job).point_id),last(job)) in reused),jobs)
    end
    jobs=[job for (i,job) in enumerate(jobs) if mod(i-1,shards)==shard]
    start_run(outdir,inputdir,@__FILE__,(purpose="independent t, loop, and their joint changes through weighted thermal rates",
        scope=scope,shard=shard,shards=shards,variant_scope=variant_scope,variants=variants,reused_evidence=reuse_evidence,table_target=2e-4,
        loop_reference_rtol=2e-9,outer_grid="128/20/36",final_observable_target=0.01,
        refined_A="semi-infinite p32 scale5 polar16; couplings retain isotropic A convention",
        explicit_loop_probe="single process, restored cache keys; never a production backend"),
        extra_sources=[joinpath(@__DIR__,"LoopQuadratureReference.jl"),
            joinpath(LOOP_REFERENCE_ROOT,"benchmark","Project.toml"),joinpath(LOOP_REFERENCE_ROOT,"benchmark","Manifest.toml"),
            [joinpath(LOOP_REFERENCE_ROOT,"scripts","analysis","relaxtime",f) for f in
                ("CrossSectionLoopProbe.jl","ThermalWeightedSampling.jl","ThermalAverageCandidates.jl","ThresholdCacheCandidates.jl")]...])
    rows=NamedTuple[];failures=NamedTuple[];sources=NamedTuple[]
    for (row,process) in jobs
        point=String(row.point_id)
        println("START ",point," ",process);flush(stdout)
        try
            source=qualified_kernel_source(dirname(inputdir),point,process)
            push!(sources,(point_id=point,channel=String(process),path=source.path,sha256=source.sha256))
            write_rows(joinpath(outdir,"table_inputs.csv"),sources)
            state=retained_state(row);bounds=ThermalAverageCandidates.process_bounds(process,state)
            measure=thermal_measure(prepare_thermal_plan(process,state),bounds)
            feature=scattering_feature_seeds(process,state,measure,bounds[1])
            seeds=sort!(unique(vcat(source.nodes,feature.seeds)))
            filter!(s->first(measure.nodes)<=s<=last(measure.nodes),seeds)
            for variant in variants
                (point,process,variant.name) in reused_variants && continue
                Main.RelaxTime.PolarizationCache.reset_cache!();GC.gc()
                variant_state=variant.inputs===:retained ? state :
                    refined_a_state(state;couplings=variant.inputs===:a_k,
                        nodes=variant.a_nodes,cos_nodes=variant.a_cos)
                variant_seeds=if variant.inputs===:retained
                    seeds
                else
                    anchors=scattering_feature_seeds(process,variant_state,measure,bounds[1])
                    sort!(unique(vcat(seeds,anchors.seeds)))
                end
                backend=ReferencePolarizationBackend()
                probe=LoopProbe(backend;channels=variant.loop===:production ? :all : variant.loop)
                direct=if variant.loop===:production
                    let nt=variant.t_nodes,qp=variant_state.quark_params,tp=variant_state.thermo_params,kc=variant_state.K_coeffs,pr=process
                        s->CS.total_cross_section(pr,s,qp,tp,kc;n_points=nt)
                    end
                else
                    let nt=variant.t_nodes,pr=process,st=variant_state,backend_probe=probe
                        s->probe_cross_section(backend_probe,pr,s,st;n_points=nt)
                    end
                end
                elapsed=@timed variant_table(direct,measure,bounds[1],variant_seeds)
                result=elapsed.value
                id=point*"_"*String(process)*"_"*variant.name
                write_rows(joinpath(outdir,"histories",id*".csv"),result.history)
                write_rows(joinpath(outdir,"tables",id*".csv"),[(s=s,sigma=v) for (s,v) in zip(result.table.nodes,result.table.sigma)])
                push!(rows,(point_id=point,channel=String(process),variant=variant.name,
                    baseline_rate=source.rate,rate=result.rate,relative_change=result.rate/source.rate-1,
                    interpolation_converged=result.converged,last_l1_relative=last(result.history).l1_relative,
                    sigma_evaluations=result.evaluations,loop_evaluations=backend.evaluations,
                    max_loop_error=backend.max_loop_error,reference_seconds=elapsed.time))
                write_rows(joinpath(outdir,"rates.csv"),rows)
                println("VARIANT ",variant.name," converged=",result.converged,
                    " relative_change=",last(rows).relative_change," seconds=",elapsed.time);flush(stdout)
                result.converged || push!(failures,(point_id=point,channel=String(process),variant=variant.name,error="variant interpolation budget exhausted"))
            end
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,channel=String(process),variant="exception",error=sprint(showerror,err,catch_backtrace())))
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
        write_rows(joinpath(outdir,"failures.csv"),failures)
        GC.gc()
    end
    finish_run(outdir;cases=length(jobs),reused_cases=length(reused),
        reused_comparisons=length(reused_variants),comparisons=length(rows),failures=length(failures))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=7 || error("Usage: benchmark_thermal_kernel_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|representatives SHARD SHARDS REUSE_DIR all|a_convergence]")
    kernel_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)>=3 ? ARGS[3] : "pilot",
        shard=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 0,shards=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 1,
        reuse_dir=length(ARGS)>=6 ? abspath(ARGS[6]) : nothing,
        variant_scope=length(ARGS)>=7 ? ARGS[7] : "all")
end
