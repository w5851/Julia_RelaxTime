#!/usr/bin/env julia
# Full native candidate versus independent, jointly refined numerical settings.
using CSV, JSON, SHA
include(joinpath(@__DIR__,"run_transport_background_tail_validation.jl"))
include(joinpath(@__DIR__,"ThermalWeightedSampling.jl"))
include(joinpath(@__DIR__,"StreamedThermalValidation.jl"))
if !isdefined(Main,:TransportValidationSelection)
    include(joinpath(@__DIR__,"TransportValidationSelection.jl"))
end
using .ThermalWeightedSampling, .ThermalAverageCandidates, .StreamedThermalValidation
using .TransportValidationSelection: requested_validation_rows

function joint_config(name)
    name=="candidate" && return (name=name,bg_limit=25.0,bg_p=96,bg_c=24,
        a_p=32,a_c=16,p=96,c=36,phi=36,t=24,tr_p=32,tr_c=16)
    name=="reference" && return (name=name,bg_limit=30.0,bg_p=128,bg_c=32,
        a_p=64,a_c=32,p=128,c=48,phi=54,t=48,tr_p=64,tr_c=32)
    error("unknown joint configuration")
end

function joint_background(row,config)
    template=Models.create_model(:PNJL)
    model=model_with_thermal_limit(template,config.bg_limit)
    seed=retained_seed(template,row)
    equilibrium=Main.EquilibriumFacade.solve_equilibrium_backend(row.T_fm,row.muq_fm;
        xi=row.xi,p_num=config.bg_p,t_num=config.bg_c,seed_state=seed,model=model,
        models_residual_norm_max=1e-8)
    params=Models.TransportWorkflow.build_equilibrium_params(equilibrium,row.T_fm,row.muq_fm;xi=row.xi)
    # build_equilibrium_params only creates the thermodynamic/chemical skeleton;
    # its mass slots are NaN until the equilibrium masses are supplied.
    masses=NamedTuple{(:u,:d,:s)}(Tuple(Float64.(equilibrium.masses)))
    all(m->isfinite(m) && m>0,masses) || error("joint equilibrium masses must be finite and positive")
    qp=merge(Main.ParameterAdapters.normalize_quark_input(params.quark_params),(m=masses,))
    tp=Main.ParameterAdapters.normalize_thermo_input(params.thermo_params)
    legacy_A=hasproperty(config,:legacy_A) && config.legacy_A
    p,w=legacy_A ? Main.RelaxTime.GaussLegendre.gauleg(0.0,20.0,config.a_p) :
        AS._cached_semi_infinite_momentum_grid(config.a_p,5.0)
    rt=Main.RelaxTime
    av=rt.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=config.a_c,use_aniso=true)
    kp,kw=legacy_A ? (rt.GaussLegendre.DEFAULT_MOMENTUM_NODES,rt.GaussLegendre.DEFAULT_MOMENTUM_WEIGHTS) : (p,w)
    iso=rt.AFieldBuilder.build_A_triplet(qp,tp,kp,kw;cos_nodes=config.a_c,use_aniso=false)
    coupling_A=legacy_A ? iso : av
    gu=rt.EffectiveCouplings.calculate_G_from_A(coupling_A.u,qp.m.u)
    gs=rt.EffectiveCouplings.calculate_G_from_A(coupling_A.s,qp.m.s)
    kc=rt.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    cache=Models.TransportWorkflow.WorkflowCache();cache.model_cache[:PNJL]=model
    return (equilibrium=equilibrium,workflow_cache=cache,
        state=(quark_params=merge(qp,(A=av,)),thermo_params=tp,K_coeffs=kc))
end

function joint_transport(row,background,config;tau_kwargs,check_density=true)
    p,w=hasproperty(config,:tr_limit) ? Main.RelaxTime.GaussLegendre.gauleg(0.0,config.tr_limit,config.tr_p) :
        AS._cached_semi_infinite_momentum_grid(config.tr_p,5.0)
    c,cw=Main.RelaxTime.GaussLegendre.gauleg(-1.0,1.0,config.tr_c)
    integration=Main.RelaxTime.TransportCoefficients.TransportIntegrationConfig(
        p_grid=p,p_w=w,cos_grid=c,cos_w=cw)
    result=Models.solve_transport_from_equilibrium(background.equilibrium,row.T_fm,row.muq_fm;
        xi=row.xi,p_num=config.bg_p,t_num=config.bg_c,compute_tau=true,compute_bulk=true,
        K_coeffs=background.state.K_coeffs,tau_kwargs=tau_kwargs,
        transport_config=integration,workflow_cache=background.workflow_cache)
    result.bulk_coeffs===nothing && error("joint bulk derivatives failed")
    values=(Tuple(result.tau)...,result.transport.eta,result.transport.zeta,result.transport.sigma,
        result.thermo_background.entropy)
    all(v->isfinite(v) && v>0,values) || error("joint output must be finite and positive")
    check_density || return (result=result,density_error=NaN)
    # The existing workflow uses its model's density grid (64/8), independently
    # of p_num/t_num. Check its tail and resolution at every new background.
    qp,tp=background.state.quark_params,background.state.thermo_params
    dp,dw=AS._cached_semi_infinite_momentum_grid(128,5.0)
    density_error=maximum(SPECIES) do s
        reference=AS.number_density(s,AS.get_mass(s,qp),AS.get_mu(s,qp),tp.T,tp.Φ,tp.Φbar,0.0;
            p_grid=dp,p_w=dw,angle_nodes=1)/sqrt(1+tp.ξ)
        abs(getproperty(result.densities,s)/reference-1)
    end
    density_error<=1e-3 || error("joint workflow density quadrature exceeds its 0.1 percent budget: $density_error")
    return (result=result,density_error=density_error)
end

function independent_joint_rate(process,state,config;history_path=nothing,table_path=nothing)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    plan=prepare_thermal_plan(process,state;p_nodes=config.p,angle_nodes=config.c,phi_nodes=config.phi)
    bounds=ThermalAverageCandidates.process_bounds(process,state)
    support=thermal_support(plan,bounds)
    measure=ThermalMeasure([support.lower,support.upper],[1.0,1.0]) # support only
    s0=first(bounds)
    features=scattering_feature_seeds(process,state,measure,s0)
    thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
    values=Dict{Float64,Float64}()
    direct=s->get!(values,s) do
        v=CS.total_cross_section(process,s,qp,tp,kc;n_points=config.t)
        isfinite(v) && v>=0 || throw(DomainError((s,v),"invalid independent joint sigma"))
        v
    end
    previous=nothing;history=NamedTuple[]
    for order in (192,384,768,1536,3072)
        grid=thermal_seed_grid(measure,s0;thresholds=thresholds,order=order)
        append!(grid,features.seeds);sort!(unique!(grid))
        filter!(s->support.lower<=s<=support.upper,grid)
        current=WeightedSigmaTable(grid,direct.(grid),s0,s0)
        if previous===nothing
            previous=current;continue
        end
        query=compile_comparison(TableComparison(current,current,previous))
        elapsed=@timed streamed_integrals(plan,query,bounds,Val(4))
        _,rate,_,defect=elapsed.value.values
        relative=defect/max(rate,1e-30)
        push!(history,(order=order,evaluations=length(values),rate=rate,l1_relative=relative,
            events=support.events,stream_seconds=elapsed.time))
        history_path===nothing || write_rows(history_path,history)
        if relative<=2e-4
            table_path===nothing || write_rows(table_path,
                [(s=s,sigma=v) for (s,v) in zip(current.nodes,current.sigma)])
            return (rate=rate,reference_l1=relative,order=order,evaluations=length(values))
        end
        previous=current
    end
    error("independent joint reference exhausted its sigma resolution budget")
end

function joint_observables(result)
    return merge(NamedTuple{Tuple(Symbol(:tau_,s) for s in SPECIES)}(Tuple(result.tau)),
        (eta_over_s=result.transport.eta/result.thermo_background.entropy,
         zeta_over_s=result.transport.zeta/result.thermo_background.entropy,
         sigma_over_T=result.transport.sigma/result.thermo_params.T))
end

function load_joint_reference(directory,inputdir,config)
    manifest=JSON.parsefile(joinpath(directory,"manifest.json"))
    status=JSON.parsefile(joinpath(directory,"status.json"))
    status["joint_target_passed"] || error("independent reference reuse requires a qualified joint run")
    manifest["parameters"]["configurations"]["reference"]==JSON.parse(JSON.json(config)) ||
        error("independent reference configuration changed")
    sha256(read(joinpath(directory,"inputs","review_points.csv")))==
        sha256(read(joinpath(inputdir,"review_points.csv"))) || error("independent reference inputs changed")
    for file in manifest["source_files"]
        path=replace(file["path"],'\\'=>'/')
        (startswith(path,"src/") || startswith(path,"config/models/") ||
            startswith(path,"config/physics/") || path in ("Project.toml","Manifest.toml") ||
            (startswith(path,"scripts/analysis/relaxtime/") &&
            path!="scripts/analysis/relaxtime/run_transport_joint_validation.jl")) || continue
        bytes2hex(sha256(read(joinpath(ROOT,file["path"]))))==file["sha256"] ||
            error("independent reference numerical source changed: $path")
    end
    # Candidate orders and reuse plumbing may change; the actual reference
    # background, quadrature and observable definitions must remain identical.
    current=replace(read(@__FILE__,String),"\r\n"=>"\n")
    original=replace(read(joinpath(directory,"source","scripts","analysis","relaxtime",
        "run_transport_joint_validation.jl"),String),"\r\n"=>"\n")
    for name in ("joint_background","joint_transport","independent_joint_rate","joint_observables")
        pattern=Regex("(?ms)^function "*name*"\\(.*?(?=^function |\\z)")
        a,b=match(pattern,current),match(pattern,original)
        a!==nothing && b!==nothing && strip(a.match)==strip(b.match) ||
            error("independent reference implementation changed: $name")
    end
    rate_path=joinpath(directory,"rates.csv");observable_path=joinpath(directory,"observables.csv")
    rates=filter(r->r.configuration=="reference",collect(CSV.File(rate_path)))
    observables=filter(r->r.configuration=="reference",collect(CSV.File(observable_path)))
    required=Set(String.(Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
    names=(Tuple(Symbol(:tau_,s) for s in SPECIES)...,:eta_over_s,:zeta_over_s,:sigma_over_T)
    points=Dict{String,Any}()
    for row in observables
        point=String(row.point_id);selected=filter(r->r.point_id==point,rates)
        length(selected)==length(required) && Set(String(r.channel) for r in selected)==required ||
            error("independent reference has an incomplete channel set: $point")
        all(r->isfinite(r.rate) && r.rate>0 && r.estimated_or_reference_l1<=2e-4,selected) ||
            error("independent reference channel qualification failed: $point")
        values=NamedTuple{names}(Tuple(Float64(getproperty(row,name)) for name in names))
        all(v->isfinite(v) && v>0,values) && row.density_relative_error<=1e-3 ||
            error("independent reference observable qualification failed: $point")
        points[point]=(rates=selected,observables=values,density_error=Float64(row.density_relative_error))
    end
    isempty(points) && error("independent reference contains no qualified points")
    evidence=(directory=abspath(directory),rates_sha256=bytes2hex(sha256(read(rate_path))),
        observables_sha256=bytes2hex(sha256(read(observable_path))),points=sort!(collect(keys(points))))
    return (points=points,evidence=evidence)
end

function joint_validation(inputdir,outdir;scope="pilot",shard=0,shards=1,reuse_dir=nothing,
    candidate_p::Union{Nothing,Int}=nothing,reference_dir=nothing,point_ids=nothing,extra_sources=String[])
    scope in ("pilot","representatives") || error("unknown joint scope")
    controls=Set([(0.0,1.0,-0.22),(0.0,1.0,0.37),(0.0,1.2,-0.47),
        (450.0,1.0,0.26),(450.0,1.1,0.35),(900.0,1.0,-0.04),
        (0.0,1.0,0.0),(450.0,1.0,0.0),(900.0,1.1,0.0),
        (900.0,1.0,-0.003),(900.0,1.0,0.003),(0.0,1.2,-0.5)])
    scope=="pilot" && (controls=Set([(0.0,1.0,-0.22)]))
    source=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    rows=point_ids===nothing ? filter(r->(r.muB_MeV,r.alpha_T,r.xi) in controls,source) :
        requested_validation_rows(source,point_ids)
    point_ids===nothing && length(rows)!=length(controls) &&
        error("joint input does not contain every requested control")
    shards>0 && 0<=shard<shards || error("invalid joint shard")
    rows=[r for (i,r) in enumerate(rows) if mod(i-1,shards)==shard]
    candidate_p===nothing || candidate_p>0 || error("candidate momentum order must be positive")
    candidate=joint_config("candidate")
    candidate_p===nothing || (candidate=merge(candidate,(p=candidate_p,)))
    configs=(candidate=candidate,reference=joint_config("reference"))
    independent_reference=reference_dir===nothing ? nothing :
        load_joint_reference(reference_dir,inputdir,configs.reference)
    reuse_evidence=nothing
    if reuse_dir!==nothing
        old=JSON.parsefile(joinpath(reuse_dir,"manifest.json"))
        status=JSON.parsefile(joinpath(reuse_dir,"status.json"))
        status["joint_target_passed"] || error("cannot reuse an unqualified joint pilot")
        old["parameters"]["configurations"]==JSON.parse(JSON.json(configs)) || error("joint reuse configuration changed")
        sha256(read(joinpath(reuse_dir,"inputs","review_points.csv")))==
            sha256(read(joinpath(inputdir,"review_points.csv"))) || error("joint reuse background inputs changed")
        for file in old["source_files"]
            path=replace(file["path"],'\\'=>'/')
            (startswith(path,"src/relaxtime/") || startswith(path,"src/models/pnjl_physics/") ||
                startswith(path,"config/models/") || startswith(path,"config/physics/") ||
                path=="src/models/workflow_apps/TransportWorkflow.jl") || continue
            bytes2hex(sha256(read(joinpath(ROOT,file["path"]))))==file["sha256"] ||
                error("joint reuse numerical source changed: $path")
        end
        path=joinpath(reuse_dir,"comparisons.csv")
        done=Set(String(r.point_id) for r in CSV.File(path))
        filter!(r->!(String(r.point_id) in done),rows)
        reuse_evidence=(path=abspath(path),sha256=bytes2hex(sha256(read(path))),points=sort!(collect(done)))
    end
    start_run(outdir,inputdir,@__FILE__,(purpose="joint native tau and transport convergence",
        scope=scope,shard=shard,shards=shards,points=length(rows),
        requested_point_ids=point_ids,configurations=configs,reused=reuse_evidence,
        reused_independent_reference=independent_reference===nothing ? nothing : independent_reference.evidence,
        target=0.01,cache_target=1e-3,reference_sigma_target=2e-4,
        physical_conventions="full-RS A in polarization and couplings; legacy_A retains the pre-switch isotropic reference; unchanged linear-xi B0 and physical sigma cutoff",
        branch="retained local condensate seed at every parameter point",
        timing="includes diagnostics and reference construction; not a production benchmark"),
        extra_sources=vcat(extra_sources,[joinpath(@__DIR__,f) for f in ("run_transport_background_tail_validation.jl",
            "run_transport_background_validation.jl","TerminalTransportQuadrature.jl",
            "ThermalWeightedSampling.jl","ThermalAverageCandidates.jl","ThresholdCacheCandidates.jl",
            "StreamedThermalValidation.jl","TransportValidationSelection.jl")]))
    rate_rows=NamedTuple[];observable_rows=NamedTuple[];comparison_rows=NamedTuple[]
    background_rows=NamedTuple[];failures=NamedTuple[]
    for row in rows
        point=String(row.point_id);outputs=Dict{String,Any}()
        println("START ",point);flush(stdout)
        try
            for config in configs
                if config.name=="reference" && independent_reference!==nothing &&
                    haskey(independent_reference.points,point)
                    source=independent_reference.points[point]
                    outputs[config.name]=source.observables
                    for r in source.rates
                        push!(rate_rows,(point_id=point,configuration=config.name,channel=String(r.channel),
                            rate=Float64(r.rate),estimated_or_reference_l1=Float64(r.estimated_or_reference_l1),
                            reference_order=Int(r.reference_order),evaluations=Int(r.evaluations),seconds=0.0))
                    end
                    push!(observable_rows,merge((point_id=point,configuration=config.name,
                        density_relative_error=source.density_error,transport_seconds=0.0),source.observables))
                    write_rows(joinpath(outdir,"rates.csv"),rate_rows)
                    write_rows(joinpath(outdir,"observables.csv"),observable_rows)
                    println("REUSED REFERENCE ",point," from ",reference_dir);flush(stdout)
                    continue
                end
                background_time=@timed joint_background(row,config)
                background=background_time.value;state=background.state
                qp,tp=state.quark_params,state.thermo_params
                push!(background_rows,(point_id=point,configuration=config.name,
                    m_u=qp.m.u,m_d=qp.m.d,m_s=qp.m.s,Phi=tp.Φ,Phibar=tp.Φbar,
                    A_u=qp.A.u,A_d=qp.A.d,A_s=qp.A.s,seconds=background_time.time))
                write_rows(joinpath(outdir,"background.csv"),background_rows)
                caches=Dict{Symbol,AS.ThermalWeightedCache}();rates=Dict{Symbol,Float64}()
                options=(p_nodes=config.p,angle_nodes=config.c,phi_nodes=config.phi,
                    n_sigma_points=config.t,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm)
                for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
                    Main.RelaxTime.PolarizationCache.reset_cache!();GC.gc()
                    id=point*"_"*String(process)*"_"*config.name
                    elapsed=@timed if config.name=="candidate"
                        cache=AS.build_thermal_weighted_cache(process,qp,tp,state.K_coeffs;options...)
                        caches[process]=cache
                        write_rows(joinpath(outdir,"tables",id*".csv"),
                            [(s=s,sigma=v) for (s,v) in zip(cache.table.nodes,cache.table.sigma)])
                        (rate=cache.rate,reference_l1=cache.estimated_error/cache.rate,order=0,
                            evaluations=cache.statistics.evaluations)
                    else
                        independent_joint_rate(process,state,config;
                            history_path=joinpath(outdir,"histories",id*".csv"),
                            table_path=joinpath(outdir,"tables",id*".csv"))
                    end
                    value=elapsed.value;rates[process]=value.rate
                    push!(rate_rows,(point_id=point,configuration=config.name,channel=String(process),
                        rate=value.rate,estimated_or_reference_l1=value.reference_l1,
                        reference_order=value.order,evaluations=value.evaluations,seconds=elapsed.time))
                    write_rows(joinpath(outdir,"rates.csv"),rate_rows)
                    println("RATE ",config.name," ",process," ",value.rate," seconds=",elapsed.time);flush(stdout)
                end
                tau_options=config.name=="candidate" ? merge(options,
                    (sigma_cache_policy=:thermal_weighted,cs_caches=caches,propagator_quark_params=qp)) :
                    (existing_rates=NamedTuple(rates),)
                elapsed=@timed joint_transport(row,background,config;tau_kwargs=tau_options)
                result=elapsed.value.result
                # This call consumes all 17 native caches through the actual
                # workflow/RelaxationTime APIs, including their identity checks.
                outputs[config.name]=joint_observables(result)
                push!(observable_rows,merge((point_id=point,configuration=config.name,
                    density_relative_error=elapsed.value.density_error,transport_seconds=elapsed.time),outputs[config.name]))
                write_rows(joinpath(outdir,"observables.csv"),observable_rows)
            end
            candidate,reference=outputs["candidate"],outputs["reference"]
            for name in keys(reference)
                c,r=getproperty(candidate,name),getproperty(reference,name)
                push!(comparison_rows,(point_id=point,observable=String(name),candidate=c,reference=r,
                    relative_change=c/r-1,passed=abs(c/r-1)<=0.01))
            end
            write_rows(joinpath(outdir,"comparisons.csv"),comparison_rows)
            println("POINT ",point," max_relative=",maximum(abs(r.relative_change) for r in comparison_rows if r.point_id==point));flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
        GC.gc()
    end
    finish_run(outdir;points=length(rows),comparisons=length(comparison_rows),failures=length(failures),
        joint_target_passed=isempty(failures) && length(comparison_rows)==9length(rows) && all(r.passed for r in comparison_rows))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=8 || error("Usage: run_transport_joint_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|representatives SHARD SHARDS REUSE_DIR|- CANDIDATE_P REFERENCE_DIR|-]")
    joint_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)>=3 ? ARGS[3] : "pilot",
        shard=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 0,shards=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 1,
        reuse_dir=length(ARGS)>=6 && ARGS[6]!="-" ? abspath(ARGS[6]) : nothing,
        candidate_p=length(ARGS)>=7 ? parse(Int,ARGS[7]) : nothing,
        reference_dir=length(ARGS)>=8 && ARGS[8]!="-" ? abspath(ARGS[8]) : nothing)
end
