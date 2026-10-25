#!/usr/bin/env julia
# Independent numerical axes at the revised full-RS production backgrounds.
# No historical isotropic-coupling rate or table is used as a reference.
include(joinpath(@__DIR__, "run_transport_joint_validation.jl"))

const FULL_RS_KERNEL_CASES = Set([
    (0.0,1.0,-0.22,:uubar_to_uubar), (0.0,1.0,0.37,:uubar_to_uubar),
    (0.0,1.2,-0.47,:ssbar_to_ssbar), (450.0,1.0,0.26,:uubar_to_uubar),
    (450.0,1.1,0.35,:usbar_to_usbar), (900.0,1.0,-0.04,:uubar_to_uubar),
    (900.0,1.0,-0.04,:usbar_to_usbar), (900.0,1.0,-0.04,:ssbar_to_ssbar)])
const FULL_RS_OUTER_CASES = Set([
    (0.0,1.0,-0.22,:udbar_to_udbar), (450.0,1.1,0.35,:udbar_to_udbar),
    (900.0,1.0,-0.04,:udbar_to_udbar)])
const FULL_RS_CONTROLS = Set([
    (0.0,1.0,-0.22), (0.0,1.0,0.37), (0.0,1.2,-0.47),
    (450.0,1.0,0.26), (450.0,1.1,0.35), (900.0,1.0,-0.04),
    (0.0,1.0,0.0), (450.0,1.0,0.0), (900.0,1.1,0.0),
    (900.0,1.0,-0.003), (900.0,1.0,0.003), (0.0,1.2,-0.5)])

function full_rs_refined_a(state, nodes, angles)
    rt=Main.RelaxTime; qp,tp=state.quark_params,state.thermo_params
    p,w=AS._cached_semi_infinite_momentum_grid(nodes,5.0)
    av=rt.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=angles,use_aniso=true)
    gu=rt.EffectiveCouplings.calculate_G_from_A(av.u,qp.m.u)
    gs=rt.EffectiveCouplings.calculate_G_from_A(av.s,qp.m.s)
    kc=rt.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
    return (quark_params=merge(qp,(A=av,)),thermo_params=tp,K_coeffs=kc)
end

function full_rs_identity(state,rebuilt)
    keys(state.K_coeffs)==keys(rebuilt.K_coeffs) || error("full-RS coupling fields differ")
    a_exact=state.quark_params.A==rebuilt.quark_params.A
    original,actual=Tuple(state.K_coeffs),Tuple(rebuilt.K_coeffs)
    finite=all(isfinite,original) && all(isfinite,actual)
    relative=finite ? maximum(abs(a-b)/max(abs(a),abs(b),floatmin(Float64))
        for (a,b) in zip(original,actual)) : Inf
    # The @fastmath coupling helpers can round differently when inlined at
    # different call sites. A 24-background audit found identical A values and
    # at most one ULP in K. This is a roundoff guard, not a convergence budget.
    limit=8eps(Float64)
    return (A_exact=a_exact,coupling_max_relative=relative,coupling_roundoff_limit=limit,
        passed=a_exact && finite && relative<=limit)
end

function full_rs_outer_configs()
    base=(p_nodes=48,angle_nodes=36,phi_nodes=36,scale=10.0,density_nodes=32)
    changes=(("base",(;)), ("p96",(p_nodes=96,)), ("p192",(p_nodes=192,)),
        ("polar48",(angle_nodes=48,)), ("polar64",(angle_nodes=64,)),
        ("phi54",(phi_nodes=54,)), ("phi72",(phi_nodes=72,)),
        ("density256",(density_nodes=256,)),
        ("scale5_refined",(p_nodes=128,angle_nodes=48,phi_nodes=54,scale=5.0,density_nodes=256)),
        ("scale20_refined",(p_nodes=128,angle_nodes=48,phi_nodes=54,scale=20.0,density_nodes=256)),
        ("joint",(p_nodes=128,angle_nodes=48,phi_nodes=54,density_nodes=256)),
        ("joint_fine",(p_nodes=192,angle_nodes=64,phi_nodes=72,density_nodes=256)))
    return [(name=name,options=merge(base,options)) for (name,options) in changes]
end

function full_rs_grid_rate(process,state,options,t_nodes,prefix)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    p,pw=AS._cached_semi_infinite_momentum_grid(options.p_nodes,options.scale)
    dp,dw=AS._cached_semi_infinite_momentum_grid(options.density_nodes,10.0)
    plan=prepare_thermal_plan(process,state;p_nodes=options.p_nodes,
        angle_nodes=options.angle_nodes,phi_nodes=options.phi_nodes,
        p_grid=p,p_w=pw,density_p_grid=dp,density_p_w=dw)
    bounds=ThermalAverageCandidates.process_bounds(process,state)
    support=thermal_support(plan,bounds); s0=first(bounds)
    measure=ThermalMeasure([support.lower,support.upper],[1.0,1.0])
    features=scattering_feature_seeds(process,state,measure,s0)
    thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
    values=Dict{Float64,Float64}()
    direct=s->get!(values,s) do
        v=CS.total_cross_section(process,s,qp,tp,kc;n_points=t_nodes)
        isfinite(v) && v>=0 || throw(DomainError((s,v),"invalid full-RS reference cross section"))
        v
    end
    previous=nothing; history=NamedTuple[]
    # Same independent fixed-grid L1 protocol as the joint reference, with
    # explicit radial scale and density axes. Every actual event is streamed.
    for order in (192,384,768,1536,3072)
        grid=thermal_seed_grid(measure,s0;thresholds=thresholds,order=order)
        append!(grid,features.seeds); sort!(unique!(grid))
        filter!(s->support.lower<=s<=support.upper,grid)
        current=WeightedSigmaTable(grid,direct.(grid),s0,s0)
        if previous!==nothing
            query=compile_comparison(TableComparison(current,current,previous))
            elapsed=@timed streamed_integrals(plan,query,bounds,Val(4))
            _,rate,_,defect=elapsed.value.values
            isfinite(rate) && rate>0 || error("invalid full-RS reference rate")
            relative=defect/rate
            push!(history,(order=order,evaluations=length(values),rate=rate,
                l1_relative=relative,events=support.events,stream_seconds=elapsed.time))
            write_rows(prefix*"_history.csv",history)
            if relative<=2e-4
                write_rows(prefix*"_table.csv",[(s=s,sigma=v) for (s,v) in zip(current.nodes,current.sigma)])
                return (rate=rate,reference_l1=relative,order=order,evaluations=length(values))
            end
        end
        previous=current
    end
    error("full-RS reference exhausted its sigma resolution budget")
end

function full_rs_cheap_checks!(density_rows,pair_rows,moment_rows,a_rows,identity_rows,point,state,background_config)
    rt=Main.RelaxTime; qp,tp=state.quark_params,state.thermo_params
    rp,rw=AS._cached_semi_infinite_momentum_grid(256,5.0)
    density_configs=((name="workflow_p64_c8",p=64,c=8,scale=background_config.bg_limit,semi=false),
        (name="thermal_p32_c36",p=32,c=36,scale=10.0,semi=true),
        (name="thermal_p32_c48",p=32,c=48,scale=10.0,semi=true),
        (name="identity_p128_c64",p=128,c=64,scale=10.0,semi=true))
    for species in SPECIES
        mass,mu=AS.get_mass(species,qp),AS.get_mu(species,qp)
        reference=AS.number_density(species,mass,mu,tp.T,tp.Φ,tp.Φbar,0.0;
            p_grid=rp,p_w=rw,angle_nodes=1)/sqrt(1+tp.ξ)
        isfinite(reference) && reference>0 || error("invalid independent RS density")
        for config in density_configs
            p,w=config.semi ? AS._cached_semi_infinite_momentum_grid(config.p,config.scale) :
                rt.GaussLegendre.gauleg(0.0,config.scale,config.p)
            value=AS.number_density(species,mass,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ;
                p_grid=p,p_w=w,angle_nodes=config.c)
            relative=value/reference-1; target=config.name=="identity_p128_c64" ? 1e-9 : 1e-3
            push!(density_rows,(point_id=point,background_configuration=background_config.name,
                species=String(species),configuration=config.name,
                density=value,reference=reference,relative_change=relative,target=target,
                passed=isfinite(relative) && abs(relative)<=target))
        end
    end
    for name in ("thermal_p32_c36","thermal_p32_c48")
        relative=[r.relative_change for r in density_rows if r.point_id==point &&
            r.background_configuration==background_config.name && r.configuration==name]
        length(relative)==6 || error("incomplete independent density normalization check")
        bound=maximum(abs(1/((1+a)*(1+b))-1) for a in relative for b in relative)
        push!(pair_rows,(point_id=point,background_configuration=background_config.name,
            configuration=name,maximum_pair_rate_error=bound,passed=isfinite(bound) && bound<=1e-3))
    end
    configs=((name="production32_c16",p=32,c=16,scale=5.0),
        (name="p64_c16",p=64,c=16,scale=5.0), (name="p32_c32",p=32,c=32,scale=5.0),
        (name="joint64_c32",p=64,c=32,scale=5.0),
        (name="reference128_c48",p=128,c=48,scale=10.0),
        (name="reference256_c64",p=256,c=64,scale=10.0))
    grams=map(configs) do config
        p,w=AS._cached_semi_infinite_momentum_grid(config.p,config.scale)
        c,cw=rt.GaussLegendre.gauleg(-1.0,1.0,config.c)
        terminal_grams(qp,tp,p,w,c,cw)
    end
    reference=last(grams)
    for (config,value) in zip(configs,grams)
        eta=maximum(abs(getproperty(value,s)[1,1]/getproperty(reference,s)[1,1]-1) for s in SPECIES)
        sigma=maximum(abs(getproperty(value,s)[1,2]/getproperty(reference,s)[1,2]-1) for s in SPECIES)
        bulk=maximum(gram_relative_change(getproperty(value,s),getproperty(reference,s)) for s in SPECIES)
        push!(moment_rows,(point_id=point,background_configuration=background_config.name,
            configuration=config.name,eta_any_positive_tau_bound=eta,
            sigma_any_positive_tau_bound=sigma,bulk_any_bracket_bound=bulk,
            passed=all(isfinite,(eta,sigma,bulk)) && max(eta,sigma,bulk)<=0.01))
    end
    states=[full_rs_refined_a(state,p,c) for (p,c) in ((32,16),(64,32),(128,64))]
    actual=background_config.a_p==32 ? first(states) : states[2]
    push!(identity_rows,merge((point_id=point,background_configuration=background_config.name),
        full_rs_identity(state,actual)))
    for (config,value) in zip(("A32_c16","A64_c32","A128_c64"),states), flavor in (:u,:d,:s)
        av=getproperty(value.quark_params.A,flavor)
        ref=getproperty(last(states).quark_params.A,flavor)
        relative=av/ref-1
        push!(a_rows,(point_id=point,background_configuration=background_config.name,
            configuration=config,flavor=String(flavor),
            A=av,reference=ref,relative_change=relative,passed=isfinite(relative) && abs(relative)<=0.01))
    end
end

function full_rs_layer_validation(inputdir,outdir;layer="cheap",point_ids=nothing,channel_cases=nothing,extra_sources=String[])
    layer in ("cheap","kernel","outer") || throw(ArgumentError("layer must be cheap, kernel, or outer"))
    source=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    cases=layer=="kernel" ? FULL_RS_KERNEL_CASES : FULL_RS_OUTER_CASES
    if point_ids===nothing
        channel_cases===nothing || throw(ArgumentError("channel_cases requires point_ids"))
        controls=layer=="cheap" ? FULL_RS_CONTROLS : Set((x[1],x[2],x[3]) for x in cases)
        rows=filter(r->(r.muB_MeV,r.alpha_T,r.xi) in controls,source)
        length(rows)==length(controls) || error("missing or duplicated full-RS layer control points")
    else
        rows=requested_validation_rows(source,point_ids)
        if layer!="cheap"
            channel_cases===nothing && throw(ArgumentError("channel_cases is required for requested kernel/outer points"))
            cases=Set((String(id),Symbol(channel)) for (id,channel) in channel_cases)
            length(cases)==length(channel_cases) && Set(first.(cases))==Set(String.(point_ids)) ||
                throw(ArgumentError("channel_cases must uniquely cover every requested point"))
            required=Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
            all(last(case) in required for case in cases) || throw(ArgumentError("unsupported channel_cases channel"))
        elseif channel_cases!==nothing
            throw(ArgumentError("cheap checks do not take channel_cases"))
        end
    end
    kernel_configs=((name="base",a_p=32,a_c=16,t=24),
        (name="A64",a_p=64,a_c=32,t=24), (name="A128",a_p=128,a_c=64,t=24),
        (name="t48",a_p=32,a_c=16,t=48), (name="t96",a_p=32,a_c=16,t=96),
        (name="joint_A64_t48",a_p=64,a_c=32,t=48))
    outer_configs=full_rs_outer_configs()
    base_options=first(outer_configs).options
    background_config=merge(joint_config("candidate"),(p=48,))
    backgrounds=layer=="cheap" ? (background_config,joint_config("reference")) : (background_config,)
    start_run(outdir,inputdir,@__FILE__,(purpose="independent full-RS A, t, thermal, density and terminal axes",
        layer=layer,points=length(rows),requested_point_ids=point_ids,requested_channel_cases=channel_cases,
        backgrounds=backgrounds,
        kernel_configurations=kernel_configs,outer_configurations=outer_configs,
        expected_cases=layer=="cheap" ? length(rows) : length(cases),target=0.01,
        reference_sigma_target=2e-4,density_target=1e-3,rs_identity_target=1e-9,
        physical_conventions="full-RS A in polarization and couplings; cm velocity; unchanged B0 and finite sigma support",
        historical_rates_or_tables_reused=false,timing="diagnostic only; not a production benchmark"),
        extra_sources=vcat(extra_sources,[joinpath(@__DIR__,f) for f in ("run_transport_joint_validation.jl",
            "run_transport_background_tail_validation.jl","run_transport_background_validation.jl",
            "TerminalTransportQuadrature.jl","ThermalWeightedSampling.jl","ThermalAverageCandidates.jl",
            "ThresholdCacheCandidates.jl","StreamedThermalValidation.jl","TransportValidationSelection.jl")]))
    background_rows=NamedTuple[]; density_rows=NamedTuple[]; pair_rows=NamedTuple[]
    moment_rows=NamedTuple[]; a_rows=NamedTuple[]; identity_rows=NamedTuple[]
    rate_rows=NamedTuple[]; comparison_rows=NamedTuple[]; failures=NamedTuple[]
    for row in rows
        point=String(row.point_id); println("START ",layer," ",point); flush(stdout)
        try
          for background_config in backgrounds
            background=joint_background(row,background_config); state=background.state
            qp,tp=state.quark_params,state.thermo_params
            push!(background_rows,(point_id=point,configuration=background_config.name,
                T_fm=tp.T,muq_fm=qp.μ.u,xi=tp.ξ,
                m_u=qp.m.u,m_d=qp.m.d,m_s=qp.m.s,Phi=tp.Φ,Phibar=tp.Φbar,
                A_u=qp.A.u,A_d=qp.A.d,A_s=qp.A.s))
            write_rows(joinpath(outdir,"background.csv"),background_rows)
            if layer=="cheap"
                full_rs_cheap_checks!(density_rows,pair_rows,moment_rows,a_rows,identity_rows,point,state,background_config)
                write_rows(joinpath(outdir,"densities.csv"),density_rows)
                write_rows(joinpath(outdir,"pair_bounds.csv"),pair_rows)
                write_rows(joinpath(outdir,"moments.csv"),moment_rows)
                write_rows(joinpath(outdir,"A.csv"),a_rows)
                write_rows(joinpath(outdir,"identity.csv"),identity_rows)
                continue
            end
            for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
                key=point_ids===nothing ? (row.muB_MeV,row.alpha_T,row.xi,process) : (String(row.point_id),process)
                key in cases || continue
                values=Dict{String,Float64}()
                configs=layer=="kernel" ? kernel_configs : outer_configs
                for config in configs
                    Main.RelaxTime.PolarizationCache.reset_cache!(); GC.gc()
                    variant=layer=="kernel" ? full_rs_refined_a(state,config.a_p,config.a_c) : state
                    options=layer=="kernel" ? base_options : config.options
                    t_nodes=layer=="kernel" ? config.t : 24
                    prefix=joinpath(outdir,"references",point*"_"*String(process)*"_"*config.name)
                    elapsed=@timed full_rs_grid_rate(process,variant,options,t_nodes,prefix)
                    value=elapsed.value; values[config.name]=value.rate
                    push!(rate_rows,(point_id=point,channel=String(process),configuration=config.name,
                        rate=value.rate,reference_l1=value.reference_l1,reference_order=value.order,
                        evaluations=value.evaluations,seconds=elapsed.time))
                    write_rows(joinpath(outdir,"rates.csv"),rate_rows)
                    println("RATE ",layer," ",process," ",config.name," ",value.rate," seconds=",elapsed.time); flush(stdout)
                end
                if layer=="kernel"
                    pairs=(("base","A128"),("A64","A128"),("base","t96"),
                        ("t48","t96"),("base","joint_A64_t48"))
                else
                    pairs=Tuple((config.name,"joint_fine") for config in configs if config.name!="joint_fine")
                end
                for (candidate,reference) in pairs
                    c,r=values[candidate],values[reference]; relative=c/r-1
                    push!(comparison_rows,(point_id=point,channel=String(process),candidate_configuration=candidate,
                        reference_configuration=reference,candidate=c,reference=r,relative_change=relative,
                        passed=isfinite(relative) && abs(relative)<=0.01))
                end
                write_rows(joinpath(outdir,"comparisons.csv"),comparison_rows)
            end
          end
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,layer=layer,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",sprint(showerror,err)); flush(stdout)
        end
        GC.gc()
    end
    expected_rates=layer=="cheap" ? 0 : length(cases)*(layer=="kernel" ? length(kernel_configs) : length(outer_configs))
    expected_comparisons=layer=="cheap" ? 0 : length(cases)*(layer=="kernel" ? 5 : length(outer_configs)-1)
    expected_backgrounds=length(rows)*length(backgrounds)
    complete=isempty(failures) && length(background_rows)==expected_backgrounds &&
        length(rate_rows)==expected_rates && length(comparison_rows)==expected_comparisons &&
        (layer!="cheap" || (length(density_rows)==24expected_backgrounds &&
            length(pair_rows)==2expected_backgrounds && length(moment_rows)==6expected_backgrounds &&
            length(a_rows)==9expected_backgrounds && length(identity_rows)==expected_backgrounds))
    passed=complete && all(r->r.passed,comparison_rows) && all(r->r.passed,density_rows) &&
        all(r->r.passed,pair_rows) && all(r->r.passed,moment_rows) && all(r->r.passed,a_rows) &&
        all(r->r.passed,identity_rows)
    finish_run(outdir;layer=layer,points=length(rows),rates=length(rate_rows),comparisons=length(comparison_rows),
        backgrounds=length(background_rows),densities=length(density_rows),pair_bounds=length(pair_rows),
        moments=length(moment_rows),A_checks=length(a_rows),identity_checks=length(identity_rows),
        failures=length(failures),all_cases_complete=complete,layer_target_passed=passed,
        final_transport_qualification=false)
    return passed
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==3 || error("Usage: run_transport_full_rs_layer_validation.jl INPUT_DIR NEW_OUTPUT_DIR cheap|kernel|outer")
    full_rs_layer_validation(abspath(ARGS[1]),abspath(ARGS[2]);layer=ARGS[3]) || exit(1)
end
