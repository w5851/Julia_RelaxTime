#!/usr/bin/env julia
# Change only the numerical thermal upper limit in the same PNJL model.
include(joinpath(@__DIR__,"run_transport_background_validation.jl"))

function model_with_thermal_limit(model,limit)
    p=model.params
    parameters=NamedTuple{fieldnames(typeof(p))}(Tuple(getfield(p,k) for k in fieldnames(typeof(p))))
    return Models.PNJLModel(typeof(p)(;merge(parameters,(thermal_p_max_inv_fm=Float64(limit),))...))
end

function background_tail_validation(inputdir,outdir;scope="controls")
    source=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    sourcekeys=Set([(0.0,1.0,-0.22),(0.0,1.0,0.37),(0.0,1.2,-0.47),
        (450.0,1.0,0.26),(450.0,1.1,0.35),(900.0,1.0,-0.04),
        (900.0,1.0,-0.003),(900.0,1.0,0.003)])
    rows=scope=="all" ? source : filter(r->r.xi in (-0.5,0.0,0.5) ||
        (r.muB_MeV,r.alpha_T,r.xi) in sourcekeys,source)
    scope in ("controls","all") && !isempty(rows) || error("invalid background tail scope")
    configs=((name="thermal10",limit=10.0,p=64,c=16),
        (name="thermal15",limit=15.0,p=64,c=16),
        (name="thermal20",limit=20.0,p=64,c=16),
        (name="thermal25",limit=25.0,p=96,c=24),
        (name="thermal30",limit=30.0,p=128,c=32))
    start_run(outdir,inputdir,@__FILE__,(purpose="same-model thermal-tail convergence of background and bulk derivatives",
        scope=scope,points=length(rows),configurations=configs,relative_target=0.01,
        branch="each solve starts from the retained condensates; all shifts retained",
        tau="retained and held fixed to isolate background response",
        physical_constants="unchanged; only thermal_p_max_inv_fm and quadrature orders vary"),
        extra_sources=[joinpath(@__DIR__,"run_transport_background_validation.jl"),
            joinpath(@__DIR__,"TerminalTransportQuadrature.jl")])
    template=Models.create_model(:PNJL);tw=Models.TransportWorkflow
    p,pw=AS._cached_semi_infinite_momentum_grid(64,5.0)
    ct,cw=Main.RelaxTime.GaussLegendre.gauleg(-1.0,1.0,32)
    config=Main.RelaxTime.TransportCoefficients.TransportIntegrationConfig(p_grid=p,p_w=pw,cos_grid=ct,cos_w=cw)
    result_rows=NamedTuple[];failures=NamedTuple[]
    for row in rows
        point=String(row.point_id);println("START ",point);flush(stdout)
        try
            seed=retained_seed(template,row)
            tau=NamedTuple{SPECIES}(Tuple(Float64(getproperty(row,Symbol(:tau_,s))) for s in SPECIES))
            for c in configs
                model=model_with_thermal_limit(template,c.limit)
                cache=tw.WorkflowCache();cache.model_cache[:PNJL]=model
                elapsed=@timed begin
                    base=Main.EquilibriumFacade.solve_equilibrium_backend(row.T_fm,row.muq_fm;
                        xi=row.xi,p_num=c.p,t_num=c.c,seed_state=seed,model=model,models_residual_norm_max=1e-8)
                    Models.solve_transport_from_equilibrium(base,row.T_fm,row.muq_fm;
                        xi=row.xi,p_num=c.p,t_num=c.c,tau=tau,compute_tau=false,compute_bulk=true,
                        transport_config=config,workflow_cache=cache)
                end
                result=elapsed.value;bc=result.bulk_coeffs
                bc===nothing && error("background-tail bulk derivative failed")
                qp=Main.ParameterAdapters.normalize_quark_input(result.quark_params)
                tp=Main.ParameterAdapters.normalize_thermo_input(result.thermo_params)
                grams=terminal_grams(qp,tp,p,pw,ct,cw)
                tr=transport_from_grams(grams,qp,tp,tau;bulk_coeffs=bc)
                all(isfinite,(tr.eta,tr.sigma,tr.zeta,result.thermo_background.entropy)) || error("nonfinite tail result")
                mass_drift=maximum(abs(getproperty(qp.m,f)/Float64(getproperty(row,Symbol(:m_,f)))-1) for f in (:u,:d,:s))
                eq=result.equilibrium;bg=result.thermo_background
                push!(result_rows,(point_id=point,configuration=c.name,T_fm=Float64(row.T_fm),
                    muq_fm=Float64(row.muq_fm),muB_MeV=Float64(row.muB_MeV),alpha_T=Float64(row.alpha_T),xi=Float64(row.xi),
                    phi_u=eq.x_state[1],phi_d=eq.x_state[2],phi_s=eq.x_state[3],
                    m_u=qp.m.u,m_d=qp.m.d,m_s=qp.m.s,Phi=tp.Φ,Phibar=tp.Φbar,
                    entropy=bg.entropy,pressure=bg.pressure,energy=bg.energy,
                    n_u=result.densities.u,n_d=result.densities.d,n_s=result.densities.s,
                    n_ubar=result.densities.ubar,n_dbar=result.densities.dbar,n_sbar=result.densities.sbar,
                    eta=tr.eta,sigma=tr.sigma,zeta=tr.zeta,eta_over_s=tr.eta/bg.entropy,
                    zeta_over_s=tr.zeta/bg.entropy,sigma_over_T=tr.sigma/tp.T,
                    v_n_sq=bc.v_n_sq,dmuB_dT=bc.dμB_dT_sigma,
                    dMu_dT=bc.dM_dT[1],dMd_dT=bc.dM_dT[2],dMs_dT=bc.dM_dT[3],
                    dMu_dmuB=bc.dM_dμB[1],dMd_dmuB=bc.dM_dμB[2],dMs_dmuB=bc.dM_dμB[3],
                    retained_mass_relative_change=mass_drift,seconds=elapsed.time))
                write_rows(joinpath(outdir,"background.csv"),result_rows)
                println("LIMIT ",c.limit," mass_drift=",mass_drift," eta_s=",tr.eta/bg.entropy,
                    " zeta_s=",tr.zeta/bg.entropy);flush(stdout)
            end
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
        GC.gc()
    end
    finish_run(outdir;points=length(rows),comparisons=length(result_rows),failures=length(failures),
        all_points_complete=isempty(failures) && length(result_rows)==length(rows)*length(configs))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=3 || error("Usage: run_transport_background_tail_validation.jl INPUT_DIR NEW_OUTPUT_DIR [controls|all]")
    background_tail_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)==3 ? ARGS[3] : "controls")
end
