#!/usr/bin/env julia
# Same-branch background and A-field quadrature checks, with tau held fixed.
using CSV, NLsolve, StaticArrays
include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
if !isdefined(Main,:TerminalTransportQuadrature)
    include(joinpath(@__DIR__, "TerminalTransportQuadrature.jl"))
end
using .ThresholdUpgradeCommon, .TerminalTransportQuadrature

function retained_seed(model, row)
    target=SVector{3}(Float64(row.m_u),Float64(row.m_d),Float64(row.m_s))
    residual=x->Vector(Models.calculate_mass_vec(model,SVector{3}(x))-target)
    solved=nlsolve(residual,[-1.0,-1.0,-1.5];xtol=1e-12,ftol=1e-12)
    maximum(abs,residual(solved.zero))<1e-10 && all(<=(0.0),solved.zero) ||
        error("could not reconstruct the retained condensate branch")
    return SVector{5}(solved.zero...,Float64(row.Phi),Float64(row.Phibar))
end

function a_field_checks(qp,tp)
    configurations=((name="retained16_c4",n=16,c=4,limit=20.0,semi=false),
        (name="finite32_c8",n=32,c=8,limit=20.0,semi=false),
        (name="finite64_c16",n=64,c=16,limit=20.0,semi=false),
        (name="semi32_c16",n=32,c=16,limit=5.0,semi=true),
        (name="semi64_c32",n=64,c=32,limit=5.0,semi=true),
        (name="semi128_c64",n=128,c=64,limit=5.0,semi=true))
    rows=NamedTuple[]
    for config in configurations
        p,w=config.semi ? AS._cached_semi_infinite_momentum_grid(config.n,config.limit) :
            Main.RelaxTime.GaussLegendre.gauleg(0.0,config.limit,config.n)
        a=Main.RelaxTime.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=config.c,use_aniso=true)
        iso=Main.RelaxTime.AFieldBuilder.build_A_triplet(qp,tp,p,w;cos_nodes=config.c,use_aniso=false)
        for flavor in (:u,:d,:s)
            push!(rows,(configuration=config.name,flavor=String(flavor),
                A=getproperty(a,flavor),A_isotropic=getproperty(iso,flavor)))
        end
    end
    # Couplings historically use isotropic A on the separate 0..10 default grid.
    p=Main.RelaxTime.GaussLegendre.DEFAULT_MOMENTUM_NODES
    w=Main.RelaxTime.GaussLegendre.DEFAULT_MOMENTUM_WEIGHTS
    a=Main.RelaxTime.AFieldBuilder.build_A_triplet(qp,tp,p,w;use_aniso=false)
    for flavor in (:u,:d,:s)
        push!(rows,(configuration="retained_coupling_A",flavor=String(flavor),
            A=NaN,A_isotropic=getproperty(a,flavor)))
    end
    return rows
end

function background_validation(inputdir,outdir;scope="pilot")
    scope in ("pilot","all") || throw(ArgumentError("unknown background scope: $scope"))
    source=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    rows=scope=="all" ? source : filter(r->(r.muB_MeV,r.alpha_T,r.xi) in
        ((0.0,1.0,0.37),(900.0,1.0,-0.003),(0.0,1.2,-0.5)),source)
    isempty(rows) && error("empty background selection")
    configs=((name="base24_8",p=24,c=8),(name="p48",p=48,c=8),
        (name="polar16",p=24,c=16),(name="joint48_16",p=48,c=16),
        (name="reference96_32",p=96,c=32))
    start_run(outdir,inputdir,@__FILE__,(purpose="same-branch background, derivative, and A quadrature",
        scope=scope,points=length(rows),configurations=configs,relative_target=0.01,
        transport="fixed retained tau; semi-infinite p64 scale5 and polar32",
        seed="invert the current mass map to recover the retained condensates",
        branch="same-point local seed; reject a mass jump greater than 5 percent"),
        extra_sources=[joinpath(@__DIR__,"TerminalTransportQuadrature.jl")])
    model=Models.create_model(:PNJL)
    p,pw=AS._cached_semi_infinite_momentum_grid(64,5.0)
    ct,cw=Main.RelaxTime.GaussLegendre.gauleg(-1.0,1.0,32)
    TC=Main.RelaxTime.TransportCoefficients
    integration=TC.TransportIntegrationConfig(p_grid=p,p_w=pw,cos_grid=ct,cos_w=cw)
    background=NamedTuple[];moments=NamedTuple[];afields=NamedTuple[];failures=NamedTuple[]
    for row in rows
        point=String(row.point_id)
        println("START ",point);flush(stdout)
        try
            seed=retained_seed(model,row)
            tau=NamedTuple{SPECIES}(Tuple(Float64(getproperty(row,Symbol(:tau_,s))) for s in SPECIES))
            for config in configs
                elapsed=@timed Models.solve_gap_and_transport(Float64(row.T_fm),Float64(row.muq_fm);
                    xi=Float64(row.xi),p_num=config.p,t_num=config.c,seed_state=seed,
                    solver_backend=:models,models_residual_norm_max=1e-8,
                    compute_tau=false,tau=tau,compute_bulk=true,transport_config=integration)
                result=elapsed.value
                bc=result.bulk_coeffs
                bc===nothing && error("bulk derivatives failed for $(config.name)")
                qp=Main.ParameterAdapters.normalize_quark_input(result.quark_params)
                tp=Main.ParameterAdapters.normalize_thermo_input(result.thermo_params)
                mass_drift=maximum(abs(getproperty(qp.m,f)/Float64(getproperty(row,Symbol(:m_,f)))-1)
                    for f in (:u,:d,:s))
                mass_drift<=0.05 || error("background left retained branch: relative mass drift=$mass_drift")
                grams=terminal_grams(qp,tp,p,pw,ct,cw)
                tr=transport_from_grams(grams,qp,tp,tau;bulk_coeffs=bc)
                all(isfinite,(tr.eta,tr.sigma,tr.zeta,bc.s,bc.v_n_sq,bc.dμB_dT_sigma)) ||
                    error("nonfinite background transport/derivatives")
                bg=result.thermo_background;eq=result.equilibrium
                push!(background,(point_id=point,muB_MeV=Float64(row.muB_MeV),alpha_T=Float64(row.alpha_T),
                    xi=Float64(row.xi),T_fm=Float64(row.T_fm),muq_fm=Float64(row.muq_fm),
                    configuration=config.name,p_nodes=config.p,polar_nodes=config.c,
                    phi_u=eq.x_state[1],phi_d=eq.x_state[2],phi_s=eq.x_state[3],
                    Phi=tp.Φ,Phibar=tp.Φbar,m_u=qp.m.u,m_d=qp.m.d,m_s=qp.m.s,
                    entropy=bg.entropy,pressure=bg.pressure,energy=bg.energy,
                    n_u=result.densities.u,n_d=result.densities.d,n_s=result.densities.s,
                    n_ubar=result.densities.ubar,n_dbar=result.densities.dbar,n_sbar=result.densities.sbar,
                    eta=tr.eta,sigma=tr.sigma,zeta=tr.zeta,
                    eta_over_s=tr.eta/bg.entropy,zeta_over_s=tr.zeta/bg.entropy,sigma_over_T=tr.sigma/tp.T,
                    v_n_sq=bc.v_n_sq,dmuB_dT=bc.dμB_dT_sigma,
                    dMu_dT=bc.dM_dT[1],dMd_dT=bc.dM_dT[2],dMs_dT=bc.dM_dT[3],
                    dMu_dmuB=bc.dM_dμB[1],dMd_dmuB=bc.dM_dμB[2],dMs_dmuB=bc.dM_dμB[3],
                    retained_mass_relative_change=mass_drift,seconds=elapsed.time))
                for s in SPECIES,i in 1:3,j in i:3
                    push!(moments,(point_id=point,configuration=config.name,species=String(s),
                        i=i,j=j,moment=getproperty(grams,s)[i,j]))
                end
                if config.name in ("base24_8","reference96_32")
                    append!(afields,[merge((point_id=point,background=config.name),a) for a in a_field_checks(qp,tp)])
                end
                write_rows(joinpath(outdir,"background.csv"),background)
                write_rows(joinpath(outdir,"grams.csv"),moments)
                write_rows(joinpath(outdir,"afields.csv"),afields)
                println("GRID ",config.name," mass_drift=",mass_drift," eta_s=",tr.eta/bg.entropy,
                    " zeta_s=",tr.zeta/bg.entropy," seconds=",elapsed.time);flush(stdout)
            end
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
        GC.gc()
    end
    finish_run(outdir;points=length(rows),comparisons=length(background),failures=length(failures),
        all_points_complete=isempty(failures) && length(background)==length(rows)*length(configs))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=3 || error("Usage: run_transport_background_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|all]")
    background_validation(abspath(ARGS[1]),abspath(ARGS[2]);scope=length(ARGS)==3 ? ARGS[3] : "pilot")
end
