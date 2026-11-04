#!/usr/bin/env julia
# quote-remote only: SOURCE_ROOT INPUT_DIR OUTPUT_DIR STAGE SELECTION PROFILE
using StaticArrays, NLsolve
include(joinpath(@__DIR__,"run_transport_mott_probe.jl"))
include(joinpath(@__DIR__,"transport_mott_thresholds.jl"))
include(joinpath(@__DIR__,"transport_mott_kaon_endpoint.jl"))
include(joinpath(SOURCE_ROOT,"scripts","relaxtime","transport_sampling_scan_config.jl"))
include(joinpath(SOURCE_ROOT,"scripts","relaxtime","gap_transport_scan_phase_equilibrium.jl"))
const SP_OTHER=Models.TransportSamplingProfiles
const PHASE_OTHER=Main.GapTransportScanPhaseEquilibrium
const TC_OTHER=Main.RelaxTime.TransportCoefficients
const SOURCE_OTHER=collect(CSV.File(joinpath(INPUT_DIR,"source_rows.csv")))
const SPECIES_OTHER=(:u,:d,:s,:ubar,:dbar,:sbar)

function original_seed_other(model,row)
    target=SVector{3}(Float64(row.m_u),Float64(row.m_d),Float64(row.m_s))
    residual=x->Vector(Models.calculate_mass_vec(model,SVector{3}(x))-target)
    solved=nlsolve(residual,[-1.0,-1.0,-1.5];xtol=1e-12,ftol=1e-12)
    maximum(abs,residual(solved.zero))<1e-10 || error("original condensate reconstruction failed")
    SVector{5}(solved.zero...,Float64(row.Phi),Float64(row.Phibar))
end

function background_other(target,xi,profile,level;seed=nothing)
    ap,ac,p,c,limit=level=="base" ? (32,16,96,24,25.0) :
        level=="A64" ? (64,32,96,24,25.0) :
        level=="background_ref" ? (64,32,128,32,30.0) : error("unknown background level")
    model=SP_OTHER.profile_model(profile;thermal_p_max_inv_fm=limit)
    opts=(model=model,p_num=p,t_num=c,residual_norm_max=1e-10,max_iter=80)
    T,mu=Float64(target.T_MeV)/HC,Float64(target.muB_MeV)/(3HC)
    eq=PHASE_OTHER.solve_with_multiseed_governance(T,mu,Float64(xi),opts;continuation_seed=seed)
    eq!==nothing && eq.residual_norm<=1e-10 || error("background failed at xi=$xi")
    row=(point_id=String(target.curve_id),window=String(target.window),
        T_MeV=Float64(target.T_MeV),muB_MeV=Float64(target.muB_MeV),xi=Float64(xi),T_fm=T,muq_fm=mu,
        m_u=eq.masses[1],m_d=eq.masses[2],m_s=eq.masses[3],Phi=eq.x_state[4],Phibar=eq.x_state[5],
        phi_u=eq.x_state[1],phi_d=eq.x_state[2],phi_s=eq.x_state[3],
        background_residual=eq.residual_norm,omega=eq.omega)
    (row=row,equilibrium=eq,state=frozen_state(row,profile,ap,ac))
end

endpoint_other(state,meson,cut;n_nodes=128)=meson==:K ? kaon_endpoint_inverse(state;n_nodes=n_nodes) :
    diagnostic_threshold_inverse(state,meson,cut;n_nodes=n_nodes)

function locate_other()
    config=JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String))
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"base.toml"))
    locations,iterations,comparisons,polerows,candidates=Any[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    function checkpoint()
        for (name,rows) in (("iterations",iterations),("kernel_controls",comparisons),("poles",polerows),("pole_candidates",candidates))
            isempty(rows) || CSV.write(joinpath(OUTPUT_DIR,name*".csv"),rows)
        end
        write_json(joinpath(OUTPUT_DIR,"locations.json"),locations)
    end
    for target in config.targets
        meson,cut=Symbol(target.meson),Symbol(target.cut)
        original=only(r for r in SOURCE_OTHER if String(r.point_id)==String(target.point_left))
        seed=original_seed_other(SP_OTHER.profile_model(profile),original)
        levels=Dict{String,Any}()
        for level in ("base","A64","background_ref")
            cache=Dict{Float64,Any}()
            getbg(x)=get!(cache,Float64(x)) do
                background_other(target,x,profile,level;seed=seed)
            end
            for n_nodes in (level=="base" ? (128,256) : (128,))
                function scalar(x)
                    bg=getbg(x);v=endpoint_other(bg.state,meson,cut;n_nodes=n_nodes)
                    push!(iterations,(curve_id=String(target.curve_id),meson=String(meson),cut=String(cut),level=level,
                        n_nodes=n_nodes,xi=Float64(x),inverse_re=real(v),inverse_im=imag(v),
                        background_residual=bg.equilibrium.residual_norm))
                    isfinite(real(v)) && isfinite(imag(v)) && abs(imag(v))<1e-7 || error("not a closed threshold")
                    real(v)
                end
                # Never expand through the already reviewed first-order gaps.
                lo,hi=Float64(target.search_left),Float64(target.search_right)
                flo,fhi=scalar(lo),scalar(hi)
                flo*fhi<0 || error("endpoint bracket missing: $(target.curve_id) $level $flo $fhi")
                for _ in 1:45
                    hi-lo<=1e-9 && break
                    mid=(lo+hi)/2;fm=scalar(mid)
                    if flo*fm<=0
                        hi=mid;fhi=fm
                    else
                        lo=mid;flo=fm
                    end
                end
                key=level=="base"&&n_nodes==256 ? "pair_ref" : level
                levels[key]=(xi=(lo+hi)/2,xi_left=lo,xi_right=hi,inverse_left=flo,inverse_right=fhi,n_nodes=n_nodes)
                println("LOCATION ",target.curve_id," ",key," xi=",(lo+hi)/2);flush(stdout)
            end
            if level=="base"
                center=levels["base"].xi;bg=getbg(center);st=bg.state
                threshold=meson in MIXED ? (cut==:uu ? 2st.qp.m.u : 2st.qp.m.s) : MT.mott_threshold_mass(meson,st.qp)
                for offset in (-0.01,-1e-3,-1e-4,-1e-6,1e-6,1e-4,1e-3,0.01)
                    native=MP.meson_mass_equation(meson,threshold+offset,0.0,0.0,st.qp,st.tp,st.kc)
                    stable=meson==:K ? kaon_pair_inverse(st,threshold+offset;n_nodes=256) :
                        diagnostic_threshold_inverse(st,meson,cut;n_nodes=256,mass_offset=offset)
                    push!(comparisons,(curve_id=String(target.curve_id),meson=String(meson),xi=center,offset_inv_fm=offset,
                        native_re=real(native),native_im=imag(native),stable_re=real(stable),stable_im=imag(stable),
                        native_finite=isfinite(real(native))&&isfinite(imag(native))))
                end
                scan=[getbg(center+d) for d in (-0.004,-0.002,-0.0002,0.0002,0.002,0.004)
                    if Float64(target.segment_left)<center+d<Float64(target.segment_right)]
                for direction in ("forward","reverse")
                    previous=nothing
                    for b in (direction=="forward" ? scan : reverse(scan))
                        row=merge(b.row,(point_id=String(target.curve_id)*"_"*string(b.row.xi),))
                        cs,sm,previous=probe_roots(row,b.state,meson,direction,previous)
                        append!(candidates,cs);push!(polerows,sm)
                    end
                end
            end
        end
        push!(locations,(curve_id=String(target.curve_id),meson=String(meson),cut=String(cut),
            window=String(target.window),T_MeV=Float64(target.T_MeV),muB_MeV=Float64(target.muB_MeV),levels=levels))
        checkpoint()
    end
    points=NamedTuple[]
    for (target,location) in zip(config.targets,locations)
        center=location.levels["base"].xi
        offsets=Float64.(config.xi_offsets)
        for (index,offset) in enumerate(offsets)
            x=center+offset
            Float64(target.segment_left)<=x<=Float64(target.segment_right) || continue
            push!(points,(point_id=String(target.curve_id)*"_p"*lpad(string(index),2,'0'),
                curve_id=String(target.curve_id),meson=String(target.meson),cut=String(target.cut),window=String(target.window),
                T_MeV=Float64(target.T_MeV),muB_MeV=Float64(target.muB_MeV),xi=x,xi_mott=center,delta_xi=offset,
                kind="mott_centered",original_point_id=""))
        end
        for (index,id) in enumerate((String(target.point_left),String(target.point_right)))
            row=only(r for r in SOURCE_OTHER if String(r.point_id)==id)
            push!(points,(point_id=String(target.curve_id)*"_c"*string(index),curve_id=String(target.curve_id),
                meson=String(target.meson),cut=String(target.cut),window=String(target.window),
                T_MeV=Float64(target.T_MeV),muB_MeV=Float64(target.muB_MeV),xi=Float64(row.xi),xi_mott=center,
                delta_xi=Float64(row.xi)-center,kind="original_control",original_point_id=id))
        end
    end
    CSV.write(joinpath(OUTPUT_DIR,"transport_plan.csv"),points)
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",targets=length(locations),
        transport_plan_points=length(points),calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        production_kernel_changed=false,targets_sha256=digest(joinpath(INPUT_DIR,"targets.json")),
        input_rows_sha256=digest(joinpath(INPUT_DIR,"source_rows.csv")),script_sha256=digest(@__FILE__),
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

function transport_other()
    curve_id,variant=ARGS[5],ARGS[6]
    config=JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String))
    target=only(t for t in config.targets if String(t.curve_id)==curve_id)
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,variant*".toml"))
    runtime=SP_OTHER.build_runtime(profile)
    plan_path=joinpath(get(ENV,"MOTT_LOCATE_DIR",joinpath(dirname(OUTPUT_DIR),"locate_v1")),"transport_plan.csv")
    points=filter(r->String(r.curve_id)==curve_id,collect(CSV.File(plan_path)))
    wanted=get(ENV,"MOTT_POINT_IDS","all")
    wanted=="all" || (points=filter(r->String(r.point_id) in split(wanted,','),points))
    isempty(points) && error("no selected points")
    sort!(points;by=r->Float64(r.xi))
    original=only(r for r in SOURCE_OTHER if String(r.point_id)==String(target.point_left))
    seed=original_seed_other(runtime.model,original)
    observables,rates,weights,contributions=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    for point in points
        began=time();bg=background_other(target,Float64(point.xi),profile,"base";seed=seed)
        seed=bg.equilibrium.x_state
        info=(point_id=String(point.point_id),curve_id=curve_id,window=String(target.window),xi=Float64(point.xi))
        println("POINT_START ",point.point_id," xi=",point.xi);flush(stdout)
        on_channel=r->begin
            push!(rates,merge(info,r));CSV.write(joinpath(OUTPUT_DIR,"rates.csv"),rates)
        end
        result=SP_OTHER.solve_profile_transport(runtime,bg.equilibrium,bg.row.T_fm,bg.row.muq_fm;
            xi=bg.row.xi,compute_bulk=true,on_channel=on_channel)
        row=merge(bg.row,info,(meson=String(target.meson),cut=String(target.cut),variant=variant,
            xi_mott=Float64(point.xi_mott),delta_xi=Float64(point.delta_xi),point_kind=String(point.kind),
            original_point_id=ismissing(point.original_point_id) ? "" : String(point.original_point_id),
            entropy=result.thermo_background.entropy,s_fm3inv=result.thermo_background.entropy,
            eta=result.transport.eta,zeta=result.transport.zeta,sigma=result.transport.sigma,
            eta_over_s=result.transport.eta/result.thermo_background.entropy,
            zeta_over_s=result.transport.zeta/result.thermo_background.entropy,
            sigma_over_T=result.transport.sigma/bg.row.T_fm,seconds=time()-began),
            NamedTuple{Tuple(Symbol("tau_",s) for s in SPECIES_OTHER)}(Tuple(result.tau)),
            NamedTuple{Tuple(Symbol("tauinv_",s) for s in SPECIES_OTHER)}(Tuple(result.tau_inv)),
            NamedTuple{Tuple(Symbol("n_",s) for s in SPECIES_OTHER)}(Tuple(result.densities)))
        push!(observables,row)
        for species in SPECIES_OTHER
            unit_tau=NamedTuple{SPECIES_OTHER}(Tuple(s==species ? 1.0 : 0.0 for s in SPECIES_OTHER))
            eta=TC_OTHER.shear_viscosity(bg.state.qp,bg.state.tp;tau=unit_tau,config=runtime.transport_config)
            sigma=TC_OTHER.electric_conductivity(bg.state.qp,bg.state.tp;tau=unit_tau,config=runtime.transport_config)
            zeta=TC_OTHER.bulk_viscosity_isentropic(bg.state.qp,bg.state.tp;tau=unit_tau,
                config=runtime.transport_config,bulk_coeffs_isentropic=result.bulk_coeffs)
            push!(weights,merge(info,(species=String(species),tau_fm=getproperty(result.tau,species),
                eta_coefficient=eta,zeta_coefficient=zeta,sigma_coefficient=sigma,
                eta_over_s_coefficient=eta/row.entropy,zeta_over_s_coefficient=zeta/row.entropy,
                sigma_over_T_coefficient=sigma/bg.row.T_fm)))
        end
        for name in (:eta,:zeta,:sigma)
            reconstructed=sum(getproperty(w,Symbol(name,"_coefficient"))*w.tau_fm for w in weights if w.point_id==info.point_id)
            abs(reconstructed/getproperty(row,name)-1)<1e-10 || error("terminal species replay failed: $name")
        end
        rr=filter(r->r.point_id==info.point_id,rates)
        length(rr)==17 || error("incomplete transport channels")
        for r in rr
            single=Dict(Symbol(v.channel)=>(v.channel==r.channel ? Float64(r.rate) : 0.0) for v in rr)
            partial=RTM.relaxation_rates(result.densities,single)
            for species in SPECIES_OTHER
                rc=getproperty(partial,species)
                push!(contributions,merge(info,(species=String(species),channel=String(r.channel),
                    rate_fm2=Float64(r.rate),collision_density_weight_inv_fm3=rc/Float64(r.rate),
                    inverse_tau_contribution_inv_fm=rc,inverse_tau_inv_fm=getproperty(result.tau_inv,species))))
            end
        end
        for (file,rows) in (("observables",observables),("species_weights",weights),("channel_contributions",contributions))
            CSV.write(joinpath(OUTPUT_DIR,file*".csv"),rows)
        end
        println("POINT_DONE ",point.point_id," seconds=",time()-began);flush(stdout);GC.gc()
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",curve_id=curve_id,variant=variant,
        points=length(observables),channels=length(rates),species_weights=length(weights),
        point_selection=wanted,production_kernel_changed=false,
        effective_profile=SP_OTHER.effective_config(profile),background_actual=(residual_norm_max=1e-10,max_iter=80),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",script_sha256=digest(@__FILE__),
        plan_sha256=digest(plan_path),targets_sha256=digest(joinpath(INPUT_DIR,"targets.json")),
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==6 || error("SOURCE_ROOT INPUT_DIR OUTPUT_DIR STAGE SELECTION PROFILE required")
    ARGS[4]=="locate" ? locate_other() : ARGS[4]=="transport" ? transport_other() : error("unknown stage")
end
