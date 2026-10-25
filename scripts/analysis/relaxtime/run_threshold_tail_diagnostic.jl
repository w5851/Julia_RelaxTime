#!/usr/bin/env julia
# Resolve the finite-x feature in the retained xi=0.37 diagnostic plot.
include(joinpath(@__DIR__,"run_threshold_cache_candidates.jl"))
const TAIL_SA=Main.RelaxTime.ScatteringAmplitude
const TAIL_TP=Main.RelaxTime.TotalPropagator

function tail_components(process,s,state;n_points=24)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    ctx=TAIL_SA.prepare_scattering_context(process,qp)
    mi,mj,mc,md=ctx.m1,ctx.m2,ctx.m3,ctx.m4
    bounds=CS.calculate_t_bounds(s,mi,mj,mc,md)
    ctx.particle_c!=ctx.particle_d || error("Diagnostic assumes distinguishable final particles")
    z,w=CS.standard_nodes_weights(n_points)
    scale=(bounds.t_max-bounds.t_min)/2
    shift=(bounds.t_max+bounds.t_min)/2
    k0=sqrt(s)
    dsS,dsP=TAIL_TP.calculate_propagator_by_channel(process,:s,k0,0.0,qp,tp,kc)
    dsMixedS=TAIL_TP.total_propagator_mixed(process,:s,:S,k0,0.0,qp,tp,kc)
    Ec,Ed=CS.calculate_final_state_energies(s,mc,md)
    muc,mud=CS.get_chemical_potential(ctx.particle_c,qp),CS.get_chemical_potential(ctx.particle_d,qp)
    integrals=zeros(9)
    min_m2=Inf
    for (zi,wi) in zip(z,w)
        t=scale*zi+shift
        u=mi^2+mj^2+mc^2+md^2-s-t
        cms=TAIL_TP.calculate_cms_momentum(process,s,t,:t,qp;u=u)
        dtS,dtP=TAIL_TP.calculate_propagator_by_channel(process,:t,cms.k0,cms.k,qp,tp,kc)
        v=TAIL_SA.calculate_mandelstam_variables(s,t,u,mi,mj,mc,md)
        a=v.s_12_plus*v.s_34_plus-v.u_14_plus*v.u_23_plus+v.t_13_plus*v.t_24_plus
        b=v.s_12_plus*v.s_34_plus-v.u_14_minus*v.u_23_minus+v.t_13_minus*v.t_24_minus
        c=v.s_12_minus*v.s_34_minus-v.u_14_minus*v.u_23_minus+v.t_13_plus*v.t_24_plus
        d=v.s_12_minus*v.s_34_minus-v.u_14_plus*v.u_23_plus+v.t_13_minus*v.t_24_minus
        scalar_term(D)=abs2(D)*v.s_12_plus*v.s_34_plus
        scalar_cross(D)=-2TAIL_SA.CROSS_TERM_FACTOR*real(D*(conj(dtS)*a-conj(dtP)*b))
        ss=scalar_term(dsS)
        sp=abs2(dsP)*v.s_12_minus*v.s_34_minus
        ts=abs2(dtS)*v.t_13_plus*v.t_24_plus
        tps=abs2(dtP)*v.t_13_minus*v.t_24_minus
        xs=scalar_cross(dsS)
        xp=-2TAIL_SA.CROSS_TERM_FACTOR*real(dsP*(-conj(dtS)*c+conj(dtP)*d))
        full=ss+sp+ts+tps+xs+xp
        min_m2=min(min_m2,full)
        no_s_scalar=sp+ts+tps+xp
        no_s_mixed_scalar=no_s_scalar+scalar_term(dsS-dsMixedS)+scalar_cross(dsS-dsMixedS)
        ct=CS.cos_theta_star(s,t,mi,mj,mc,md)
        block=CS.combined_final_state_factor(ctx.particle_c,ctx.particle_d,Ec,Ed,mc,md,
            muc,mud,tp.T,tp.Φ,tp.Φbar,tp.ξ,ct)
        factor=scale*wi*block*CS.differential_cross_section(v.s_12_plus,v.s_12_minus,1.0)
        integrals .+= factor .* (ss,sp,ts,tps,xs,xp,max(full,0.0),max(no_s_scalar,0.0),max(no_s_mixed_scalar,0.0))
    end
    return (s_scalar=integrals[1],s_pseudoscalar=integrals[2],t_scalar=integrals[3],
        t_pseudoscalar=integrals[4],interference_scalar=integrals[5],
        interference_pseudoscalar=integrals[6],reconstructed=integrals[7],
        no_s_scalar=integrals[8],no_s_mixed_scalar=integrals[9],min_amplitude_squared=min_m2,
        s_scalar_propagator_abs2=abs2(dsS),s_mixed_scalar_abs2=abs2(dsMixedS))
end

function tail_denominators(state,s,plans)
    p=paired_polarizations(state,s,plans)
    kc=state.K_coeffs
    ms=Main.RelaxTime.EffectiveCouplings.mixing_matrix_elements(p.u.S,p.s.S,kc,:S)
    mp=Main.RelaxTime.EffectiveCouplings.mixing_matrix_elements(p.u.P,p.s.P,kc,:P)
    return (scalar_mixed=(ms.M00*ms.M88-ms.M08^2)/kc.det_K_minus,
            pseudoscalar_mixed=(mp.M00*mp.M88-mp.M08^2)/kc.det_K_plus,
            scalar_triplet=1-4kc.K123_minus*p.u.S,
            pi=1-4kc.K123_plus*p.u.P)
end

function run_tail_diagnostic(ablationdir,outdir)
    start_run(outdir,joinpath(ablationdir,"inputs"),@__FILE__,(
        purpose="finite-x peak: direct quadrature, loop sensitivity and amplitude decomposition",
        xi=0.37,channel="uubar_to_uubar",x_range=[0.12,0.34],x_step=0.001,
        t_orders=[24,48],paired_orders=[64,128],
        counterfactual="remove s scalar or only s mixed scalar including its interference; diagnostic only"),
        extra_sources=[joinpath(@__DIR__,"ThresholdCacheCandidates.jl"),
                       joinpath(@__DIR__,"run_threshold_cache_candidates.jl")])
    case=load_upgrade_case(ablationdir,0.37,:uubar_to_uubar)
    state=case.state; qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    plans(n)=(u=PairLoopPlan(qp.m.u,qp.μ.u,tp;n_nodes=n),s=PairLoopPlan(qp.m.s,qp.μ.s,tp;n_nodes=n))
    p64,p128=plans(64),plans(128)
    result=NamedTuple[]
    for (i,x) in enumerate(range(0.12,0.34;step=0.001))
        s=case.s0+x
        old=CS.total_cross_section(case.process,s,qp,tp,kc;n_points=24)
        old48=CS.total_cross_section(case.process,s,qp,tp,kc;n_points=48)
        new64=paired_cross_section(case.process,s,state,p64;n_points=24)
        new128=paired_cross_section(case.process,s,state,p128;n_points=24)
        parts=tail_components(case.process,s,state)
        isapprox(parts.reconstructed,old;rtol=1e-10,atol=1e-12) ||
            error("Amplitude decomposition does not reconstruct the production kernel at x=$x")
        den=tail_denominators(state,s,p128)
        push!(result,merge((x=x,s=s,energy_inv_fm=sqrt(s),sigma=old,
            sigma_t48=old48,sigma_paired64=new64,sigma_paired128=new128,
            t_relative_change=old48/old-1,paired_relative_change=new64/old-1,
            paired_order_relative_change=new128/new64-1,
            scalar_mixed_re=real(den.scalar_mixed),scalar_mixed_im=imag(den.scalar_mixed),
            scalar_triplet_re=real(den.scalar_triplet),scalar_triplet_im=imag(den.scalar_triplet),
            pseudoscalar_mixed_re=real(den.pseudoscalar_mixed),pseudoscalar_mixed_im=imag(den.pseudoscalar_mixed),
            pi_re=real(den.pi),pi_im=imag(den.pi)),parts))
        if i%20==0
            write_rows(joinpath(outdir,"tail_feature.csv"),result)
            println("TAIL ",i,"/221");flush(stdout)
        end
    end
    write_rows(joinpath(outdir,"tail_feature.csv"),result)
    # Refine sign-changing real denominators. This is a real-axis resonance
    # diagnostic, not a complex meson pole or a new Mott threshold calculation.
    roots=NamedTuple[]
    for field in (:scalar_mixed,:scalar_triplet,:pseudoscalar_mixed,:pi)
        re_field=Symbol(field,"_re")
        for i in 1:length(result)-1
            a,b=result[i].x,result[i+1].x
            fa,fb=getproperty(result[i],re_field),getproperty(result[i+1],re_field)
            fa*fb<0 || continue
            for _ in 1:42
                m=(a+b)/2
                fm=real(getproperty(tail_denominators(state,case.s0+m,p128),field))
                if signbit(fm)==signbit(fa)
                    a=m;fa=fm
                else
                    b=m
                end
            end
            x=(a+b)/2
            d=getproperty(tail_denominators(state,case.s0+x,p128),field)
            push!(roots,(family=String(field),x=x,energy_inv_fm=sqrt(case.s0+x),
                         real_denominator=real(d),imag_denominator=imag(d)))
        end
    end
    write_rows(joinpath(outdir,"denominator_roots.csv"),roots)
    finish_run(outdir;points=length(result),max_t_relative_change=maximum(abs(r.t_relative_change) for r in result),
        max_paired_relative_change=maximum(abs(r.paired_relative_change) for r in result),
        max_paired_order_relative_change=maximum(abs(r.paired_order_relative_change) for r in result),
        maximum_reconstruction_relative_error=maximum(abs(r.reconstructed/r.sigma-1) for r in result))
    println("COMPLETE ",outdir)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: run_threshold_tail_diagnostic.jl ABLATION_DIR NEW_OUTPUT_DIR")
    run_tail_diagnostic(abspath(ARGS[1]),abspath(ARGS[2]))
end
