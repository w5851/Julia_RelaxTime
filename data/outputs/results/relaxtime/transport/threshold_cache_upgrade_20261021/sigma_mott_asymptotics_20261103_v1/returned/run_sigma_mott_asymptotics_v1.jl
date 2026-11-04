#!/usr/bin/env julia
# Isolated quote-remote diagnostic. The frozen production source is never edited.
include(joinpath(@__DIR__, "run_sigma_mott_thermal_audit.jl"))

const SIG_PROP = Main.RelaxTime.MesonPropagator

function asym_background(target, xi, profile, seed; level="base")
    p, c, ap, ac, limit = level == "base" ? (96,24,32,16,25.0) : (128,32,64,32,30.0)
    model = SP_OTHER.profile_model(profile; thermal_p_max_inv_fm=limit)
    raw = Models.solve_constraint(model, Models.FixedMu(), Float64(target.T_MeV)/HC;
        μ_fm=Float64(target.muB_MeV)/(3HC), seed_guess=Vector(seed),
        seed_candidates=[Vector(seed)], xi=xi, p_num=p, t_num=c,
        solver=Models.NLsolveGapSolver(method=:newton,jacobian=:forward,xtol=1e-14,ftol=2e-14),
        residual_norm_max=1e-12, iterations=100, physicality_check=Models.is_physical_solution)
    raw.converged || error("Strict local background failed at xi=$xi")
    eq = PHASE_OTHER._normalize_equilibrium_result(raw;solver_backend=:models_strict_local)
    eq === nothing && error("Invalid strict local background")
    row=(point_id="asym",window=String(target.window),T_MeV=Float64(target.T_MeV),
        muB_MeV=Float64(target.muB_MeV),xi=xi,T_fm=Float64(target.T_MeV)/HC,
        muq_fm=Float64(target.muB_MeV)/(3HC),m_u=eq.masses[1],m_d=eq.masses[2],m_s=eq.masses[3],
        Phi=eq.x_state[4],Phibar=eq.x_state[5],phi_u=eq.x_state[1],phi_d=eq.x_state[2],
        phi_s=eq.x_state[3],background_residual=eq.residual_norm,omega=eq.omega)
    (row=row,equilibrium=eq,state=frozen_state(row,profile,ap,ac))
end

asym_plans(state,n)=(u=PAIR_OTHER.PairLoopPlan(state.qp.m.u,state.qp.μ.u,state.tp;n_nodes=n),
                    s=PAIR_OTHER.PairLoopPlan(state.qp.m.s,state.qp.μ.s,state.tp;n_nodes=n))

function asym_matrix(state,x,plans)
    s=4state.qp.m.u^2+x
    pol=map((:u,:s)) do f
        mass=getproperty(state.qp.m,f)
        coefficient=f==:u ? -x : 4mass^2-s
        loop=coefficient==0 ? complex(0.0) : coefficient*PAIR_OTHER.pair_b0_pv(getproperty(plans,f),s)
        -Main.Constants_PNJL.N_color/(8pi^2)*(2getproperty(state.qp.A,f)+loop)
    end
    EC.mixing_matrix_elements(pol...,state.kc,:S)
end

function asym_mixed(state,x,plans;regularized=false)
    e=asym_matrix(state,x,plans)
    det=e.M00*e.M88-e.M08^2
    regularized && abs(det)<SIG_PROP.PROPAGATOR_DENOM_EPS &&
        (det+=complex(0.0,SIG_PROP.PROPAGATOR_DENOM_EPS))
    # uubar and ddbar have the same diagonal lambda_0/lambda_8 currents.
    j0,j8=SIG_PROP._current_components(:u,:ubar)
    numerator=j0*(e.M00*j0+e.M08*j8)+j8*(e.M08*j0+e.M88*j8)
    2e.det_K*numerator/det
end

function asym_scalar_sigma(state,x,plans;regularized=false,n=24)
    s=4state.qp.m.u^2+x
    d=asym_mixed(state,x,plans;regularized=regularized)
    z,w=SIG_CS.standard_nodes_weights(n)
    m,tp=state.qp.m.u,state.tp
    en=sqrt(s)/2;mu=state.qp.μ.u
    block=sum(wi*SIG_CS.combined_final_state_factor(:u,:ubar,en,en,m,m,mu,mu,
        tp.T,tp.Φ,tp.Φbar,tp.ξ,zi) for (zi,wi) in zip(z,w))/2
    abs2(d)*x^2/(16pi*s)*block
end

function asym_root(state,plans)
    f=x->real(MP._mixed_inverse(:S,:light,asym_matrix(state,x,plans)))
    f0=f(0.0); b=(f(1e-5)-f0)/1e-5
    guess=-f0/b
    if guess<=0
        return (x=guess,width=0.0,det_abs=abs(let e=asym_matrix(state,0.0,plans);e.M00*e.M88-e.M08^2 end),
                inverse_zero=f0,bracket=0.0)
    end
    lo,hi=0.0,max(2guess,1e-9)
    flo=f0
    f(lo)*f(hi)<=0 || error("Root is not bracketed")
    for _ in 1:65
        mid=(lo+hi)/2
        lo<mid<hi || break
        fm=f(mid)
        if signbit(fm)==signbit(flo)
            lo=mid;flo=fm
        else
            hi=mid
        end
        hi-lo<1e-16 && break
    end
    x=(lo+hi)/2
    dx=max(x*1e-3,1e-9)
    slope=(f(x+dx)-f(x-dx))/(2dx)
    e=asym_matrix(state,x,plans)
    inv=MP._mixed_inverse(:S,:light,e)
    (x=x,width=abs(imag(inv)/slope),det_abs=abs(e.M00*e.M88-e.M08^2),
        inverse_zero=f0,bracket=hi-lo)
end

function asym_gl(f,a,b,n)
    z,w=SIG_CS.standard_nodes_weights(n)
    sum(((b-a)/2*wi).*f((a+b)/2+(b-a)/2*zi) for (zi,wi) in zip(z,w))
end

function asym_integrate(f,a,b;atol=1e-21,rtol=2e-11,depth=0)
    coarse=asym_gl(f,a,b,16)
    fine=asym_gl(f,a,b,32)
    err=abs.(fine.-coarse)
    if all(err .<= atol .+ rtol.*abs.(fine))
        return (value=fine,error=err,leaves=1)
    end
    depth<18 || error("Unresolved continuous integral [$a,$b], defect=$err")
    left=asym_integrate(f,a,(a+b)/2;atol=atol/2,rtol=rtol,depth=depth+1)
    right=asym_integrate(f,(a+b)/2,b;atol=atol/2,rtol=rtol,depth=depth+1)
    (value=left.value.+right.value,error=left.error.+right.error,leaves=left.leaves+right.leaves)
end

function asym_integral(f,root,xmax;atol=1e-20,rtol=1e-9)
    value=zeros(2);err=zeros(2);leaves=0
    # A finite central peak map and geometric tails retain the full [0,xmax] integral.
    # No assumed h exponent or fitted singular term enters this quadrature.
    if root.x>0 && root.x<xmax
        r,w=root.x,root.width
        w>0 || error("Nonpositive resonance width")
        a,b=max(r/2,r-64w),min(xmax,r+64w)
        tlo,thi=atan((a-r)/w),atan((b-r)/w)
        mapped=t->begin
            x=r+w*tan(t)
            f(x).*(w/cos(t)^2)
        end
        result=asym_integrate(mapped,tlo,thi;atol=atol/4,rtol=rtol)
        value .+= result.value;err .+= result.error;leaves+=result.leaves
        left=sort!(unique(vcat(0.,a,[r-w*2.0^j for j in 6:50 if 0<r-w*2.0^j<a])))
        right=sort!(unique(vcat(b,xmax,[r+w*2.0^j for j in 6:50 if b<r+w*2.0^j<xmax])))
        intervals=vcat(collect(zip(left[1:end-1],left[2:end])),collect(zip(right[1:end-1],right[2:end])))
    else
        scale=max(abs(root.x),1e-12)
        edges=sort!(unique(vcat(0.,xmax,[scale*2.0^j for j in -3:50 if 0<scale*2.0^j<xmax])))
        intervals=collect(zip(edges[1:end-1],edges[2:end]))
    end
    for (a,b) in intervals
        result=asym_integrate(f,a,b;atol=atol/(2length(intervals)),rtol=rtol)
        value .+= result.value;err .+= result.error;leaves+=result.leaves
    end
    (value=value,error=err,leaves=leaves)
end

function run_asymptotics()
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"s6.toml"))
    runtime=SP_OTHER.build_runtime(profile)
    target=only(t for t in JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String)).targets if String(t.curve_id)=="other07")
    oldcenter=0.015327819995582102
    original=only(r for r in SOURCE_OTHER if String(r.point_id)==String(target.point_left))
    seed=original_seed_other(runtime.model,original)
    initial=background_other(target,oldcenter,profile,"base";seed=seed)
    seed=initial.equilibrium.x_state
    bgcache=Dict{Float64,Any}()
    bg(x)=get!(bgcache,x) do
        asym_background(target,x,profile,seed)
    end
    endpoint(x)=real(MP._mixed_inverse(:S,:light,asym_matrix(bg(x).state,0.0,asym_plans(bg(x).state,128))))
    lo,hi=oldcenter-2e-6,oldcenter+2e-6;flo=endpoint(lo)
    flo*endpoint(hi)<0 || error("Missing strict Mott center bracket")
    for _ in 1:40
        mid=(lo+hi)/2;fm=endpoint(mid)
        if signbit(fm)==signbit(flo)
            lo=mid;flo=fm
        else
            hi=mid
        end
        hi-lo<3e-13 && break
    end
    center=(lo+hi)/2
    centerbg=bg(center)
    endpoint_slope=(endpoint(center+1e-6)-endpoint(center-1e-6))/2e-6
    write_json(joinpath(OUTPUT_DIR,"center.json"),(xi=center,left=lo,right=hi,old_xi=oldcenter,
        shift=center-oldcenter,inverse=endpoint(center),inverse_slope=endpoint_slope,
        background_residual=centerbg.row.background_residual))
    println("ASYM_CENTER ",center," shift=",center-oldcenter);flush(stdout)
    mode=get(ENV,"SIGMA_ASYM_MODE","pilot")
    steps=mode=="pilot" ? [2e-4,2e-5,2e-6,2e-7] : [2e-4/2.0^j for j in 0:10]
    offsets=mode=="pilot" ? sort!(unique(vcat(-steps,0.,steps))) :
        sort!(unique([sign*h*factor for sign in (-1.,1.),h in steps,factor in (.75,.875,1.,1.125,1.25)]))
    rows,checks=NamedTuple[],NamedTuple[]
    xmax=0.01
    for (index,u) in enumerate(offsets)
        started=time();b=bg(center+u);state=b.state
        plans=asym_plans(state,128);root=asym_root(state,plans)
        problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,state.qp,state.tp,state.kc;runtime.tau_kwargs...)
        prefactor=problem.plan.prefactor
        g=sigma_weight_table(state,prefactor,xmax;order=8)
        evals=Ref(0)
        f=x->begin
            evals[]+=1
            a=asym_scalar_sigma(state,x,plans)
            [x*g(x)*a,0.0]
        end
        result=asym_integral(f,root,xmax;atol=1e-20,rtol=2e-9)
        row=(point_id="asym_"*lpad(index,3,'0'),xi=center+u,delta_xi=u,
            mixed_square_rate_fm2=result.value[1],quadrature_error_fm2=result.error[1],
            root_x_inv_fm2=root.x,half_width_inv_fm2=root.width,root_det_abs=root.det_abs,
            guard_active_at_peak=root.det_abs<SIG_PROP.PROPAGATOR_DENOM_EPS,
            root_bracket_inv_fm2=root.bracket,inverse_at_threshold=root.inverse_zero,
            background_residual=b.row.background_residual,leaves=result.leaves,evaluations=evals[],seconds=time()-started)
        push!(rows,row)
        if mode=="pilot"
            for x in unique(vcat([1e-5,1e-3],root.x>0 ? [root.x] : Float64[]))
                SIG_PC.reset_cache!()
                native=sigma_components(:uubar_to_uubar,4state.qp.m.u^2+x,state).s_mixed_scalar_square_fm2
                alternate=asym_scalar_sigma(state,x,plans;regularized=true)
                raw=asym_scalar_sigma(state,x,plans)
                push!(checks,(delta_xi=u,x_inv_fm2=x,native_mixed_sigma_fm2=native,
                    pair_regularized_sigma_fm2=alternate,pair_unfloored_sigma_fm2=raw,
                    relative_pair_native=alternate/native-1,regularized_over_unfloored=alternate/raw))
            end
        end
        CSV.write(joinpath(OUTPUT_DIR,"integrals.csv"),rows)
        isempty(checks) || CSV.write(joinpath(OUTPUT_DIR,"kernel_checks.csv"),checks)
        println("ASYM_POINT ",index,"/",length(offsets)," u=",u," I=",result.value[1],
            " err=",result.error[1]," guard=",row.guard_active_at_peak," seconds=",round(time()-started;digits=2));flush(stdout)
        SIG_PC.reset_cache!();GC.gc()
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",mode=mode,points=length(rows),
        xi_mott=center,steps=steps,xmax_inv_fm2=xmax,calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        scope="Continuum near-threshold mixed-scalar square; both mixed poles; no other-amplitude interference or full transport reconstruction.",
        numerical_representation="Pair PV quadrature, exact x^2 numerator and continuous weight; determinant guard omitted only in this isolated diagnostic.",
        production_kernel_changed=false,files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_asymptotics()
end
