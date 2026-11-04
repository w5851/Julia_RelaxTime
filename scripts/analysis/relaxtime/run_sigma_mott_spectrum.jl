#!/usr/bin/env julia
# quote-remote only. Uses the unchanged frozen kernel and explicit native roots.
const SIG_LIBRARY_PATH = let archived=joinpath(@__DIR__, "run_transport_other_mott_refinement_library_v2.jl")
    isfile(archived) ? archived : joinpath(@__DIR__, "run_transport_other_mott_refinement.jl")
end
include(SIG_LIBRARY_PATH)
const SIG_LIBRARY_SHA = digest(SIG_LIBRARY_PATH)
const SIG_CS = Main.RelaxTime.TotalCrossSection
const SIG_SA = Main.RelaxTime.ScatteringAmplitude
const SIG_TP = Main.RelaxTime.TotalPropagator
const SIG_PC = Main.RelaxTime.PolarizationCache
const SIG_TS = AS.ThermalWeightedSampling

function sigma_inverse_pair(state, x, plans)
    s = 4state.qp.m.u^2 + x
    pol = map((:u, :s)) do flavor
        mass = getproperty(state.qp.m, flavor)
        coefficient = flavor == :u ? -x : 4mass^2 - s
        loop = coefficient == 0 ? complex(0.0) : coefficient * PAIR_OTHER.pair_b0_pv(getproperty(plans, flavor), s)
        -Main.Constants_PNJL.N_color / (8pi^2) * (2getproperty(state.qp.A, flavor) + loop)
    end
    MP._mixed_inverse(:S, :light, EC.mixing_matrix_elements(pol..., state.kc, :S))
end

sigma_inverse_native(state, x) = MP.meson_mass_equation(:sigma, sqrt(4state.qp.m.u^2+x),
    0.0, 0.0, state.qp, state.tp, state.kc)

function sigma_real_roots(f, grid, s0)
    roots = NamedTuple[]
    for (left, right) in zip(grid[1:end-1], grid[2:end])
        fa, fb = real(f(left)), real(f(right))
        isfinite(fa) && isfinite(fb) || error("Nonfinite sigma root screen")
        signbit(fa) == signbit(fb) && fa != 0 && fb != 0 && continue
        a, b = left, right
        for _ in 1:60
            m = (a+b)/2
            a < m < b || break
            fm = real(f(m))
            if signbit(fm) == signbit(fa)
                a=m; fa=fm
            else
                b=m
            end
            b-a <= 16eps(s0) && break
        end
        x = (a+b)/2
        h = max(256eps(s0), min(x*1e-3, 1e-6))
        lo, hi = max(x-h, x/2), x+h
        slope = (real(f(hi))-real(f(lo)))/(hi-lo)
        value = f(x)
        width = abs(imag(value)/slope)
        all(isfinite, (x, width, slope, real(value), imag(value))) || error("Invalid local sigma root")
        push!(roots, (x_inv_fm2=x, half_width_inv_fm2=width, slope=slope,
            inverse_re=real(value), inverse_im=imag(value), bracket_width_inv_fm2=b-a,
            width_resolved=width>256eps(s0)))
    end
    roots
end

"""Separate the two scalar s-channel terms, retaining their interference.

The mixed term is the complete sigma/sigma-prime flavor propagator. Removing it
is a diagnostic counterfactual, not a positive single-meson rate decomposition.
"""
function sigma_components(process, s, state; n_points=24)
    qp,tp,kc=state.qp,state.tp,state.kc
    ctx=SIG_SA.prepare_scattering_context(process,qp)
    mi,mj,mc,md=ctx.m1,ctx.m2,ctx.m3,ctx.m4
    all(m->isapprox(m,mi;rtol=1e-12,atol=1e-14),(mj,mc,md)) || error("Equal light-mass diagnostic only")
    ctx.particle_c != ctx.particle_d || error("Distinguishable final particles required")
    bounds=SIG_CS.calculate_t_bounds(s,mi,mj,mc,md)
    z,w=SIG_CS.standard_nodes_weights(n_points)
    scale=(bounds.t_max-bounds.t_min)/2
    shift=(bounds.t_max+bounds.t_min)/2
    dsS,dsP=SIG_TP.calculate_propagator_by_channel(process,:s,sqrt(s),0.0,qp,tp,kc)
    mixed=SIG_TP.total_propagator_mixed(process,:s,:S,sqrt(s),0.0,qp,tp,kc)
    Ec,Ed=SIG_CS.calculate_final_state_energies(s,mc,md)
    muc,mud=SIG_CS.get_chemical_potential(ctx.particle_c,qp),SIG_CS.get_chemical_potential(ctx.particle_d,qp)
    totals=zeros(10)
    for (zi,wi) in zip(z,w)
        t=scale*zi+shift
        u=mi^2+mj^2+mc^2+md^2-s-t
        cms=SIG_TP.calculate_cms_momentum(process,s,t,:t,qp;u=u)
        dtS,dtP=SIG_TP.calculate_propagator_by_channel(process,:t,cms.k0,cms.k,qp,tp,kc)
        v=SIG_SA.calculate_mandelstam_variables(s,t,u,mi,mj,mc,md)
        a=v.s_12_plus*v.s_34_plus-v.u_14_plus*v.u_23_plus+v.t_13_plus*v.t_24_plus
        b=v.s_12_plus*v.s_34_plus-v.u_14_minus*v.u_23_minus+v.t_13_minus*v.t_24_minus
        c=v.s_12_minus*v.s_34_minus-v.u_14_minus*v.u_23_minus+v.t_13_plus*v.t_24_plus
        d=v.s_12_minus*v.s_34_minus-v.u_14_plus*v.u_23_plus+v.t_13_minus*v.t_24_minus
        scalar_term(D)=abs2(D)*v.s_12_plus*v.s_34_plus
        scalar_cross(D)=-2SIG_SA.CROSS_TERM_FACTOR*real(D*(conj(dtS)*a-conj(dtP)*b))
        ss=scalar_term(dsS)
        sp=abs2(dsP)*v.s_12_minus*v.s_34_minus
        ts=abs2(dtS)*v.t_13_plus*v.t_24_plus
        tps=abs2(dtP)*v.t_13_minus*v.t_24_minus
        xs=scalar_cross(dsS)
        xp=-2SIG_SA.CROSS_TERM_FACTOR*real(dsP*(-conj(dtS)*c+conj(dtP)*d))
        full=max(ss+sp+ts+tps+xs+xp,0.0)
        no_scalar=max(sp+ts+tps+xp,0.0)
        no_mixed=max(sp+ts+tps+xp+scalar_term(dsS-mixed)+scalar_cross(dsS-mixed),0.0)
        ct=SIG_CS.cos_theta_star(s,t,mi,mj,mc,md)
        block=SIG_CS.combined_final_state_factor(ctx.particle_c,ctx.particle_d,Ec,Ed,mc,md,
            muc,mud,tp.T,tp.Φ,tp.Φbar,tp.ξ,ct)
        factor=scale*wi*block*SIG_CS.differential_cross_section(v.s_12_plus,v.s_12_minus,1.0)
        totals .+= factor .* (full,ss,sp,ts,tps,xs,xp,no_scalar,no_mixed,scalar_term(mixed))
    end
    (sigma_fm2=totals[1],s_scalar_fm2=totals[2],s_pseudoscalar_fm2=totals[3],
        t_scalar_fm2=totals[4],t_pseudoscalar_fm2=totals[5],interference_scalar_fm2=totals[6],
        interference_pseudoscalar_fm2=totals[7],no_s_scalar_fm2=totals[8],no_s_mixed_scalar_fm2=totals[9],
        s_mixed_scalar_square_fm2=totals[10],mixed_delta_fm2=totals[1]-totals[9],
        mixed_propagator_abs2_fm4=abs2(mixed))
end

function sigma_peak_area(process,state,root,order)
    xr,width=root.x_inv_fm2,root.half_width_inv_fm2
    root.width_resolved || return nothing
    edges=sort!(unique(max.(xr/1000,xr .+ width .* [-32.,-16.,-8.,-4.,-2.,-1.,-.5,0.,.5,1.,2.,4.,8.,16.,32.])))
    z,w=SIG_CS.standard_nodes_weights(order)
    totals=zeros(4)
    s0=4state.qp.m.u^2
    for (a,b) in zip(edges[1:end-1],edges[2:end]),(zi,wi) in zip(z,w)
        x=(a+b)/2+(b-a)/2*zi
        SIG_PC.reset_cache!()
        part=sigma_components(process,s0+x,state)
        totals .+= (b-a)/2*wi .* (part.sigma_fm2,part.s_mixed_scalar_square_fm2,
                                 part.mixed_delta_fm2,x*part.s_mixed_scalar_square_fm2)
    end
    (order=order,x_left_inv_fm2=first(edges),x_right_inv_fm2=last(edges),
        full_area_dimensionless=totals[1],mixed_square_area_dimensionless=totals[2],
        mixed_delta_area_dimensionless=totals[3],x_mixed_square_area_inv_fm2=totals[4])
end

function run_sigma_spectrum()
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,ARGS[6]*".toml"))
    runtime=SP_OTHER.build_runtime(profile)
    config=JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String))
    target=only(t for t in config.targets if String(t.curve_id)=="other07")
    original=only(r for r in SOURCE_OTHER if String(r.point_id)==String(target.point_left))
    seed=original_seed_other(runtime.model,original)
    points=collect(CSV.File(joinpath(INPUT_DIR,"transport_plan.csv")))
    selected=get(ENV,"MOTT_POINT_IDS","all")
    selected=="all" || filter!(r->String(r.point_id) in split(selected,','),points)
    isempty(points) && error("No spectrum points")
    backgrounds,inverses,roots,spectra,areas,thermal=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    max_reconstruction=0.0
    function checkpoint()
        for (name,rows) in (("background",backgrounds),("inverse_spectrum",inverses),("real_axis_roots",roots),
                           ("cross_section_spectrum",spectra),("peak_areas",areas),("thermal_support",thermal))
            isempty(rows) || CSV.write(joinpath(OUTPUT_DIR,name*".csv"),rows)
        end
    end
    for point in points
        started=time()
        bg=background_other(target,Float64(point.xi),profile,"base";seed=seed)
        seed=bg.equilibrium.x_state
        state=bg.state
        info=(point_id=String(point.point_id),xi=bg.row.xi,delta_xi=Float64(point.delta_xi))
        push!(backgrounds,merge(bg.row,info))
        s0=4state.qp.m.u^2
        pair128=(u=PAIR_OTHER.PairLoopPlan(state.qp.m.u,state.qp.μ.u,state.tp;n_nodes=128),
                 s=PAIR_OTHER.PairLoopPlan(state.qp.m.s,state.qp.μ.s,state.tp;n_nodes=128))
        pair256=(u=PAIR_OTHER.PairLoopPlan(state.qp.m.u,state.qp.μ.u,state.tp;n_nodes=256),
                 s=PAIR_OTHER.PairLoopPlan(state.qp.m.s,state.qp.μ.s,state.tp;n_nodes=256))
        grid=10.0.^range(-10,log10(0.05);length=241)
        native=x->sigma_inverse_native(state,x)
        p128=x->sigma_inverse_pair(state,x,pair128)
        p256=x->sigma_inverse_pair(state,x,pair256)
        found=Dict{String,Any}()
        for (tag,f) in (("native",native),("pair128",p128),("pair256",p256))
            found[tag]=sigma_real_roots(f,grid,s0)
            for r in found[tag]
                push!(roots,merge(info,(backend=tag,),r))
            end
        end
        for x in grid
            n,a,b=native(x),p128(x),p256(x)
            push!(inverses,merge(info,(x_inv_fm2=x,s_inv_fm2=s0+x,native_re=real(n),native_im=imag(n),
                pair128_re=real(a),pair128_im=imag(a),pair256_re=real(b),pair256_im=imag(b))))
        end
        detailed=collect(grid)
        for r in found["native"]
            width=max(r.half_width_inv_fm2,256eps(s0))
            append!(detailed,[r.x_inv_fm2+y*width for y in range(-32.,32.;length=257)
                              if 1e-10<r.x_inv_fm2+y*width<0.05])
        end
        sort!(unique!(detailed))
        problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,state.qp,state.tp,state.kc;runtime.tau_kwargs...)
        measure=SIG_TS.thermal_measure(problem.plan,problem.context.bounds;relative_velocity=:cm)
        for r in found["native"]
            for multiple in (1.,4.,32.)
                lo=max(0.,r.x_inv_fm2-multiple*r.half_width_inv_fm2)
                hi=r.x_inv_fm2+multiple*r.half_width_inv_fm2
                ia=searchsortedfirst(measure.nodes,s0+lo)
                ib=searchsortedlast(measure.nodes,s0+hi)
                push!(thermal,merge(info,(kind="native_peak_band",width_multiple=multiple,x_left_inv_fm2=lo,
                    x_right_inv_fm2=hi,events=max(0,ib-ia+1),weight=sum(@view measure.weights[ia:ib]),
                    minimum_x_inv_fm2=first(measure.nodes)-s0,total_unique_events=length(measure.nodes))))
            end
        end
        for (lo,hi) in zip((0.,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2),
                           (1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.05))
            ia=searchsortedfirst(measure.nodes,s0+lo)
            ib=searchsortedlast(measure.nodes,s0+hi)
            push!(thermal,merge(info,(kind="fixed_x_band",width_multiple=0.,x_left_inv_fm2=lo,x_right_inv_fm2=hi,
                events=max(0,ib-ia+1),weight=sum(@view measure.weights[ia:ib]),
                minimum_x_inv_fm2=first(measure.nodes)-s0,total_unique_events=length(measure.nodes))))
        end
        for process in (:uubar_to_uubar,:uubar_to_ddbar)
            for (i,x) in enumerate(detailed)
                SIG_PC.reset_cache!()
                actual=SIG_CS.total_cross_section(process,s0+x,state.qp,state.tp,state.kc;n_points=24)
                part=sigma_components(process,s0+x,state)
                discrepancy=abs(actual-part.sigma_fm2)/max(abs(actual),1e-20)
                max_reconstruction=max(max_reconstruction,discrepancy)
                discrepancy<1e-8 || error("Scalar decomposition fails at $(point.point_id), $process, x=$x")
                push!(spectra,merge(info,(channel=String(process),x_inv_fm2=x,s_inv_fm2=s0+x,
                    native_sigma_fm2=actual,reconstruction_relative=discrepancy),part))
            end
            for r in found["native"]
                r.width_resolved || continue
                for order in (16,32)
                    push!(areas,merge(info,(channel=String(process),root_x_inv_fm2=r.x_inv_fm2,
                        half_width_inv_fm2=r.half_width_inv_fm2),sigma_peak_area(process,state,r,order)))
                end
            end
        end
        checkpoint()
        println("SPECTRUM_DONE ",point.point_id," native_roots=",length(found["native"]),
            " thermal_min_x=",first(measure.nodes)-s0," seconds=",round(time()-started;digits=2));flush(stdout)
        SIG_PC.reset_cache!();GC.gc()
    end
    digest(SIG_LIBRARY_PATH)==SIG_LIBRARY_SHA || error("Diagnostic library changed during calculation")
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",points=length(points),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",production_kernel_changed=false,
        max_reconstruction_relative=max_reconstruction,
        interpretation="Real-axis roots are resonance anchors, not complex meson poles. Peak areas integrate dx and are not thermal rates. Mixed scalar removal includes sigma and sigma-prime and their interference.",
        profile=SP_OTHER.effective_config(profile),script_sha256=digest(@__FILE__),library_sha256=SIG_LIBRARY_SHA,
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_sigma_spectrum()
end
