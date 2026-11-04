#!/usr/bin/env julia
# quote-remote verification of side curvature, kernel representation and regulator.
include(joinpath(@__DIR__, "run_sigma_mott_asymptotics_v3.jl"))

function check_asymptotics()
    source=joinpath(dirname(OUTPUT_DIR),"asym_full_v2")
    baseline=collect(CSV.File(joinpath(source,"integrals.csv")))
    meta=JSON3.read(read(joinpath(source,"summary.json"),String))
    center=Float64(meta.xi_mott)
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"s6.toml"))
    runtime=SP_OTHER.build_runtime(profile)
    target=only(t for t in JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String)).targets if String(t.curve_id)=="other07")
    reference_steps=[2e-4/2.0^j for j in (0,4,8,10)]
    rows,checks=NamedTuple[],NamedTuple[]
    variants=("loop256","weight_fine","guarded")
    for old in baseline
        u=Float64(old.delta_xi)
        any(abs(abs(u)-h*f)<1e-13*h for h in reference_steps,f in (.75,.875,1.,1.125,1.25)) || continue
        seed=SVector{5}(Float64(old.phi_u),Float64(old.phi_d),Float64(old.phi_s),Float64(old.Phi),Float64(old.Phibar))
        bg=asym_background(target,center+u,profile,seed)
        state=bg.state
        problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,state.qp,state.tp,state.kc;runtime.tau_kwargs...)
        prefactor=problem.plan.prefactor
        k=asym_smooth_kernel(state;n=128,degree=16)
        kfine=asym_smooth_kernel(state;n=256,degree=20)
        g=sigma_weight_table(state,prefactor,0.01;order=8)
        gf=sigma_weight_table(state,prefactor,0.01;fine=true,order=12)
        gt=sigma_threshold_g_oracle(state,prefactor;p_nodes=128,c_nodes=32)
        push!(checks,(delta_xi=u,check="continuous_weight_at_threshold",relative=gf(0.)/gt-1))
        push!(checks,(delta_xi=u,check="kernel_degree_and_loop",relative=(kfine.root-k.root)/max(abs(k.root),1e-30)))
        for tag in variants
            ker=tag=="loop256" ? kfine : k
            weight=tag=="weight_fine" ? gf : g
            root=(x=ker.root,width=ker.width)
            f=x->[x*weight(x)*asym_smooth_sigma(state,x,ker;
                regularized=tag=="guarded",n=tag=="weight_fine" ? 48 : 24),0.0]
            result=if tag=="guarded"
                # The safeguard is discontinuous at |det|=eps: split each crossing explicitly.
                edges=Float64[0.,0.01]
                if root.x>0 && ker.root_det_abs<SIG_PROP.PROPAGATOR_DENOM_EPS
                    fn=x->abs(complex(asym_polynomial(ker.shifted,(x-ker.root)/ker.xmax),
                        x^1.5*asym_polynomial(ker.coefficients[:,2],x/ker.xmax)))-SIG_PROP.PROPAGATOR_DENOM_EPS
                    width_floor=SIG_PROP.PROPAGATOR_DENOM_EPS/abs(ker.shifted[2]/ker.xmax)
                    for side in (-1.,1.)
                        a,b=sort([max(0.,root.x+side*8width_floor),root.x])
                        if fn(a)*fn(b)<0
                            fa=fn(a)
                            for _ in 1:60
                                m=(a+b)/2
                                if signbit(fn(m))==signbit(fa);a=m;fa=fn(m);else;b=m;end
                                b-a<1e-21 && break
                            end
                            push!(edges,(a+b)/2)
                        end
                    end
                end
                sort!(unique!(edges))
                val=zeros(2);err=zeros(2);leaves=0
                for (a,b) in zip(edges[1:end-1],edges[2:end])
                    # A separate map on each interval keeps all regulator shoulders.
                    localroot=(x=root.x-a,width=root.width)
                    part=asym_integral(y->f(a+y),localroot,b-a;atol=1e-21,rtol=2e-9)
                    val .+= part.value;err .+= part.error;leaves+=part.leaves
                end
                (value=val,error=err,leaves=leaves)
            else
                asym_integral(f,root,0.01;atol=1e-22,rtol=1e-12)
            end
            push!(rows,(point_id=String(old.point_id),variant=tag,xi=center+u,delta_xi=u,
                mixed_square_rate_fm2=result.value[1],quadrature_error_fm2=result.error[1],
                baseline_rate_fm2=Float64(old.mixed_square_rate_fm2),
                rate_difference_fm2=result.value[1]-Float64(old.mixed_square_rate_fm2),
                root_x_inv_fm2=ker.root,half_width_inv_fm2=ker.width,
                guard_active_at_peak=ker.root>0 && ker.root_det_abs<SIG_PROP.PROPAGATOR_DENOM_EPS,
                background_residual=bg.row.background_residual,leaves=result.leaves))
        end
        CSV.write(joinpath(OUTPUT_DIR,"reference_integrals.csv"),rows)
        CSV.write(joinpath(OUTPUT_DIR,"checks.csv"),checks)
        println("ASYM_CHECK u=",u," variants=",length(variants));flush(stdout)
        SIG_PC.reset_cache!();GC.gc()
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",xi_mott=center,
        source_summary_sha256=digest(joinpath(source,"summary.json")),steps=reference_steps,variants=variants,
        scope=String(meta.scope),production_kernel_changed=false,
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    check_asymptotics()
end
