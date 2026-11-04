#!/usr/bin/env julia
# Final independent background-grid and smaller derivative-step diagnostics.
include(joinpath(@__DIR__, isfile(joinpath(@__DIR__,"run_sigma_mott_asymptotics_v3.jl")) ?
    "run_sigma_mott_asymptotics_v3.jl" : "run_sigma_mott_asymptotics.jl"))

function background_asymptotic_checks()
    source=joinpath(dirname(OUTPUT_DIR),"asym_full_v2")
    rows=collect(CSV.File(joinpath(source,"integrals.csv")))
    meta=JSON3.read(read(joinpath(source,"summary.json"),String))
    center=Float64(meta.xi_mott)
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"s6.toml"))
    runtime=SP_OTHER.build_runtime(profile)
    target=only(t for t in JSON3.read(read(joinpath(INPUT_DIR,"targets.json"),String)).targets if String(t.curve_id)=="other07")
    near=rows[argmin(abs.(Float64.(getproperty.(rows,:delta_xi))))]
    seed=SVector{5}(Float64(near.phi_u),Float64(near.phi_d),Float64(near.phi_s),Float64(near.Phi),Float64(near.Phibar))
    fine0=asym_background(target,center,profile,seed;level="fine")
    k0=asym_smooth_kernel(fine0.state;n=256,degree=20)
    left=asym_background(target,center-1e-6,profile,seed;level="fine")
    right=asym_background(target,center+1e-6,profile,seed;level="fine")
    dr=(asym_smooth_kernel(right.state).root-asym_smooth_kernel(left.state).root)/2e-6
    fine_center=center-k0.root/dr
    fine_at=asym_background(target,fine_center,profile,seed;level="fine")
    fine_residual=asym_smooth_kernel(fine_at.state).root
    step_h=[2e-4/2.0^j for j in (0,4,8,10)]
    output=NamedTuple[]
    for variant in ("same_background_small_step","fine_background")
        c=variant=="fine_background" ? fine_center : center
        level=variant=="fine_background" ? "fine" : "base"
        offsets=sort!(unique([sg*h*f for sg in (-1.,1.),h in step_h,f in (.875,.9375,1.,1.0625,1.125)]))
        for u in offsets
            bg=asym_background(target,c+u,profile,seed;level=level)
            st=bg.state
            kernel=asym_smooth_kernel(st;n=256,degree=20)
            problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,st.qp,st.tp,st.kc;runtime.tau_kwargs...)
            weight=sigma_weight_table(st,problem.plan.prefactor,0.01;fine=true,order=12)
            root=(x=kernel.root,width=kernel.width)
            f=x->[x*weight(x)*asym_smooth_sigma(st,x,kernel;n=48),0.0]
            result=asym_integral(f,root,.01;atol=1e-22,rtol=1e-12)
            push!(output,(point_id=variant*"_"*string(length(output)+1),variant=variant,xi=c+u,delta_xi=u,
                mixed_square_rate_fm2=result.value[1],quadrature_error_fm2=result.error[1],
                root_x_inv_fm2=kernel.root,background_residual=bg.row.background_residual))
            CSV.write(joinpath(OUTPUT_DIR,"background_reference.csv"),output)
        end
        println("BACKGROUND_ASYM_VARIANT ",variant," completed");flush(stdout)
    end
    # Independent standard-form coefficients from the center kernel and thermal weight.
    b=asym_background(target,center,profile,seed)
    st=b.state;k=asym_smooth_kernel(st;n=256,degree=20)
    problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,st.qp,st.tp,st.kc;runtime.tau_kwargs...)
    g0=sigma_threshold_g_oracle(st,problem.plan.prefactor;p_nodes=128,c_nodes=32)
    m,tp=st.qp.m.u,st.tp;mu=st.qp.μ.u
    block=SIG_CS.combined_final_state_factor(:u,:ubar,m,m,m,m,mu,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ,0.)
    F=g0*abs2(k.coefficients[1,3])/(16pi*4m*m)*block
    br=k.coefficients[2,1]/k.xmax;ci=abs(k.coefficients[1,2])
    kl=asym_smooth_kernel(asym_background(target,center-1e-6,profile,seed).state)
    kr=asym_smooth_kernel(asym_background(target,center+1e-6,profile,seed).state)
    drdu=(kr.root-kl.root)/2e-6
    left_coefficient=.75*pi*F/(abs(br)*ci)*abs(drdu)^1.5
    right_log_coefficient=6F/br^2*drdu^2
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",steps=step_h,
        original_center=center,fine_center=fine_center,fine_center_shift=fine_center-center,
        fine_center_root_residual_inv_fm2=fine_residual,
        left_sqrt_h_coefficient=left_coefficient,right_log_coefficient=right_log_coefficient,
        root_slope=drdu,coefficient_origin="Independent threshold kernel/weight values; not fitted to curvature or rates.",
        derivative_steps=[.125,.0625],backgrounds=(base=(p=96,c=24,A_p=32,A_c=16,thermal_cutoff=25.),
            fine=(p=128,c=32,A_p=64,A_c=32,thermal_cutoff=30.)),
        production_kernel_changed=false,files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    background_asymptotic_checks()
end
