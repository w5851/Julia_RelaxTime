"""Fixed-background PNJL, routing and RPA pole-health audit.

Numerical completion is separate from physical acceptance. Every failure is
retained. No new equilibrium/density scan, old default or baseline is changed.
"""
module CausalGBUPNJLStabilityAudit
include("causal_gbu_pnjl_stability.jl")
const H=CausalGBUPNJLStability
const R=H.R
const G=H.G
using CSV,JSON3

function settings()
    mesh=parse(Int,get(ENV,"GBU_HEALTH_MESH","256"))
    np=parse(Int,get(ENV,"GBU_HEALTH_NODES","256"))
    qs=sort!(unique!(parse.(Float64,split(get(ENV,"GBU_HEALTH_Q","0,1.4,3.2,6.2"),','))))
    channels=unique(Symbol.(split(get(ENV,"GBU_HEALTH_CHANNELS","pi_plus,pi_minus,K_plus,K_minus"),',')))
    variants=unique(split(get(ENV,"GBU_HEALTH_VARIANTS","all_hard,thermal24"),','))
    min(mesh,np)>=32 && iseven(mesh) && iseven(np) && !isempty(qs) &&
        first(qs)==0 && all(isfinite,qs) && !isempty(channels) && !isempty(variants) &&
        all(c->c in R.CHANNELS,channels) && all(v->v in ("all_hard","thermal24"),variants) ||
        throw(ArgumentError("invalid health settings"))
    return (mesh=mesh,np=np,qs=qs,channels=channels,variants=variants,etas=(0.03,0.01,0.003),ne=128)
end

function q0_contact(bg,ch,thermal,n)
    a,b=R.charged_rpa_spec(ch).pair
    ps,ws=R.gauleg(0.,thermal,n)
    independent=Main.RelaxTime.OneLoopIntegrals.A(bg.m[a],bg.mu[a],bg.T,bg.Phi,bg.PhiBar,ps,ws)+
        Main.RelaxTime.OneLoopIntegrals.A(bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,ps,ws)
    # A's vacuum term is always Lambda; match the full thermal domain explicitly.
    grid=R.build_spectral_bubble(0.,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=thermal,momentum_nodes=n,angle_nodes=32)
    error=abs(grid.contact_inv_fm2-independent)
    return (channel=String(ch),thermal_cutoff_inv_fm=thermal,independent_A_sum=independent,
        contact=grid.contact_inv_fm2,error=error,passed=error<1e-9,production_authorized=false)
end

function background_gap_rows(bg,n)
    c=Main.Constants_PNJL
    bare=(c.m_ud0_inv_fm,c.m_ud0_inv_fm,c.m_s0_inv_fm)
    flavors=(:u,:d,:s)
    rows=NamedTuple[]
    # The default upstream thermal ceiling is recorded separately from the
    # all-hard and extended-loop candidates. No coupling is written back.
    for thermal in unique([bg.vacuum,c.thermal_p_max_inv_fm,24.]),nodes in (div(n,2),n)
        ps,ws=R.gauleg(0.,thermal,nodes)
        A=ntuple(i->Main.RelaxTime.OneLoopIntegrals.A(bg.m[flavors[i]],bg.mu[flavors[i]],
            bg.T,bg.Phi,bg.PhiBar,ps,ws),3)
        state=H.tadpole_gap_state(Tuple(bg.m),A,bare,c.G_fm2,c.K_fm5;Nc=c.N_color)
        for (i,flavor) in enumerate(flavors)
            push!(rows,(flavor=String(flavor),thermal_cutoff_inv_fm=thermal,nodes=nodes,
                A_inv_fm2=A[i],condensate_inv_fm3=state.condensates_inv_fm3[i],
                frozen_mass_inv_fm=bg.m[flavor],predicted_mass_inv_fm=state.masses_inv_fm[i],
                gap_residual_inv_fm=state.masses_inv_fm[i]-bg.m[flavor],
                K12_difference_fm2=state.K12_fm2-bg.coupling[:pi_plus],
                K45_difference_fm2=state.K45_fm2-bg.coupling[:K_plus],
                upstream_default_thermal=thermal==c.thermal_p_max_inv_fm,
                full_stationarity_certified=false,solver_called=false,production_authorized=false))
        end
    end
    return rows
end

function health_case(bg,ch,q,thermal,s)
    a,b=R.charged_rpa_spec(ch).pair
    shift=bg.mu[a]-bg.mu[b]
    K=bg.coupling[ch]
    g=R.build_spectral_bubble(q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=thermal,momentum_nodes=s.np,angle_nodes=div(s.np,2))
    pole_rows,contour_rows,mesh_rows=NamedTuple[],NamedTuple[],NamedTuple[]
    roots_by_mesh=Vector{Vector{Float64}}()
    for mesh in (div(s.mesh,2),s.mesh)
        p=R.build_bubble_dispersion(g;segment_nodes=mesh,energy_nodes=s.ne)
        roots=G.count_all_real_gaps(p,K,q;root_nodes=128)
        lambdas=sort!([parse(Float64,l) for r in roots.rows for l in split(r.roots_lambda_inv_fm,';') if !isempty(l)])
        push!(roots_by_mesh,lambdas)
        f(z)=1-4K*R.cauchy_transform(p,z)
        residues_ok=true
        derivative_difference=0.
        for lambda in lambdas
            # Choose a finite-difference step INSIDE this certified open gap.
            row=only(filter(r->r.left_lambda_inv_fm<lambda<r.right_lambda_inv_fm,roots.rows))
            exact=-4K*real(H.spectral_derivative(p,lambda))
            difference=H.gap_derivative_check(f,lambda,row.left_lambda_inv_fm,row.right_lambda_inv_fm,exact)
            numerical,error=difference.value,difference.error
            derivative_difference=max(derivative_difference,error)
            weight=H.pole_weight(exact,K,lambda-shift)
            residues_ok &= weight.sign_passed && difference.passed
            push!(pole_rows,(mesh=mesh,lambda_inv_fm=lambda,k0_inv_fm=lambda-shift,
                inverse_residual=abs(f(lambda)),derivative=exact,finite_difference_derivative=numerical,
                derivative_relative_error=error,weight_fm=weight.weight_fm,
                derivative_initial_error=difference.initial_error,
                derivative_step_change=difference.step_change,derivative_halvings=difference.halvings,
                derivative_step_inv_fm=difference.step_inv_fm,derivative_check_passed=difference.passed,
                simple=weight.simple,sign_passed=weight.sign_passed,production_authorized=false))
        end
        passivity=H.passivity_audit(p,K,shift)
        anchors=collect(p.energy[1:4:end])
        for l in lambdas
            append!(anchors,[l+d for d in (-0.03,-0.01,-0.003,0.,0.003,0.01,0.03)])
        end
        uhp_clear=true
        for eta in s.etas
            c=H.upper_rectangle_count(f,roots.tail.radius_inv_fm,eta;anchors=anchors)
            uhp_clear &= c.resolved && c.count==0
            push!(contour_rows,merge((mesh=mesh,tail_radius_lambda_inv_fm=roots.tail.radius_inv_fm,
                exterior_inverse_bound=roots.tail.inverse_deviation_bound),c,(production_authorized=false,)))
        end
        static=R.spectral_bubble_cut(g,0.;energy_nodes=s.ne)
        push!(mesh_rows,(mesh=mesh,real_gap_checks_passed=roots.passed,real_root_count=roots.count,
            pole_weight_checks_passed=residues_ok,max_derivative_relative_error=derivative_difference,
            UHP_tested_bands_clear=uhp_clear,static_cut_imag=static.imaginary,
            static_cut_passed=abs(static.imaginary)<1e-10,
            minimum_k0_rho=passivity.minimum_k0_rho,negative_spectral_cells=passivity.negative_segments,
            numerical_passivity=passivity.numerical_passivity,
            sign_sufficient_condition_float64=passivity.sign_sufficient_condition_float64,
            static_inverse_real=passivity.static_inverse_real,
            static_inverse_imag=passivity.static_inverse_imag,
            full_continuum_UHP_certified=false,production_authorized=false))
    end
    stable=length(roots_by_mesh[1])==length(roots_by_mesh[2])
    drift=stable ? maximum(abs.(roots_by_mesh[1].-roots_by_mesh[2]);init=0.) : Inf
    fine=last(mesh_rows)
    diagnostic_ok=stable && all(r.real_gap_checks_passed && r.pole_weight_checks_passed &&
        r.UHP_tested_bands_clear && r.static_cut_passed && r.static_inverse_real>0 for r in mesh_rows)
    summary=(real_root_counts_stable=stable,max_root_mesh_drift_inv_fm=drift,
        diagnostic_checks_passed=diagnostic_ok,numerical_passivity=fine.numerical_passivity,
        full_continuum_stability_certified=false,
        production_blockers=fine.numerical_passivity ?
            "near_axis_strip;continuous_kernel;Mott_Levinson;density_convergence" :
            "signed_thermal_extension_spectrum;near_axis_strip;continuous_kernel;Mott_Levinson;density_convergence",
        density_computed=false,production_authorized=false)
    return (;summary,pole_rows,contour_rows,mesh_rows)
end

function routing_case(bg,ch,q,s)
    a,b=R.charged_rpa_spec(ch).pair
    zs=ComplexF64[-1.0+0.6im,0.4+0.6im,3.2+0.8im] # INTERNAL lambda
    args=(q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,bg.vacuum,24.,zs)
    angular=H.single_sphere_loop(args...;nodes=s.np,coordinates=:angular)
    radial=H.single_sphere_loop(args...;nodes=s.np,coordinates=:radial)
    coarse=H.single_sphere_loop(args...;nodes=div(s.np,2),coordinates=:angular)
    reverse=H.single_sphere_loop(q,bg.m[b],bg.mu[b],bg.m[a],bg.mu[a],bg.T,
        bg.Phi,bg.PhiBar,bg.vacuum,24.,-conj.(zs);nodes=s.np)
    ps,ws=R.gauleg(0.,24.,s.np)
    A1=Main.RelaxTime.OneLoopIntegrals.A(bg.m[a],bg.mu[a],bg.T,bg.Phi,bg.PhiBar,ps,ws)
    A2=Main.RelaxTime.OneLoopIntegrals.A(bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,ps,ws)
    numerical=max(maximum(abs.(angular.values.-radial.values)),abs(angular.contact-radial.contact))
    node=maximum(abs.(angular.values.-coarse.values))
    reflection=maximum(abs.(angular.values.-conj.(reverse.values)))
    zero=q==0 ? R.build_spectral_bubble(0.,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=24.,
        momentum_nodes=s.np,angle_nodes=div(s.np,2)) : nothing
    q0error=zero===nothing ? NaN : maximum(abs(angular.values[i]-
        R.spectral_bubble(zero,real(z)-bg.mu[a]+bg.mu[b];eta_inv_fm=imag(z)).value) for (i,z) in enumerate(zs))
    return (coordinate_error=numerical,node_change=node,flavor_reflection_error=reflection,
        A1_minus_equilibrium=angular.A1-A1,A2_minus_equilibrium=angular.A2-A2,
        q0_two_line_error=q0error,numerics_passed=numerical<1e-6 && node<1e-6,
        flavor_reflection_passed=reflection<1e-9,
        status=reflection<1e-9 ? "probe_passed_not_density_certified" : "routing_symmetry_failed",
        is_legacy_implementation=false,density_computed=false,production_authorized=false)
end

function main()
    s=settings()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_HEALTH_OUTPUT",joinpath(base,"pnjl_routing_stability_20260906"))
    R.method_contract()
    hashes=R.start_output(output)
    summaries,poles,contours,meshes,qzeros,routing,failures=(NamedTuple[] for _ in 1:7)
    save(name,rows)=isempty(rows) ? nothing : CSV.write(joinpath(output,name*".csv"),rows)
    save("background_gap",background_gap_rows(bg,s.np))
    for ch in s.channels
        for thermal in (bg.vacuum,24.)
            push!(qzeros,q0_contact(bg,ch,thermal,s.np))
            save("q0_contacts",qzeros)
        end
        for q in s.qs
            meta=(channel=String(ch),q_inv_fm=q)
            try
                push!(routing,merge(meta,routing_case(bg,ch,q,s)))
                println("[health-routing] $(ch) q=$(q) status=$(last(routing).status)")
            catch err
                err isa InterruptException && rethrow()
                push!(failures,merge(meta,(variant="single_sphere",reason=sprint(showerror,err))))
            end
            save("routing",routing)
            for (variant,thermal) in (("all_hard",bg.vacuum),("thermal24",24.))
                variant in s.variants || continue
                try
                    result=health_case(bg,ch,q,thermal,s)
                    header=merge(meta,(variant=variant,thermal_cutoff_inv_fm=thermal))
                    push!(summaries,merge(header,result.summary))
                    append!(poles,[merge(header,r) for r in result.pole_rows])
                    append!(contours,[merge(header,r) for r in result.contour_rows])
                    append!(meshes,[merge(header,r) for r in result.mesh_rows])
                    println("[health] $(ch) q=$(q) $(variant) checks=$(result.summary.diagnostic_checks_passed) passivity=$(result.summary.numerical_passivity)")
                catch err
                    err isa InterruptException && rethrow()
                    push!(failures,merge(meta,(variant=variant,reason=sprint(showerror,err))))
                    println(stderr,"[health-failed] $(ch) q=$(q) $(variant): $(sprint(showerror,err))")
                end
                for (name,rows) in (("summary",summaries),("poles",poles),("contours",contours),
                                    ("mesh_checks",meshes),("failures",failures))
                    save(name,rows)
                end
                flush(stdout)
            end
        end
    end
    complete=isempty(failures) && length(summaries)==length(s.channels)*length(s.variants)*length(s.qs) &&
        length(routing)==length(s.channels)*length(s.qs)
    R.finish_output(output,hashes,Dict("status"=>complete ? "diagnostic_completed_not_production_accepted" : "evaluation_failed",
        "all_evaluations_completed"=>complete,"settings"=>s,"background"=>bg,
        "q0_passed"=>all(r.passed for r in qzeros),
        "routing_numerics_passed"=>all(r.numerics_passed for r in routing),
        "routing_flavor_reflection_passed"=>all(r.flavor_reflection_passed for r in routing),
        "main_diagnostic_checks_passed"=>all(r.diagnostic_checks_passed for r in summaries),
        "all_numerical_passivity_passed"=>all(r.numerical_passivity for r in summaries),
        "full_continuum_stability_certified"=>false,"solver_called"=>false,
        "density_computed"=>false,"production_authorized"=>false,
        "limitations"=>["All-hard is not the thermal-extended upstream equilibrium",
            "UHP contours leave 0<Im(lambda)<0.003 unresolved and test spectral interpolants",
            "Signed spectral contribution is not the same as signed GBU continuum density",
            "First-line sphere is an explicit routing probe, not the historical legacy B0",
            "No clipping, flavor averaging, phase folding or production promotion"] ))
    println("[health] output=$(output) completed=$(complete); production=false")
    complete || error("health audit evaluations failed; evidence retained")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
