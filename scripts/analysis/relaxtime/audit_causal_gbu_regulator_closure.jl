"""Same-background regulator closure audit, BEFORE any alternative density claim.

The centered regulator is a sensitivity model, not a promoted provider. A weak
Poisson-kernel test checks a distributional eta limit, not pointwise PV parity.
"""
module CausalGBURegulatorClosure
include("causal_gbu_regulator_checks.jl")
const C=CausalGBURegulatorChecks
const R=C.R
using CSV, JSON3

function main(;configurations=vcat([(q,mesh) for q in (0.0,1.0) for mesh in (128,256,512)],[(3.0,256),(6.0,256)]),
              np=256,nx=128)
    np>=8 && nx>=8 && iseven(np) && iseven(nx) || throw(ArgumentError("even loop nodes >=8 required"))
    !isempty(configurations) && all(isfinite(q) && q>=0 && mesh>=8 for (q,mesh) in configurations) ||
        throw(ArgumentError("invalid closure configurations"))
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=joinpath(base,"negative_density_phase_fig2_like")
    output=get(ENV,"GBU_CLOSURE_OUTPUT",joinpath(base,"regulator_routing_closure_20260905"))
    bg=R.frozen_background(input)
    hashes=R.start_output(output)
    comparisons,values,gaps,weak,failures=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    save(name,rows)=isempty(rows) ? nothing : CSV.write(joinpath(output,name*".csv"),rows)
    for channel in R.CHANNELS
        a,b=R.charged_rpa_spec(channel).pair
        m1,m2,mu1,mu2=bg.m[a],bg.m[b],bg.mu[a],bg.mu[b]
        shift=mu1-mu2
        zero=C.centered_atoms(0.0,m1,mu1,m2,mu2,bg.T,bg.Phi,bg.PhiBar,bg.vacuum,24.;np=np,nx=nx)
        for (q,mesh) in configurations
            s=R.Settings(mesh=mesh,ne=128,np=np,nx=nx,thermal=24.,nw=2400)
            for regulator in (:two_line,:centered)
                try
                    bb=regulator===:two_line ? R.bubble_at(bg,channel,q,s) : C.centered_profile(bg,channel,q;mesh=mesh,ne=128)
                    cuts=C.support_intervals(q,m1,m2,bg.vacuum,24.,regulator)
                    atom=regulator===:two_line ? bb.grid : C.centered_atoms(q,m1,mu1,m2,mu2,bg.T,bg.Phi,bg.PhiBar,bg.vacuum,24.;np=np,nx=nx)
                    atom_coarse=regulator===:two_line ? R.build_spectral_bubble(q,m1,mu1,m2,mu2,bg.T;Phi=bg.Phi,PhiBar=bg.PhiBar,
                        vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=24.,momentum_nodes=div(np,2),angle_nodes=div(nx,2)) :
                        C.centered_atoms(q,m1,mu1,m2,mu2,bg.T,bg.Phi,bg.PhiBar,bg.vacuum,24.;np=div(np,2),nx=div(nx,2))
                    function direct(g,w,eta)
                        if regulator===:two_line
                            v=R.spectral_bubble(g,w;eta_inv_fm=eta)
                            return (value=v.value,residual=v.contact_identity_residual,contact=g.contact_inv_fm2)
                        end
                        v=C.atom_value(g,complex(w+shift,eta))
                        return (value=v.value,residual=v.contact_residual,contact=g.contact)
                    end
                    contact=direct(atom,0.0,0.4).contact
                    max_causal,max_node,max_contact=0.0,0.0,0.0
                    probes=(0.0,0.5*(bb.landau+bb.threshold),bb.threshold+0.5,bb.threshold+2.0,12.0)
                    for w in probes
                        v=direct(atom,w,0.4)
                        vc=direct(atom_coarse,w,0.4)
                        sp=R.cauchy_transform(bb.profile,w+shift+0.4im)
                        max_causal=max(max_causal,abs(v.value-sp))
                        max_node=max(max_node,abs(v.value-vc.value))
                        max_contact=max(max_contact,v.residual)
                        inv=bb.inverse(w)
                        push!(values,(channel=String(channel),regulator=String(regulator),q_inv_fm=q,mesh=mesh,k0_inv_fm=w,
                            lambda_inv_fm=w+shift,pv_inverse_real=real(inv),pv_inverse_imag=imag(inv),phase_over_pi=-angle(inv)/pi,
                            eta_inv_fm=0.4,direct_pi_real=real(v.value),direct_pi_imag=imag(v.value),
                            spectral_pi_real=real(sp),spectral_pi_imag=imag(sp),complex_difference=abs(v.value-sp),
                            atom_node_change=abs(v.value-vc.value),contact_identity_residual=v.residual,production_authorized=false))
                    end
                    full=C.count_support_gaps(bb,s,cuts)
                    append!(gaps,[merge((channel=String(channel),regulator=String(regulator),q_inv_fm=q,mesh=mesh),r) for r in full])
                    # The independently counted geometric gaps already cover this
                    # normal gap; do not repeat an identical contour evaluation.
                    normal=R.gap_audit(bb,s;count_contour=false)
                    phase=R.threshold_phase_limit(bb.inverse,bb.threshold)
                    all_count=sum(r.root_count for r in full)
                    count_ok=all(r.passed for r in full) && normal.passed && all_count==normal.count
                    tail=bb.inverse(64.0)
                    phase_ok=phase.passed && abs(phase.phase/pi-normal.count)<s.phase_tol && abs(angle(tail))<s.phase_tol
                    root=isempty(normal.roots) ? NaN : first(normal.roots).omega_inv_fm
                    # Small fixed representation target, not a relaxed physical phase gate.
                    causal_ok=max_causal<1e-4 && max_node<1e-6 && max_contact<1e-10
                    push!(comparisons,(channel=String(channel),regulator=String(regulator),q_inv_fm=q,mesh=mesh,
                        normal_gap_root_inv_fm=root,normal_count=normal.count,all_support_gap_count=all_count,
                        all_support_gaps_passed=all(r.passed for r in full),normal_gap_passed=normal.passed,
                        threshold_phase_over_pi=phase.phase/pi,phase_passed=phase_ok,
                        static_inverse_real=real(bb.inverse(0.)),static_inverse_imag=imag(bb.inverse(0.)),tail_phase=-angle(tail),
                        contact_inv_fm2=contact,contact_q0_inv_fm2=zero.contact,
                        wrong_q0_contact_inverse_shift=-4bg.coupling[channel]*C.NC/(8pi^2)*(contact-zero.contact),
                        max_causal_difference=max_causal,max_atom_node_change=max_node,max_contact_residual=max_contact,
                        causal_passed=causal_ok,conditional_window_passed=causal_ok && count_ok && phase_ok,
                        density_computed=false,production_authorized=false))
                    if mesh==512 && q==1.0
                        # Poisson semigroup: int h_gamma(l-c) Im Pi(l+i eta) dl
                        # = Im Pi(c+i(gamma+eta)); target eta=0 has gamma>0.
                        center=bb.threshold+shift+1.0
                        gamma=0.7
                        target=imag(R.cauchy_transform(bb.profile,center+gamma*im))
                        for eta in (0.2,0.1,0.05,0.025,0.0125,0.00625)
                            v=imag(R.cauchy_transform(bb.profile,center+(gamma+eta)*im))
                            oracle=imag(direct(atom,center-shift,gamma+eta).value)
                            push!(weak,(channel=String(channel),regulator=String(regulator),q_inv_fm=q,eta_inv_fm=eta,
                                test_center_lambda_inv_fm=center,test_width_inv_fm=gamma,pv_weak_integral=target,
                                eta_weak_integral=v,absolute_eta_error=abs(v-target),direct_loop_weak_value=oracle,
                                representation_difference=abs(v-oracle),test="whole_axis_Poisson_weight_not_GBU_density",production_authorized=false))
                        end
                    end
                    println("[gbu-closure] $(channel) $(regulator) q=$(q) mesh=$(mesh) root=$(root) causal=$(max_causal) gaps=$(count_ok)")
                catch err
                    err isa InterruptException && rethrow()
                    push!(failures,(channel=String(channel),regulator=String(regulator),q_inv_fm=q,mesh=mesh,reason=sprint(showerror,err)))
                    println(stderr,"[gbu-closure-failed] $(channel) $(regulator) q=$(q): $(sprint(showerror,err))")
                end
                save("closure_summary",comparisons);save("complex_profiles",values);save("support_gap_counts",gaps)
                save("weak_eta_convergence",weak);save("failures",failures);flush(stdout)
            end
        end
    end
    R.finish_output(output,hashes,Dict("status"=>"diagnostic_closure_checks_completed","background"=>bg,
        "configurations"=>configurations,"solver_called"=>false,"density_computed"=>false,
        "momentum_nodes"=>np,"angle_nodes"=>nx,"coarse_momentum_nodes"=>div(np,2),"coarse_angle_nodes"=>div(nx,2),
        "regulators"=>["two_line","centered"],"representation_target"=>1e-4,"direct_node_target"=>1e-6,
        "failed_evaluations"=>length(failures),"gap_count_window_inv_fm"=>[0.0,64.0],"gap_endpoint_margin_inv_fm"=>1e-6,
        "full_spectrum_certified"=>false,"limitations"=>[
            "Root counts cover geometric support gaps with endpoint margins, not cut-embedded or UHP zeros",
            "Weak eta test is for the bubble spectrum, not an alternative GBU density",
            "No centered density is authorized before all relevant q shells and endpoint/mesh gates close"],
        "production_authorized"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
