"""First-three-stage dependency audit: retain a finite-cutoff counting obstruction."""
module CausalGBUEndpointClosureAudit
include("causal_gbu_endpoint_closure.jl")
const E=CausalGBUEndpointClosure
const A=E.A
const R=E.R
using CSV,JSON3

function regular_part(rho,edges,S,c,width,z;nodes=12,atol=1e-8)
    left=S-width
    panels=sort!(unique!(vcat(edges,[left])))
    return A.integrate_panels(x->(rho(x)-(x>left ? c : 0.))/(x-z)/pi,
        panels;nodes=nodes,atol=atol)
end

function run_case(bg,ch,thermal)
    i,j=R.charged_rpa_spec(ch).pair
    a,u,b,v=bg.m[i],bg.mu[i],bg.m[j],bg.mu[j]
    ep=E.thermal_endpoint(a,u,b,v,bg.T,bg.Phi,bg.PhiBar,bg.vacuum,thermal)
    S,c,width=ep.lambda_endpoint_inv_fm,ep.rho_inside_inv_fm2,1.
    bound=E.projected_exterior(0.,a,b,bg.vacuum,thermal,bg.coupling[ch])
    grids=[R.build_spectral_bubble(0.,a,u,b,v,bg.T;Phi=bg.Phi,PhiBar=bg.PhiBar,
        vacuum_cutoff_inv_fm=bg.vacuum,thermal_cutoff_inv_fm=thermal,
        momentum_nodes=n,angle_nodes=16) for n in (128,256)]
    edges=A.cut_panels(last(grids))
    rhos=[x->A.direct_cut(bg,ch,0.,x,thermal;nodes=n).imaginary for n in (32,64)]
    regs=[regular_part(rhos[k],edges,S,c,width,S;nodes=k==1 ? 8 : 12,
        atol=k==1 ? 1e-7 : 1e-8) for k in 1:2]
    h=last(regs).value
    witness=E.endpoint_obstruction(c,S,width,h,bg.coupling[ch],bg.T,u-v)
    # Verify the endpoint subtraction away from the singular point independently.
    delta=0.01
    rz=regular_part(last(rhos),edges,S,c,width,S+delta)
    direct=A.continuous_cauchy(last(rhos),edges,complex(S+delta);atol=1e-8)
    reconstructed=rz.value+c/pi*log(delta/(width+delta))
    error=abs(reconstructed-direct.value)
    far=A.continuous_cauchy(last(rhos),edges,complex(bound.exterior_radius_inv_fm);atol=1e-8)
    farf=real(1-4bg.coupling[ch]*far.value)
    inside=[last(rhos)(S-d) for d in (1e-3,1e-5,1e-7)]
    # Only the ordinary, conservative low-energy analytic gap is used here.
    left=max(0.,abs(a-b)-u+v)
    right=a+b-u+v
    results=map(grids) do g
        inv(w)=1-4bg.coupling[ch]*R.spectral_bubble(g,w;eta_inv_fm=0.).value
        R.certify_gap_roots((w,_)->inv(w),0.,[(left,right)];physical_sheet=true,
            real_axis=true,omega_nodes=128)
    end
    roots=last(results).roots
    nodal=first(results).count==last(results).count ?
        maximum((abs(first(results).roots[k].omega_inv_fm-r.omega_inv_fm) for (k,r) in enumerate(roots));init=0.) : Inf
    p=R.build_bubble_dispersion(last(grids);segment_nodes=256,energy_nodes=64)
    f_interpolant=real(1-4bg.coupling[ch]*R.cauchy_transform(p,S))
    normalrows=NamedTuple[]
    for r in roots
        slope=real(-4bg.coupling[ch]*R.spectral_bubble(last(grids),r.omega_inv_fm;eta_inv_fm=0.).derivative)
        weight=-2bg.coupling[ch]/slope
        push!(normalrows,(channel=String(ch),thermal_cutoff_inv_fm=thermal,
            root_k0_inv_fm=r.omega_inv_fm,slope_inv_fm_negative=slope<0,
            pole_spectral_weight_fm=weight,positive_spectral_weight=weight>0,
            root_internal_node_change_inv_fm=nodal,production_authorized=false))
    end
    passed=all(r.converged for r in regs) && rz.converged && direct.converged && far.converged &&
        abs(regs[1].value-regs[2].value)<1e-6 && error<1e-6 && farf>0 &&
        all(r.passed for r in results) && nodal<1e-7 && c<0
    row=merge((channel=String(ch),q_inv_fm=0.,thermal_cutoff_inv_fm=thermal,
        lambda_endpoint_inv_fm=S,k0_endpoint_inv_fm=ep.k0_endpoint_inv_fm,
        rho_inside_limit_inv_fm2=c,rho_offset_1e3=inside[1],rho_offset_1e5=inside[2],
        rho_offset_1e7=inside[3],regular_part_at_endpoint_inv_fm2=h,
        regular_part_node_change_inv_fm2=abs(regs[1].value-regs[2].value),
        subtraction_identity_error_inv_fm2=error,far_inverse_real=farf,
        interpolant_endpoint_inverse=f_interpolant,normal_gap_positive_roots=length(roots),
        endpoint_witness_numerics_passed=passed,first_three_steps_completed=false,
        full_spectrum_counting_accepted=false,meson_density_computed=false),bound,witness)
    return (row=row,normalrows=normalrows)
end

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_ENDPOINT_OUTPUT",joinpath(base,"direction_b_endpoint_closure_20260907"))
    R.method_contract()
    hashes=R.start_output(output)
    rows,normal,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    for ch in R.CHANNELS,thermal in (10.,20.,24.)
        try
            r=run_case(bg,ch,thermal)
            push!(rows,r.row);append!(normal,r.normalrows)
            println("[endpoint] $(ch) Lth=$(thermal) log10_distance=$(r.row.asymptotic_log10_distance_inv_fm) normal=$(r.row.normal_gap_positive_roots) passed=$(r.row.endpoint_witness_numerics_passed)")
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(channel=String(ch),thermal_cutoff_inv_fm=thermal,reason=sprint(showerror,err)))
            println(stderr,"[endpoint-failed] $(last(failures))")
        end
        for (name,rs) in (("endpoint_witnesses",rows),("normal_gap_roots",normal),("failures",failures))
            isempty(rs) || CSV.write(joinpath(output,name*".csv"),rs)
        end
        flush(stdout)
    end
    all(R.hashfile(joinpath(bg.input_directory,f))==h for (f,h) in bg.input_hashes) || error("input drift")
    complete=isempty(failures) && length(rows)==12
    R.finish_output(output,hashes,Dict("status"=>complete ? "finite_thermal_endpoint_counting_obstruction" : "evaluation_failed",
        "background"=>bg,"thermal_cutoffs_inv_fm"=>[10.,20.,24.],
        "witness_rows"=>length(rows),"normal_root_rows"=>length(normal),
        "witness_numerics_passed"=>complete && all(r.endpoint_witness_numerics_passed for r in rows),
        "stage_status"=>E.stage_status(analytic_structure=false,physical_counting=false,integral_convergence=false),
        "stage1_failure_reason"=>"one_sided_negative_thermal_hard_edge_requires_exterior_zero_missing_from_smoothed_interpolant_count",
        "new_upper_half_plane_instability_found"=>false,"solver_called"=>false,"meson_density_computed"=>false,
        "cut_nodes"=>[32,64],"regular_part_atols_inv_fm2"=>[1e-7,1e-8],
        "normal_root_loop_nodes"=>[128,256],"normal_root_scan_nodes"=>128,
        "limitations"=>["Endpoint theorem is q=0; q=0 has zero external phase space",
            "Existence is exact under stated endpoint conditions; log distances are asymptotic estimates",
            "No uniqueness or total density bound is implied by a per-root Bose bound",
            "A real cutoff-induced zero is not a physical bound state or evidence of UHP instability",
            "Finite Lth approximates an intended infinite thermal integral, but root counting is nonuniform in this limit",
            "No cutoff change, spectral clipping, phase folding or production promotion"]))
    complete || error("endpoint audit incomplete; evidence retained")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
