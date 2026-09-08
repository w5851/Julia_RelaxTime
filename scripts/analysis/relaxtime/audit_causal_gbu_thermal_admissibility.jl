"""Diagnostic-only KMS sign witnesses and direct-continuum PV/near-axis probes.

No equilibrium solves, new density integration or automatic regulator choice.
Finite sampled kernel agreement cannot transfer a UHP zero count by Rouche.
"""
module CausalGBUThermalAdmissibilityAudit
include("causal_gbu_thermal_admissibility.jl")
const A=CausalGBUThermalAdmissibility
const R=A.R
using CSV,JSON3

function settings()
    return (channels=R.CHANNELS,occupations=("fermi","pnjl"),
        witness_q=(0.,1.4,3.2,6.2),thermal_inv_fm=24.,
        probe_channels=(:pi_plus,:K_plus),probe_q=(0.,3.2),
        meshes=(128,256,512),etas=(0.,0.003,0.001,0.0003),
        cut_nodes=(32,64),quad_nodes=(8,12),atols=(1e-7,1e-8),
        direct_node_tolerance_inv_fm2=1e-6,cut_tolerance_inv_fm2=1e-9)
end

function grid(bg,ch,q,thermal;Phi=bg.Phi,PhiBar=bg.PhiBar)
    a,b=R.charged_rpa_spec(ch).pair
    return R.build_spectral_bubble(q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
        Phi=Phi,PhiBar=PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=thermal,momentum_nodes=128,angle_nodes=64)
end

function witnesses(bg,s)
    rows=NamedTuple[]
    for ch in s.channels,occupation in s.occupations,q in s.witness_q
        a,b=R.charged_rpa_spec(ch).pair
        phi,bar=occupation=="fermi" ? (1.,1.) : (bg.Phi,bg.PhiBar)
        g=grid(bg,ch,q,s.thermal_inv_fm;Phi=phi,PhiBar=bar)
        h=grid(bg,ch,q,bg.vacuum;Phi=phi,PhiBar=bar)
        # q=0 has an exact radial on-shell formula. At q>0, pick an open
        # interval above the entire vacuum pair support and unitary onset.
        floor=max(hypot(bg.vacuum,bg.m[a])+hypot(bg.vacuum,bg.m[b]),hypot(q,bg.m[a]+bg.m[b]))
        for offset in (0.05,0.25,1.)
            analytic=q==0 ? A.q0_pair(bg.vacuum+offset,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],
                bg.T,phi,bar,bg.vacuum,s.thermal_inv_fm) : nothing
            lambda=analytic===nothing ? floor+offset : analytic.lambda_inv_fm
            k0=lambda-bg.mu[a]+bg.mu[b]
            direct=A.direct_cut(bg,ch,q,lambda,s.thermal_inv_fm;Phi=phi,PhiBar=bar,nodes=64)
            coarse=A.direct_cut(bg,ch,q,lambda,s.thermal_inv_fm;Phi=phi,PhiBar=bar,nodes=32)
            current=R.spectral_bubble_cut(g,k0;energy_nodes=64)
            hard=R.spectral_bubble_cut(h,k0;energy_nodes=64)
            discrepancy=abs(direct.imaginary-current.imaginary)
            analytic_error=analytic===nothing ? NaN : abs(analytic.imaginary_inv_fm2-direct.imaginary)
            tested_negative=k0>0 && direct.vacuum_pair==0 && direct.landau==0 &&
                direct.thermal_pair<0 && current.imaginary<0 && hard.imaginary==0
            push!(rows,(channel=String(ch),occupation=occupation,q_inv_fm=q,offset_inv_fm=offset,
                lambda_inv_fm=lambda,k0_inv_fm=k0,vacuum_pair_inv_fm2=direct.vacuum_pair,
                thermal_pair_inv_fm2=direct.thermal_pair,landau_inv_fm2=direct.landau,
                original_cut_inv_fm2=current.imaginary,hard_cut_inv_fm2=hard.imaginary,
                cut_error_inv_fm2=discrepancy,node_change_inv_fm2=abs(direct.imaginary-coarse.imaginary),
                analytic_q0_error_inv_fm2=analytic_error,negative_tail_witness=tested_negative,
                numerics_passed=discrepancy<s.cut_tolerance_inv_fm2 &&
                    abs(direct.imaginary-coarse.imaginary)<s.cut_tolerance_inv_fm2 &&
                    (analytic===nothing || analytic_error<s.cut_tolerance_inv_fm2),
                physical_full_correlator_sign_passed=false,production_authorized=false))
        end
    end
    return rows
end

function probes(bg,ch,q,thermal,s)
    a,b=R.charged_rpa_spec(ch).pair
    g=grid(bg,ch,q,thermal)
    edges=A.cut_panels(g)
    profiles=[R.build_bubble_dispersion(g;segment_nodes=mesh,energy_nodes=128) for mesh in s.meshes]
    low,threshold=hypot(q,bg.m[a]-bg.m[b]),hypot(q,bg.m[a]+bg.m[b])
    vacuum_top=hypot(bg.vacuum,bg.m[a])+hypot(bg.vacuum,bg.m[b])
    centers=(gap=(low+threshold)/2,unitary=(threshold+vacuum_top)/2,tail=vacuum_top+0.25)
    rows=NamedTuple[]
    for (region,center) in pairs(centers),eta in s.etas
        z=complex(center,eta)
        results=map(1:2) do i
            rho=x->A.direct_cut(bg,ch,q,x,thermal;nodes=s.cut_nodes[i]).imaginary
            A.continuous_cauchy(rho,edges,z;nodes=s.quad_nodes[i],atol=s.atols[i])
        end
        coarse,fine=results
        interps=[R.cauchy_transform(p,z) for p in profiles]
        errors=abs.(interps.-fine.value)
        finite=R.spectral_bubble(g,center-bg.mu[a]+bg.mu[b];eta_inv_fm=max(eta,0.03))
        spectral=A.rpa_spectral(fine.value,bg.coupling[ch])
        # These are sampled discrepancies, not a uniform contour bound.
        margin=A.count_transfer_margin(spectral.inverse_abs,4bg.coupling[ch]*last(errors))
        push!(rows,(channel=String(ch),q_inv_fm=q,variant=thermal==bg.vacuum ? "all_hard" : "thermal24",
            region=String(region),lambda_inv_fm=center,k0_inv_fm=center-bg.mu[a]+bg.mu[b],eta_inv_fm=eta,
            direct_real_inv_fm2=real(fine.value),direct_imag_inv_fm2=imag(fine.value),
            error_estimate_inv_fm2=fine.error_estimate,node_change_inv_fm2=abs(fine.value-coarse.value),
            quadrature_converged=coarse.converged && fine.converged,
            direct_numerics_passed=coarse.converged && fine.converged &&
                abs(fine.value-coarse.value)<s.direct_node_tolerance_inv_fm2,
            mesh128_error_inv_fm2=errors[1],mesh256_error_inv_fm2=errors[2],mesh512_error_inv_fm2=errors[3],
            mesh128_to256_inv_fm2=abs(interps[1]-interps[2]),mesh256_to512_inv_fm2=abs(interps[2]-interps[3]),
            imag_D_fm2=spectral.imaginary,rpa_identity_error_fm2=abs(spectral.imaginary-spectral.identity_value),
            contact_algebra_error_inv_fm2=finite.contact_identity_residual,
            contact_check_eta_inv_fm=max(eta,0.03),sampled_inverse_error_margin=margin.margin,
            continuous_UHP_count_certified=margin.count_transfer_certified,
            rigorous_error_bound=false,evaluations=fine.evaluations,production_authorized=false))
    end
    return rows
end

function main()
    s=settings()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_THERMAL_OUTPUT",joinpath(base,"thermal_admissibility_20260906"))
    literature=[joinpath(R.ROOT,"tmp","thermal_spectral_review_20260906","laine_vuorinen.pdf"),
        raw"D:\Desktop\PNJL论文\2026-4-11\2024Pereira PhysRevC.109.025206 New approach to the 3-momentum regularization of the in-medium one- and two-fermion line integrals with application to cross sections in the NJL model.pdf"]
    papers=Dict(p=>R.hashfile(p) for p in literature)
    R.method_contract()
    hashes=R.start_output(output)
    witness=witnesses(bg,s)
    CSV.write(joinpath(output,"spectral_witnesses.csv"),witness)
    println("[thermal] witnesses=$(length(witness)); negative=$(count(r->r.negative_tail_witness,witness))")
    rows,failures=NamedTuple[],NamedTuple[]
    for ch in s.probe_channels,q in s.probe_q,thermal in (bg.vacuum,s.thermal_inv_fm)
        try
            result=probes(bg,ch,q,thermal,s)
            append!(rows,result)
            println("[thermal] $(ch) q=$(q) thermal=$(thermal): direct=$(count(r->r.direct_numerics_passed,result))/$(length(result))")
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(channel=String(ch),q_inv_fm=q,thermal_inv_fm=thermal,reason=sprint(showerror,err)))
            println(stderr,"[thermal-failed] $(last(failures))")
        end
        isempty(rows) || CSV.write(joinpath(output,"continuous_probes.csv"),rows)
        isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
        flush(stdout)
    end
    complete=isempty(failures) && length(rows)==length(s.probe_channels)*length(s.probe_q)*2*3*length(s.etas)
    all(R.hashfile(p)==h for (p,h) in papers) || error("paper changed during audit")
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("background changed")
    R.finish_output(output,hashes,Dict("status"=>complete ? "diagnostic_completed_physical_review_required" : "evaluation_failed",
        "settings"=>s,"background"=>bg,"paper_hashes"=>papers,"all_evaluations_completed"=>complete,
        "witness_numerics_passed"=>all(r.numerics_passed for r in witness),
        "negative_tail_witnesses"=>count(r->r.negative_tail_witness,witness),
        "direct_probe_numerics_passed"=>!isempty(rows) && all(r.direct_numerics_passed for r in rows),
        "full_physical_correlator_interpretation_accepted"=>false,
        "continuous_UHP_stability_certified"=>false,"Mott_Levinson_accepted"=>false,
        "solver_called"=>false,"density_computed"=>false,"production_authorized"=>false,
        "limitations"=>["Positive-density-matrix KMS condition is distinct from signed correlation yield",
            "PNJL is a mean-field extension; Fermi control shares the negative thermal-only tail",
            "Continuous probes use existing geometric partition hints, independent cut ordinates",
            "Embedded quadrature errors and sampled interpolation drifts are not uniform error bounds",
            "No UHP pole discovery or full-continuum zero-count certification is claimed",
            "No clipping, regulator switch, equilibrium solve or new GBU integral"]))
    println("[thermal] output=$(output); completed=$(complete); physical review required")
    complete || error("thermal audit incomplete; evidence retained")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
