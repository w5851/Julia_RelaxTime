"""B2/B3 finite probes only: no new background, pole count, or meson density."""
module CausalGBUObservableClosureAudit
include("causal_gbu_observable_closure.jl")
const O=CausalGBUObservableClosure
const A=O.A
const R=O.R
using CSV,JSON3

function grid(bg,ch,q,cutoff)
    i,j=R.charged_rpa_spec(ch).pair
    return R.build_spectral_bubble(q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=cutoff,momentum_nodes=128,angle_nodes=64)
end

function real_probes(bg,ch,q,thermal)
    i,j=R.charged_rpa_spec(ch).pair
    shift=bg.mu[i]-bg.mu[j]
    low,threshold=hypot(q,bg.m[i]-bg.m[j]),hypot(q,bg.m[i]+bg.m[j])
    top=hypot(bg.vacuum,bg.m[i])+hypot(bg.vacuum,bg.m[j])
    lower=max(0.,shift)
    low>lower && top>threshold || error("representative windows unavailable")
    centers=(landau=lower+0.63*(low-lower),gap=(low+threshold)/2,
        pair=(threshold+top)/2,negative_thermal_tail=top+0.25)
    edges=A.cut_panels(grid(bg,ch,q,thermal))
    rows=NamedTuple[]
    for (region,lambda) in pairs(centers)
        vals=map(1:2) do k
            rho=x->A.direct_cut(bg,ch,q,x,thermal;nodes=k==1 ? 32 : 64).imaginary
            A.continuous_cauchy(rho,edges,complex(lambda);nodes=k==1 ? 8 : 12,
                atol=k==1 ? 1e-7 : 1e-8)
        end
        coarse,fine=vals
        p=fine.value
        parts=O.scalar_gbu_parts(p,0.,bg.coupling[ch])
        cut=A.direct_cut(bg,ch,q,lambda,thermal;nodes=64)
        f=1-4bg.coupling[ch]*p
        row=(channel=String(ch),q_inv_fm=q,region=String(region),
            lambda_inv_fm=lambda,k0_inv_fm=lambda-shift,thermal_cutoff_inv_fm=thermal,
            pi_real_inv_fm2=real(p),pi_imaginary_inv_fm2=imag(p),
            direct_cut_inv_fm2=cut.imaginary,cut_error_inv_fm2=abs(imag(p)-cut.imaginary),
            node_change_inv_fm2=abs(fine.value-coarse.value),
            quadrature_error_estimate_inv_fm2=fine.error_estimate,
            numerics_passed=coarse.converged && fine.converged &&
                abs(fine.value-coarse.value)<1e-6 && abs(imag(p)-cut.imaginary)<1e-9,
            inverse_real=real(f),inverse_imaginary=imag(f),
            principal_phase=parts.phase,gbu_local_weight=parts.weight,
            selfenergy_correction=parts.selfenergy_correction,
            weight_identity_error=parts.weight_error,optical_identity_error_fm2=parts.optical_error,
            algebra_passed=max(parts.weight_error,parts.optical_error)<1e-12,
            negative_spectrum_retained=imag(p)<0,global_phase_certified=false,
            derivative_evaluated=false,density_computed=false,production_authorized=false)
        push!(rows,row)
    end
    return rows
end

function log_probes(bg,ch,q,thermal)
    i,j=R.charged_rpa_spec(ch).pair
    phi,bar=bg.Phi,bg.PhiBar
    args=(q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T)
    zs=ComplexF64[0.8im,2+0.8im,6+1.1im] # external k0, not lambda
    lambdas=zs .+ (bg.mu[i]-bg.mu[j])
    pieces=[O.P.cylindrical_reference(args...,L,lambdas;Phi=phi,PhiBar=bar,
        component=component,nz=128,ny=128) for (L,component) in
        ((bg.vacuum,:vacuum),(thermal,:full),(thermal,:vacuum))]
    rows=NamedTuple[]
    for (k,k0) in enumerate(zs)
        z=lambdas[k]
        ps=[p.pi_p[k] for p in pieces]
        r=O.subtraction_logs(ps...,bg.coupling[ch])
        push!(rows,(channel=String(ch),q_inv_fm=q,k0_real_inv_fm=real(k0),
            eta_inv_fm=imag(k0),lambda_real_inv_fm=real(z),
            combined_log_real=real(r.combined),combined_log_imaginary=imag(r.combined),
            piecewise_log_real=real(r.separately_resummed),
            piecewise_log_imaginary=imag(r.separately_resummed),
            noncommutation_abs=abs(r.difference),
            physical_equivalence_claimed=false,production_authorized=false))
    end
    return rows
end

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_OBSERVABLE_OUTPUT",joinpath(base,"direction_b_observable_20260907"))
    literature=[joinpath(R.ROOT,"tmp","direction_b_gbu_review_20260907",p*".pdf")
        for p in ("2512.03876","1612.09556")]
    append!(literature,[joinpath(R.ROOT,"work","charged_phase_low_energy_20260905",p*".pdf")
        for p in ("1305.3907","1912.13162")])
    paper_hashes=Dict(p=>R.hashfile(p) for p in literature)
    R.method_contract()
    hashes=R.start_output(output)
    rows,logs,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    for ch in (:pi_plus,:K_plus),q in (0.,1.4,3.2)
        try
            append!(rows,real_probes(bg,ch,q,10.))
            append!(logs,log_probes(bg,ch,q,10.))
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(channel=String(ch),q_inv_fm=q,reason=sprint(showerror,err)))
            println(stderr,"[observable-failed] $(last(failures))")
        end
        isempty(rows) || CSV.write(joinpath(output,"real_axis_algebra.csv"),rows)
        isempty(logs) || CSV.write(joinpath(output,"subtraction_logs.csv"),logs)
        isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
        println("[observable] $(ch) q=$(q): real rows=$(length(rows)); log rows=$(length(logs))")
        flush(stdout)
    end
    all(R.hashfile(p)==h for (p,h) in paper_hashes) || error("paper drift")
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    complete=isempty(failures) && length(rows)==24 && length(logs)==18
    R.finish_output(output,hashes,Dict("status"=>complete ? "B2_B3_local_algebra_review" : "evaluation_failed",
        "background"=>bg,"paper_hashes"=>paper_hashes,"real_rows"=>length(rows),"log_rows"=>length(logs),
        "all_evaluations_completed"=>complete,"algebra_passed"=>complete && all(r.algebra_passed for r in rows),
        "numerics_passed"=>complete && all(r.numerics_passed for r in rows),
        "selfenergy_normalization"=>"Sigma_M=2Pi; D0=2K; F=1-4KPi",
        "cut_nodes"=>[32,64],"quadrature_nodes"=>[8,12],"quadrature_atol_inv_fm2"=>[1e-7,1e-8],
        "node_tolerance_inv_fm2"=>1e-6,"cut_tolerance_inv_fm2"=>1e-9,"algebra_tolerance"=>1e-12,
        "full_stationarity_certified"=>false,"global_phase_certified"=>false,
        "continuous_UHP_count_certified"=>false,"GBU_observable_derived_from_stationarity"=>false,
        "solver_called"=>false,"meson_density_computed"=>false,
        "limitations"=>["Principal local phase only; not a density branch",
            "No positivity or whole-contour certification follows from point probes",
            "No differentiation of the physical profile is performed in this runner",
            "GBU weight is not Gaussian logdet or a piecewise determinant subtraction",
            "No change to upstream, physical prescription, density provider or production defaults"]))
    complete || error("observable audit incomplete; all evidence retained")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
