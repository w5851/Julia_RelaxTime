#!/usr/bin/env julia
# Optional oracle, confined to the benchmark environment (ADR-0006).
const THRESHOLD_BENCH_ROOT=normpath(joinpath(@__DIR__,"..",".."))
push!(LOAD_PATH,joinpath(THRESHOLD_BENCH_ROOT,"benchmark"))
using QuadGK, CSV
include(joinpath(THRESHOLD_BENCH_ROOT,"scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"))
using .ThresholdCacheCandidates, .ThresholdUpgradeCommon

function ordered_loop_reference(plan,s)
    q2=(s-4plan.m^2)/4; L=plan.cutoff
    f(p)=pair_occupation(plan,p*p)
    if q2==0
        re,err=quadgk(p->2f(p),0.0,L;rtol=1e-11,atol=2e-13)
        return complex(re,0.0),err
    elseif q2<0
        kappa=sqrt(-q2)
        # Independent direct integral, resolving the narrow below-threshold layer.
        edges=sort!(unique(vcat([0.0,L],[p for p in (kappa,10kappa) if 0<p<L])))
        re,err=quadgk(p->2p*p*f(p)/(p*p+kappa*kappa),edges...;rtol=1e-11,atol=2e-13)
        return complex(re,0.0),err
    else
        q=sqrt(q2)
        q<L || error("reference study is below the upper cutoff")
        fq=f(q)
        hq=q*fq/2
        function regular(p)
            if p==q
                return 0.75fq+q2*ThresholdCacheCandidates.occupation_derivative(plan,q2)
            end
            return (p*p*f(p)/(p+q)-hq)/(p-q)
        end
        re,err=quadgk(regular,0.0,q,L;rtol=1e-11,atol=2e-13)
        return complex(2*(re+hq*(log(L-q)-log(q))),π*q*fq),2err
    end
end

function old_pair_loop(state,flavor,s)
    qp,tp=state.quark_params,state.thermo_params
    m=getproperty(qp.m,flavor); mu=getproperty(qp.μ,flavor)
    lambda=sqrt(s)
    b=Main.RelaxTime.OneLoopIntegrals.B0(lambda,0.0,m,mu,m,mu,tp.T;Φ=tp.Φ,Φbar=tp.Φbar)
    c=tp.ξ==0 ? (0.0,0.0) : Main.RelaxTime.OneLoopIntegralsCorrection.B0_correction(
        lambda,0.0,m,m,mu,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ)
    return complex(b...)+complex(c...)
end

function run_endpoint_audit(inputdir,outdir)
    rows=NamedTuple[]
    for xi in (-0.25,-0.2,0.0,0.36,0.37,0.38)
        state=retained_state(retained_row(inputdir,xi))
        for flavor in (:u,:s)
            qp,tp=state.quark_params,state.thermo_params
            m=getproperty(qp.m,flavor); mu=getproperty(qp.μ,flavor)
            plans=Dict(n=>PairLoopPlan(m,mu,tp;n_nodes=n) for n in (32,64,128))
            s0=4m*m
            for epsilon in (-1e-3,-1e-6,-1e-9,-1e-12,0.0,1e-12,1e-9,1e-6,1e-3)
                s=s0*(1+epsilon)
                reference,oracle_error=ordered_loop_reference(plans[128],s)
                old=complex(NaN,NaN); status="finite"
                try
                    old=old_pair_loop(state,flavor,s)
                    isfinite(old) || (status="nonfinite")
                catch err
                    status=string(typeof(err))
                end
                b32=pair_b0_pv(plans[32],s)
                b64=pair_b0_pv(plans[64],s)
                b128=pair_b0_pv(plans[128],s)
                push!(rows,(xi=xi,flavor=String(flavor),epsilon=epsilon,s=s,
                    reference_re=real(reference),reference_im=imag(reference),
                    oracle_error=oracle_error,old_status=status,old_re=real(old),old_im=imag(old),
                    old_absolute_error=abs(old-reference),
                    n32_error=abs(b32-reference),n64_error=abs(b64-reference),n128_error=abs(b128-reference)))
            end
        end
        write_rows(joinpath(outdir,"endpoint_audit.csv"),rows)
        println("ENDPOINT xi=",xi);flush(stdout)
    end
    maximum(r.n64_error for r in rows)<1e-8 ||
        error("Paired k0 candidate did not meet the endpoint oracle tolerance")
end

function main(inputdir,outdir)
    start_run(outdir,inputdir,@__FILE__,(
        purpose="independent endpoint oracle for equal-mass k0, same linear-xi loop",
        oracle="QuadGK in benchmark environment; alternate p-Cauchy subtraction",
        xis=[-0.25,-0.2,0.0,0.36,0.37,0.38],orders=[32,64,128],
        tolerance_absolute=1e-8),extra_sources=[
            joinpath(THRESHOLD_BENCH_ROOT,"scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"),
            joinpath(THRESHOLD_BENCH_ROOT,"benchmark","Project.toml"),
            joinpath(THRESHOLD_BENCH_ROOT,"benchmark","Manifest.toml")])
    run_endpoint_audit(joinpath(outdir,"inputs"),outdir)
    finish_run(outdir;endpoint_comparisons=108)
    println("COMPLETE ",outdir)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: benchmark_threshold_endpoint.jl INPUT_DIR NEW_OUTPUT_DIR")
    main(abspath(ARGS[1]),abspath(ARGS[2]))
end
