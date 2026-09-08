"""Finite-eta GBU integral versus independent discrete-pole plus PV spectrum.

The eta limit is taken at a fixed positive lower endpoint. No pointwise
equality at a bound-state jump or a continuum threshold is required.
"""
module CausalChargedEtaAudit
using CSV, JSON3, SHA
include("audit_causal_charged_bu.jl")
using .CausalChargedBUAudit: ROOT, triplet, hashfile
using Main.RelaxTime.CausalSpectralBubble
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_inverse
using Main.RelaxTime.ChargedRPAProvider: charged_pair_continuum_thresholds
using Main.RelaxTime.BUPhaseGates: certify_gap_roots
using Main.RelaxTime.ChargedPhaseBackend: bu_phase_integral

function main()
    base = joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    output = abspath(get(ENV,"CAUSAL_ETA_OUTPUT",joinpath(base,"causal_eta_integral_review")))
    isdir(output) && !isempty(readdir(output)) && error("refusing to overwrite $(output)")
    input = joinpath(base,"negative_density_phase_fig2_like")
    mp,cp = joinpath(input,"plot_manifest.json"),joinpath(input,"charged_phase_profile_detail.csv")
    m = JSON3.read(read(mp,String))
    retained = collect(CSV.File(cp))
    length(retained) == m.rows || error("input row count mismatch")
    sources = [@__FILE__;CausalChargedBUAudit.SOURCES]
    hashes = Dict(relpath(p,ROOT)=>hashfile(p) for p in sources)
    ih = Dict(p=>hashfile(p) for p in (mp,cp))
    masses,mu = triplet(m.masses_inv_fm),triplet(m.chemical_potentials_inv_fm)
    T,phi,phibar = Float64.((m.thermo.T_inv_fm,m.thermo.Phi,m.thermo.PhiBar))
    q,lower,upper = 1.0,1e-4,40.0
    rows = NamedTuple[]
    for channel in (:pi_plus,:pi_minus,:K_plus,:K_minus)
        spec = charged_rpa_spec(channel)
        f1,f2 = spec.pair
        coupling = Float64(first(filter(r->r.channel==String(channel),retained)).coupling_fm2)
        shift = mu[f1]-mu[f2]
        g = build_spectral_bubble(q,masses[f1],mu[f1],masses[f2],mu[f2],T;
            Phi=phi,PhiBar=phibar,thermal_cutoff_inv_fm=16.0)
        p = build_bubble_dispersion(g;segment_nodes=256)
        inverse(w,eta=0.0) = charged_rpa_inverse(spec,coupling,cauchy_transform(p,complex(w+shift,eta)))
        b = charged_pair_continuum_thresholds(q,masses[f1],masses[f2],mu[f1],mu[f2])
        ld,thr = b.k0_landau_upper_inv_fm,b.k0_threshold_inv_fm
        roots = certify_gap_roots((w,_)->inverse(w),q,[(max(0.0,ld),thr)];physical_sheet=true,real_axis=true)
        roots.passed || error("gap certification failed")
        unitary_w = sort!(unique(vcat(thr .+ 10.0 .^ range(-8,-2;length=60),range(thr+0.011,upper;length=2400))))
        landau_w = exp.(range(log(lower),log(ld-1e-8);length=2400))
        bound = sum((1/expm1(r.omega_inv_fm/T) for r in roots.roots);init=0.0)
        target = bound + bu_phase_integral(unitary_w,[-angle(inverse(w)) for w in unitary_w],T;weight=:gbu)+
            bu_phase_integral(landau_w,[-angle(inverse(w)) for w in landau_w],T;weight=:gbu)
        for eta in (0.02,0.01,0.005,0.0025,0.00125,0.000625)
            ws = sort!(unique(vcat(exp.(range(log(lower),log(0.1);length=160)),
                collect(range(0.101,upper;length=2400)),unitary_w,landau_w,
                [r.omega_inv_fm+x*eta for r in roots.roots for x in range(-32,32;length=513)
                 if lower < r.omega_inv_fm+x*eta < upper])))
            ds = [-angle(inverse(w,eta)) for w in ws]
            maximum(abs.(diff(ds))) < π || error("unresolved finite-eta phase branch")
            value = bu_phase_integral(ws,ds,T;weight=:gbu)
            push!(rows,(channel=String(channel),q_inv_fm=q,eta_inv_fm=eta,k0_min_inv_fm=lower,
                k0_max_inv_fm=upper,nodes=length(ws),pv_split_weight=target,finite_eta_weight=value,
                relative_integral_error=abs(value-target)/abs(target),
                phase_at_zero=-angle(inverse(0.0,eta)),production_authorized=false))
        end
        println("[causal-eta] $(channel): "*string(last(rows)))
    end
    all(p->hashfile(p)==hashes[relpath(p,ROOT)],sources) || error("source changed during audit")
    all(p->hashfile(p)==ih[p],keys(ih)) || error("input changed during audit")
    mkpath(output)
    for source in sources
        target = joinpath(output,"source_snapshot",relpath(source,ROOT))
        mkpath(dirname(target))
        Base.cp(source,target)
    end
    CSV.write(joinpath(output,"eta_integral_comparison.csv"),rows)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,Dict("status"=>"diagnostic_only_not_production","solver_called"=>false,
            "source_hashes"=>hashes,"input_hashes"=>ih,"rows"=>length(rows),"segment_nodes"=>256,
            "order_of_limits"=>"eta to zero at fixed lower, then lower to zero in PV GBU",
            "background"=>m.background,"production_authorized"=>false))
    end
    return rows
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
end
