"""Frozen BQS causal-dispersion and split-BU research diagnostic.

No equilibrium solve, old PV real part, phase fold, sign flip, density clip,
or Landau deletion. Every density remains diagnostic, never production.
"""
module CausalChargedBUAudit
using CSV, JSON3, SHA
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
const SOURCES = [@__FILE__; [joinpath(ROOT,"src","relaxtime",n*".jl") for n in
    ("CausalSpectralBubble","OneLoopIntegrals","ChargedRPAKernel","ChargedRPAProvider","AFieldBuilder","ChargedPhaseBackend","BUPhaseGates")];
    joinpath(ROOT,"src","models","pnjl_physics","QuarkDistribution.jl")]
hashfile(p) = bytes2hex(sha256(read(p)))
const SOURCE_HASHES = Dict(relpath(p,ROOT)=>hashfile(p) for p in SOURCES)
isdefined(Main,:RelaxTime) || Base.include(Main,joinpath(ROOT,"src","relaxtime","RelaxTime.jl"))
using Main.RelaxTime.CausalSpectralBubble
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_inverse
using Main.RelaxTime.ChargedRPAProvider: charged_pair_continuum_thresholds
using Main.RelaxTime.BUPhaseGates: certify_gap_roots
using Main.RelaxTime.ChargedPhaseBackend: bu_phase_integral_parts
using Main.RelaxTime.AFieldBuilder: build_A_triplet
using Main.RelaxTime.GaussLegendre: gauleg

triplet(d) = (u=Float64(d.u),d=Float64(d.d),s=Float64(d.s))
function main()
    base = joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input = abspath(get(ENV,"CAUSAL_BU_INPUT",joinpath(base,"negative_density_phase_fig2_like")))
    output = abspath(get(ENV,"CAUSAL_BU_OUTPUT",joinpath(base,"causal_bu_review_v1")))
    isdir(output) && !isempty(readdir(output)) && error("refusing to overwrite $(output)")
    mp,csv_path = joinpath(input,"plot_manifest.json"),joinpath(input,"charged_phase_profile_detail.csv")
    ih = Dict(p=>hashfile(p) for p in (mp,csv_path))
    m = JSON3.read(read(mp,String))
    retained = collect(CSV.File(csv_path))
    length(retained) == m.rows || error("input row count mismatch")
    masses,mu = triplet(m.masses_inv_fm),triplet(m.chemical_potentials_inv_fm)
    T,phi,phibar = Float64.((m.thermo.T_inv_fm,m.thermo.Phi,m.thermo.PhiBar))
    thermo = (T=T,Φ=phi,Φbar=phibar,ξ=0.0)
    A = build_A_triplet((m=masses,μ=mu),thermo;p_nodes=128,p_max=16.0,use_aniso=false)
    qs = parse.(Float64,split(get(ENV,"CAUSAL_BU_Q","0,0.5,1"),','))
    qn = parse(Int,get(ENV,"CAUSAL_BU_Q_NODES","0"))
    qmax = parse(Float64,get(ENV,"CAUSAL_BU_Q_MAX","4"))
    qweights = Float64[]
    if qn > 0
        qn >= 4 && qmax > 0 || error("q quadrature requires nodes>=4 and qmax>0")
        qs,qweights = gauleg(0.0,qmax,qn)
    end
    mesh = parse(Int,get(ENV,"CAUSAL_BU_MESH","64"))
    np = parse(Int,get(ENV,"CAUSAL_BU_P_NODES","80"))
    nx = parse(Int,get(ENV,"CAUSAL_BU_X_NODES","48"))
    nw = parse(Int,get(ENV,"CAUSAL_BU_W_NODES","600"))
    tc = parse(Float64,get(ENV,"CAUSAL_BU_THERMAL_MAX","16"))
    wmax = parse(Float64,get(ENV,"CAUSAL_BU_W_MAX","16"))
    lows = parse.(Float64,split(get(ENV,"CAUSAL_BU_LOWERS","0.01,0.001,0.0001"),','))
    all(diff(qs) .> 0) && all(qs .>= 0) || error("q grid must be increasing and nonnegative")
    nw >= 32 && all(lows .> 0) && wmax > maximum(lows) || error("invalid omega grid")
    checks,roots,shells,profiles = NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    for channel in (:pi_plus,:pi_minus,:K_plus,:K_minus)
        spec = charged_rpa_spec(channel)
        f1,f2 = spec.pair
        subset = filter(r->r.channel==String(channel) && r.variant=="pv_cut",retained)
        isempty(subset) && error("missing frozen channel $(channel)")
        coupling = Float64(first(subset).coupling_fm2)
        shift = mu[f1]-mu[f2]
        for q in qs
            g = build_spectral_bubble(q,masses[f1],mu[f1],masses[f2],mu[f2],T;
                Phi=phi,PhiBar=phibar,thermal_cutoff_inv_fm=tc,momentum_nodes=np,angle_nodes=nx)
            dispersion = build_bubble_dispersion(g;segment_nodes=mesh)
            inverse(w) = charged_rpa_inverse(spec,coupling,cauchy_transform(dispersion,w+shift))
            bounds = charged_pair_continuum_thresholds(q,masses[f1],masses[f2],mu[f1],mu[f2])
            thr,ld = bounds.k0_threshold_inv_fm,bounds.k0_landau_upper_inv_fm
            thr+0.1 < wmax || error("omega_max must include unitary continuum")
            result = certify_gap_roots((w,_)->inverse(w),q,[(max(0.0,ld),thr)];
                physical_sheet=true,real_axis=true,omega_nodes=128)
            for r in result.roots
                direct = spectral_bubble(g,r.omega_inv_fm;eta_inv_fm=0).value
                push!(roots,(channel=String(channel),q_inv_fm=q,k0_inv_fm=r.omega_inv_fm,
                    lambda_inv_fm=r.omega_inv_fm+shift,dispersion_inverse_residual=r.residual,
                    direct_inverse_residual=abs(charged_rpa_inverse(spec,coupling,direct)),
                    gap_lower_distance=r.distance_to_gap_lower,gap_upper_distance=r.distance_to_gap_upper,
                    physical_mott_certified=false,production_authorized=false))
            end
            for w in (0.0,0.3,2.1,thr+1)
                direct = spectral_bubble(g,w;eta_inv_fm=0.4)
                dispersive = cauchy_transform(dispersion,complex(w+shift,0.4))
                cut = spectral_bubble_cut(g,w)
                push!(checks,(channel=String(channel),q_inv_fm=q,k0_inv_fm=w,
                    direct_real=real(direct.value),direct_imag=imag(direct.value),
                    dispersion_real=real(dispersive),dispersion_imag=imag(dispersive),
                    dispersion_direct_error=abs(direct.value-dispersive),
                    contact_identity_residual=direct.contact_identity_residual,
                    contact_A_difference=q==0 ? g.contact_inv_fm2-A[f1]-A[f2] : NaN,
                    exact_cut_imag=cut.imaginary,static_inverse_real=real(inverse(0.0)),
                    static_inverse_imag=imag(inverse(0.0)),production_authorized=false))
            end
            # These are non-pole intervals. No unwrap determines the bound count.
            unitary_w = sort!(unique(vcat(thr .+ [1e-7,1e-6,1e-5,1e-4,1e-3],
                collect(range(thr+0.002,wmax;length=nw)))))
            unitary_phase = [-angle(inverse(w)) for w in unitary_w]
            maximum(abs.(diff(unitary_phase))) < π || error("unitary branch unresolved")
            threshold_phase = first(unitary_phase)
            for lower in lows
                landau_w = ld > lower ? exp.(range(log(lower),log(ld);length=nw)) : Float64[]
                # Exclude the exact endpoint to avoid signed-zero principal-angle ambiguity.
                isempty(landau_w) || (landau_w[end] = ld-1e-8)
                landau_phase = [-angle(inverse(w)) for w in landau_w]
                isempty(landau_phase) || maximum(abs.(diff(landau_phase))) < π || error("Landau branch unresolved")
                for weight in (:current,:gbu)
                    unitary = bu_phase_integral_parts(unitary_w,unitary_phase,T;weight=weight)
                    landau = isempty(landau_w) ? nothing : bu_phase_integral_parts(landau_w,landau_phase,T;weight=weight)
                    bound = sum((1/expm1(r.omega_inv_fm/T) for r in result.roots);init=0.0)
                    lw = landau === nothing ? 0.0 : landau.derivative
                    prefactor = q^2/(2π^2)
                    ld_count = isempty(landau_phase) ? 0.0 : (last(landau_phase)-first(landau_phase))/π
                    push!(shells,(channel=String(channel),q_inv_fm=q,weight=String(weight),
                        k0_min_inv_fm=lower,k0_max_inv_fm=wmax,gap_count=result.count,gap_passed=result.passed,
                        bound_weight=bound,unitary_weight=unitary.derivative,landau_weight=lw,
                        total_shell_inv_fm2=prefactor*(bound+unitary.derivative+lw),
                        unitary_boundary=unitary.boundary,landau_boundary=landau===nothing ? 0.0 : landau.boundary,
                        threshold_phase_over_pi=threshold_phase/π,
                        levinson_residual=threshold_phase/π-result.count,
                        total_state_balance=result.count+(last(unitary_phase)-threshold_phase)/π+ld_count,
                        infrared_phase=isempty(landau_phase) ? 0.0 : first(landau_phase),
                        upper_phase=last(unitary_phase),production_authorized=false))
                end
                if lower == minimum(lows)
                    for (region,ws,ds) in (("landau",landau_w,landau_phase),("unitary",unitary_w,unitary_phase))
                        for i in eachindex(ws)
                            z = inverse(ws[i])
                            push!(profiles,(channel=String(channel),q_inv_fm=q,region=region,
                                k0_inv_fm=ws[i],phase=ds[i],inverse_real=real(z),inverse_imag=imag(z)))
                        end
                    end
                end
            end
            println("[causal-bu] $(channel) q=$(q) roots=$(result.count) status=$(result.status) rejected=$(result.rejected) static=$(inverse(0.0)) phase_thr/pi=$(threshold_phase/π)")
        end
    end
    all(p->hashfile(p)==SOURCE_HASHES[relpath(p,ROOT)],SOURCES) || error("source changed during audit")
    all(p->hashfile(p)==ih[p],keys(ih)) || error("input changed during audit")
    mkpath(output)
    for p in SOURCES
        target = joinpath(output,"source_snapshot",relpath(p,ROOT))
        mkpath(dirname(target))
        cp(p,target)
    end
    CSV.write(joinpath(output,"causal_checks.csv"),checks)
    CSV.write(joinpath(output,"gap_roots.csv"),isempty(roots) ?
        (channel=String[],q_inv_fm=Float64[],k0_inv_fm=Float64[],production_authorized=Bool[]) : roots)
    CSV.write(joinpath(output,"split_shells.csv"),shells)
    CSV.write(joinpath(output,"phase_profiles.csv"),profiles)
    densities = NamedTuple[]
    if qn > 0
        for channel in (:pi_plus,:pi_minus,:K_plus,:K_minus), weight in (:current,:gbu), lower in lows
            rows = filter(r->r.channel==String(channel) && r.weight==String(weight) && r.k0_min_inv_fm==lower,shells)
            length(rows) == qn || error("q quadrature row mismatch")
            all(rows[i].q_inv_fm==qs[i] for i in eachindex(qs)) || error("q quadrature order mismatch")
            gate = all(r.gap_passed && abs(r.levinson_residual)<0.005 && abs(r.upper_phase)<0.005 for r in rows)
            value = sum(qweights[i]*rows[i].total_shell_inv_fm2 for i in eachindex(qs))
            push!(densities,(channel=String(channel),weight=String(weight),k0_min_inv_fm=lower,
                density_inv_fm3=value,q_nodes=qn,qmax_inv_fm=qmax,conditional_state_gate=gate,
                status=gate ? "conditional_window_integral" : "state_gate_failed",
                production_authorized=false))
        end
        CSV.write(joinpath(output,"density_integrals.csv"),densities)
    end
    record = Dict("status"=>"research_diagnostic_not_production","solver_called"=>false,
        "source_hashes"=>SOURCE_HASHES,"input_hashes"=>ih,"background"=>m.background,
        "git_head"=>readchomp(`git -C $ROOT rev-parse HEAD`),"q_values_inv_fm"=>qs,"q_weights"=>qweights,
        "segment_nodes"=>mesh,"momentum_nodes"=>np,"angle_nodes"=>nx,"omega_nodes"=>nw,
        "thermal_cutoff_inv_fm"=>tc,"vacuum_cutoff_inv_fm"=>Main.Constants_PNJL.Λ_inv_fm,
        "omega_max_inv_fm"=>wmax,"lowers_inv_fm"=>lows,"check_rows"=>length(checks),
        "root_rows"=>length(roots),"shell_rows"=>length(shells),"profile_rows"=>length(profiles),
        "phase"=>"-arg inverse; no fold, shift, unwrap or deletion; poles counted independently",
        "bose_coordinate"=>"external k0 with mu=0","production_authorized"=>false)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,record)
    end
    println("[causal-bu] output=$(output)")
    return record
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
end
