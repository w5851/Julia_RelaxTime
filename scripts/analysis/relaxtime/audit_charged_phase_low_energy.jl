"""Frozen-input infrared and spectral-consistency diagnostic, never production.

No equilibrium solve, phase fold, sign flip, or Landau deletion is performed.
The independent cut oracle is not spliced into the retained PV real part.
"""
module ChargedPhaseLowEnergyAudit
using CSV, JSON3, SHA
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
const SOURCES = [@__FILE__; [joinpath(ROOT,"src","relaxtime",name*".jl") for name in
    ("OneLoopIntegrals","ChargedRPAProvider","ChargedPhaseBackend","BUPhaseGates","PhaseNormalization")]]
const SOURCE_HASHES = Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCES)
if !isdefined(Main,:RelaxTime)
    Base.include(Main,joinpath(ROOT,"src","relaxtime","RelaxTime.jl"))
end
using Main.RelaxTime.OneLoopIntegrals: B0_pv_cut, B0_retarded, B0_spectral_cut
using Main.RelaxTime.ChargedRPAProvider: charged_polarization
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_inverse
using Main.RelaxTime.ChargedPhaseBackend: bu_phase_integral_parts
using Main.RelaxTime.AFieldBuilder: build_A_triplet

triplet(x) = (u=Float64(x.u),d=Float64(x.d),s=Float64(x.s))
filehash(p) = bytes2hex(sha256(read(p)))

function main()
    base = joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input = abspath(get(ENV,"CHARGED_IR_INPUT",joinpath(base,"negative_density_phase_fig2_like")))
    output = abspath(get(ENV,"CHARGED_IR_OUTPUT",joinpath(base,"low_energy_cut_audit_20260905")))
    isdir(output) && !isempty(readdir(output)) && error("refusing to overwrite $(output)")
    manifest_path, csv_path = joinpath(input,"plot_manifest.json"),joinpath(input,"charged_phase_profile_detail.csv")
    input_hashes = Dict(p=>filehash(p) for p in (manifest_path,csv_path))
    manifest = JSON3.read(read(manifest_path,String))
    retained = collect(CSV.File(csv_path))
    length(retained)==manifest.rows || error("retained row count mismatch")
    masses,mu = triplet(manifest.masses_inv_fm),triplet(manifest.chemical_potentials_inv_fm)
    thermo = (T=Float64(manifest.thermo.T_inv_fm),Φ=Float64(manifest.thermo.Phi),
              Φbar=Float64(manifest.thermo.PhiBar),ξ=0.0)
    A = build_A_triplet((m=masses,μ=mu),thermo;p_nodes=64,p_max=16.0,use_aniso=false)
    nodes = parse(Int,get(ENV,"CHARGED_IR_NODES","257"))
    nodes >= 32 || error("CHARGED_IR_NODES must be at least 32")
    probe_rows, integral_rows, vacuum_rows = NamedTuple[],NamedTuple[],NamedTuple[]
    parity = 0.0
    for lambda in (0.1,0.5,0.9,3.0)
        old = B0_pv_cut(lambda,1.0,1.0,0.0,1.0,0.0,0.001)
        spectral = B0_spectral_cut(lambda,1.0,1.0,0.0,1.0,0.0,0.001)
        finite_eta = B0_retarded(lambda,1.0,1.0,0.0,1.0,0.0,0.001;eta_inv_fm=1e-4,energy_nodes=512)
        push!(vacuum_rows,(lambda_inv_fm=lambda,q_inv_fm=1.0,mass_inv_fm=1.0,T_inv_fm=0.001,
            region=lambda<1 ? "spacelike" : "unitary",old_pv_imag=imag(old),
            old_finite_eta_imag=imag(finite_eta),spectral_pair=spectral.pair,spectral_landau=spectral.landau,
            spectral_imag=spectral.imaginary,cutoff_regulator="both_line_momenta",
            eta_probe_inv_fm=1e-4,eta_convergence_certified=false))
    end
    for channel in (:pi_plus,:pi_minus,:K_plus,:K_minus)
        spec = charged_rpa_spec(channel)
        rows = filter(r->r.channel==String(channel) && r.variant=="pv_cut",retained)
        isempty(rows) && error("missing frozen channel $(channel)")
        coupling = Float64(first(rows).coupling_fm2)
        inverse(w,q) = charged_rpa_inverse(spec,coupling,
            charged_polarization(spec,w,q,masses,mu,thermo,A;prescription=:ordered_pv_cut).value)
        for i in unique(round.(Int,range(1,length(rows);length=24)))
            r = rows[i]
            parity = max(parity,abs(inverse(r.omega_inv_fm,r.q_inv_fm)-complex(r.inverse_real,r.inverse_imag)))
        end
        parity <= 1e-10 || error("frozen inverse parity failed")
        f1,f2 = spec.pair
        for q in (0.0,0.5,1.0)
            z0 = inverse(0.0,q)
            delta0 = -angle(z0)
            slope0 = (-angle(inverse(1e-6,q))-delta0)/1e-6
            for w in (-0.01,-0.001,0.0,1e-5,1e-4,1e-3,0.01,0.05)
                z = inverse(w,q)
                lambda = w+mu[f1]-mu[f2]
                spectral = B0_spectral_cut(lambda,q,masses[f1],mu[f1],masses[f2],mu[f2],thermo.T;Φ=thermo.Φ,Φbar=thermo.Φbar)
                old = B0_pv_cut(lambda,q,masses[f1],mu[f1],masses[f2],mu[f2],thermo.T;Φ=thermo.Φ,Φbar=thermo.Φbar)
                push!(probe_rows,(channel=String(channel),q_inv_fm=q,k0_inv_fm=w,lambda_inv_fm=lambda,
                    inverse_real=real(z),inverse_imag=imag(z),raw_phase=-angle(z),
                    old_b0_imag=imag(old),spectral_b0_imag=spectral.imaginary,
                    spectral_pair=spectral.pair,spectral_landau=spectral.landau,
                    zero_frequency_check=w==0.0,full_propagator_from_oracle=false))
            end
            for lower in (0.01,0.001,0.0001,0.00001)
                w = exp.(range(log(lower),log(0.05);length=nodes))
                phase = [-angle(inverse(x,q)) for x in w]
                maximum(abs.(diff(phase))) < π || error("unexpected principal branch crossing in infrared window")
                for weight in (:current,:gbu)
                    parts = bu_phase_integral_parts(w,phase,thermo.T;weight=weight)
                    derivative0 = weight===:current ? slope0 : 2sin(delta0)^2*slope0
                    push!(integral_rows,(channel=String(channel),q_inv_fm=q,weight=String(weight),
                        k0_min_inv_fm=lower,k0_max_inv_fm=0.05,nodes=nodes,
                        raw_phase_at_zero=delta0,raw_phase_slope_at_zero_fm=slope0,
                        predicted_log_coefficient=thermo.T*derivative0/π,
                        derivative_weight=parts.derivative,bulk_weight=parts.bulk,boundary_weight=parts.boundary,
                        lower_boundary=parts.lower_boundary,upper_boundary=parts.upper_boundary,
                        identity_residual=parts.identity_residual,
                        coordinate="external_k0_mu_zero",physical_cut_certified=false,production_authorized=false))
                end
            end
        end
    end
    all(p->filehash(p)==SOURCE_HASHES[relpath(p,ROOT)],SOURCES) || error("source changed during audit")
    all(p->filehash(p)==input_hashes[p],keys(input_hashes)) || error("input changed during audit")
    mkpath(output)
    for p in SOURCES
        target = joinpath(output,"source_snapshot",relpath(p,ROOT))
        mkpath(dirname(target))
        cp(p,target)
    end
    CSV.write(joinpath(output,"vacuum_cut_checks.csv"),vacuum_rows)
    CSV.write(joinpath(output,"frozen_low_energy_probes.csv"),probe_rows)
    CSV.write(joinpath(output,"infrared_boundary_terms.csv"),integral_rows)
    record = Dict("status"=>"diagnostic_only_not_production","solver_called"=>false,
        "git_head"=>readchomp(`git -C $ROOT rev-parse HEAD`),"input_hashes"=>input_hashes,
        "source_hashes"=>SOURCE_HASHES,"background"=>manifest.background,
        "A_nodes"=>64,"A_pmax_inv_fm"=>16.0,"inverse_parity_max_abs"=>parity,
        "vacuum_rows"=>length(vacuum_rows),"probe_rows"=>length(probe_rows),"integral_rows"=>length(integral_rows),
        "infrared_nodes"=>nodes,"cut_oracle_boundary"=>"independent two-line regulator; imaginary part only, no hybrid propagator",
        "physical_closure"=>"failed vacuum-spacelike and zero-k0 checks in retained provider; no density promotion",
        "sources"=>["https://arxiv.org/abs/1305.3907","https://arxiv.org/abs/1912.13162","https://arxiv.org/abs/2301.09882"])
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,record)
    end
    println("[charged-ir] output=$(output), inverse_parity=$(parity)")
    return record
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
end
