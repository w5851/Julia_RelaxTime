"""Solver-free, frozen-input analytic-gap and split-spectrum diagnostic.

Consumes the retained Fig2 CSV/manifest without modifying them. Rebuilds A
from explicit settings, verifies inverse-value parity against retained rows,
then evaluates independent gap roots and disjoint Landau/unitary profiles.
No production acceptance follows from the conditional one-cut Levinson test.
"""
module ChargedGapSpectrumAudit

using CSV, JSON3, SHA
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SOURCE_PATHS = [@__FILE__,joinpath(PROJECT_ROOT,"src","relaxtime","ChargedRPAProvider.jl"),
    joinpath(PROJECT_ROOT,"src","relaxtime","ChargedPhaseBackend.jl"),
    joinpath(PROJECT_ROOT,"src","relaxtime","BUPhaseGates.jl"),
    joinpath(PROJECT_ROOT,"src","relaxtime","OneLoopIntegrals.jl")]
const SOURCE_HASHES = Dict(relpath(p,PROJECT_ROOT)=>bytes2hex(sha256(read(p))) for p in SOURCE_PATHS)
if !isdefined(Main, :RelaxTime)
    Base.include(Main, joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))
end
using Main.RelaxTime.AFieldBuilder: build_A_triplet
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_inverse
using Main.RelaxTime.ChargedRPAProvider: charged_polarization, charged_pair_continuum_thresholds
using Main.RelaxTime.BUPhaseGates: continue_gap_roots, levinson_phase_gate, mott_phase_gate
using Main.RelaxTime.ChargedPhaseBackend: strict_phase_profile, bu_phase_integral, split_bu_shell

const INPUT_DEFAULT = joinpath(PROJECT_ROOT,"data","outputs","results","relaxtime","analysis",
    "charged_rpa_phase_backend","negative_density_phase_fig2_like")

_triplet(d) = (u=Float64(d.u),d=Float64(d.d),s=Float64(d.s))
_hash(path) = bytes2hex(sha256(read(path)))

function _write_root_rows(path, rows)
    # CSV cannot infer columns from an empty Vector{NamedTuple}. A no-root
    # result must still have a readable schema, not a missing output file.
    empty_columns = (channel=String[],q_inv_fm=Float64[],track_id=Int[],event=String[],
        omega_inv_fm=Float64[],lambda_inv_fm=Float64[],residual=Float64[],slope=Float64[],
        landau_distance_inv_fm=Float64[],search_lower_distance_inv_fm=Float64[],
        unitary_distance_inv_fm=Float64[],count_scope=String[],status=String[])
    return CSV.write(path,isempty(rows) ? empty_columns : rows)
end

function compare_mott(before_dir, after_dir, output)
    isdir(output) && !isempty(readdir(output)) && error("refusing to overwrite $(output)")
    manifests = [JSON3.read(read(joinpath(d,"manifest.json"),String)) for d in (before_dir,after_dir)]
    before_manifest, after_manifest = manifests
    before_manifest.background.T_MeV < after_manifest.background.T_MeV || error("Mott temperatures must increase")
    for key in (:muB_MeV,:Q_over_B,:rhoS_target)
        before_manifest.background[key] == after_manifest.background[key] || error("Mott background mismatch: $(key)")
    end
    for key in (:q_values_inv_fm,:A_nodes,:A_pmax_inv_fm,:omega_max_inv_fm,:bose_coordinate)
        before_manifest[key] == after_manifest[key] || error("Mott configuration mismatch: $(key)")
    end
    tables = [filter(r->r.weight=="current",collect(CSV.File(joinpath(d,"split_shells.csv")))) for d in (before_dir,after_dir)]
    length(tables[1]) == length(tables[2]) || error("Mott row count mismatch")
    rows = NamedTuple[]
    for b in tables[1]
        matches = filter(a->a.channel==b.channel && a.q_inv_fm==b.q_inv_fm,tables[2])
        length(matches)==1 || error("Mott pair must be unique")
        a = only(matches)
        input(r) = (passed=r.conditional_levinson_passed,bound_state_count=r.count_gap,
                    threshold_phase=r.threshold_phase_over_pi*π)
        gate = mott_phase_gate(input(b),input(a))
        push!(rows,(channel=b.channel,q_inv_fm=b.q_inv_fm,
            T_before_MeV=before_manifest.background.T_MeV,T_after_MeV=after_manifest.background.T_MeV,
            gap_count_before=b.count_gap,gap_count_after=a.count_gap,
            bound_count_drop=gate.bound_state_count_drop,
            threshold_phase_drop_over_pi=gate.threshold_phase_drop/π,
            phase_residual=gate.phase_residual,conditional_mott_passed=gate.passed,
            physical_mott_certified=false,transition_location="not_located",
            production_authorized=false))
    end
    mkpath(output)
    CSV.write(joinpath(output,"conditional_mott_pairs.csv"),rows)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,Dict("status"=>"diagnostic_only_not_production","solver_called"=>false,
            "before"=>before_dir,"after"=>after_dir,
            "input_hashes"=>Dict(joinpath(d,f)=>_hash(joinpath(d,f)) for d in (before_dir,after_dir) for f in ("manifest.json","split_shells.csv")),
            "source_sha256"=>_hash(@__FILE__),"pairs"=>length(rows),
            "conditional_passed"=>count(r->r.conditional_mott_passed,rows),
            "boundary"=>"two endpoints only; full cut count, root entry event and Mott temperature not certified"))
    end
    println("[gap-mott] conditional=$(count(r->r.conditional_mott_passed,rows))/$(length(rows)); physical_mott_certified=false")
    return rows
end

function main()
    input = abspath(get(ENV,"CHARGED_GAP_INPUT_DIR",INPUT_DEFAULT))
    output = abspath(get(ENV,"CHARGED_GAP_OUTPUT_DIR",joinpath(dirname(input),"gap_spectrum_audit_20260905")))
    isdir(output) && !isempty(readdir(output)) && error("refusing to overwrite $(output)")
    manifest_path = joinpath(input,"plot_manifest.json")
    csv_path = joinpath(input,"charged_phase_profile_detail.csv")
    manifest = JSON3.read(read(manifest_path,String))
    source = collect(CSV.File(csv_path))
    length(source) == manifest.rows || error("input row count does not match manifest")
    masses, mu = _triplet(manifest.masses_inv_fm), _triplet(manifest.chemical_potentials_inv_fm)
    thermo = (T=Float64(manifest.thermo.T_inv_fm),Φ=Float64(manifest.thermo.Phi),
              Φbar=Float64(manifest.thermo.PhiBar),ξ=0.0)
    a_nodes = parse(Int,get(ENV,"CHARGED_GAP_A_NODES","64"))
    a_max = parse(Float64,get(ENV,"CHARGED_GAP_A_PMAX","16.0"))
    A = build_A_triplet((m=masses,μ=mu),thermo;p_nodes=a_nodes,p_max=a_max,use_aniso=false)
    qvalues = parse.(Float64,split(get(ENV,"CHARGED_GAP_Q_VALUES","0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1"),','))
    nodes = parse(Int,get(ENV,"CHARGED_GAP_ROOT_NODES","128"))
    omega_nodes = parse(Int,get(ENV,"CHARGED_GAP_OMEGA_NODES","360"))
    omega_max = parse(Float64,get(ENV,"CHARGED_GAP_OMEGA_MAX","8.0"))
    omega_min = parse(Float64,get(ENV,"CHARGED_GAP_OMEGA_MIN","0.05"))
    omega_nodes >= 16 && omega_max > omega_min > 0 || error("invalid omega window/nodes")
    root_rows, shell_rows, profile_rows, event_rows = NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    parity_max = 0.0
    for channel in (:pi_plus,:pi_minus,:K_plus,:K_minus)
        retained = filter(r->r.channel==String(channel) && r.variant=="pv_cut",source)
        isempty(retained) && error("missing input channel $(channel)")
        spec = charged_rpa_spec(channel)
        coupling = Float64(first(retained).coupling_fm2)
        inverse(w,q) = charged_rpa_inverse(spec,coupling,
            charged_polarization(spec,w,q,masses,mu,thermo,A;prescription=:ordered_pv_cut).value)
        # Check actual inverse values, not merely manifest labels or gap residuals.
        for i in unique(round.(Int,range(1,length(retained);length=24)))
            r = retained[i]
            parity_max = max(parity_max,abs(inverse(r.omega_inv_fm,r.q_inv_fm)-complex(r.inverse_real,r.inverse_imag)))
        end
        parity_max <= 1e-10 || error("frozen-input inverse parity failed: $(parity_max)")
        thresholds(q) = charged_pair_continuum_thresholds(q,masses[spec.pair[1]],masses[spec.pair[2]],
                                                         mu[spec.pair[1]],mu[spec.pair[2]])
        gaps(q) = [(max(omega_min,thresholds(q).analytic_gap_inv_fm[1]),thresholds(q).analytic_gap_inv_fm[2])]
        continuation = continue_gap_roots(inverse,qvalues,gaps;physical_sheet=true,real_axis=true,
                                          omega_nodes=nodes,max_motion=0.25)
        for result in continuation
            q, bounds = result.q, thresholds(result.q)
            thr = bounds.k0_threshold_inv_fm
            thr+1e-4 < omega_max || error("unitary threshold outside omega window")
            for root in result.roots
                push!(root_rows,(channel=String(channel),q_inv_fm=q,track_id=root.track_id,
                    event=String(root.event),omega_inv_fm=root.omega_inv_fm,
                    lambda_inv_fm=root.omega_inv_fm+bounds.chemical_potential_shift_inv_fm,
                    residual=root.residual,slope=root.slope,
                    landau_distance_inv_fm=root.omega_inv_fm-bounds.k0_landau_upper_inv_fm,
                    search_lower_distance_inv_fm=root.distance_to_gap_lower,
                    unitary_distance_inv_fm=root.distance_to_gap_upper,
                    count_scope=String(result.count_scope),status=String(result.status)))
            end
            for event in result.lost_tracks
                push!(event_rows,(channel=String(channel),q_inv_fm=q,track_id=event.track_id,event=String(event.event)))
            end
            w = sort(unique(vcat(thr .+ [1e-5,1e-4,1e-3],collect(range(thr+0.01,omega_max;length=omega_nodes)))))
            profile = strict_phase_profile(w,ComplexF64[inverse(x,q) for x in w])
            gate = levinson_phase_gate(w,profile.raw_phase,first(w);bound_state_count=result.count)
            unitary = [(omega=w,phase=profile.unwrapped_phase)]
            for i in eachindex(w)
                push!(profile_rows,(channel=String(channel),q_inv_fm=q,region="unitary",omega_inv_fm=w[i],
                    phase=profile.unwrapped_phase[i],inverse_real=real(profile.inverse_values[i]),inverse_imag=imag(profile.inverse_values[i])))
            end
            # Landau envelope includes cutoff-created gaps. It is never used for
            # bound-state counting; zero Im on this region is not certification.
            landau_hi = bounds.k0_landau_upper_inv_fm
            landau_w, landau_phase = Float64[],Float64[]
            if landau_hi > omega_min + 1e-4
                landau_w = collect(range(omega_min,landau_hi-1e-5;length=omega_nodes))
                lp = strict_phase_profile(landau_w,ComplexF64[inverse(x,q) for x in landau_w])
                landau_phase = lp.unwrapped_phase
                for i in eachindex(landau_w)
                    push!(profile_rows,(channel=String(channel),q_inv_fm=q,region="landau_envelope",omega_inv_fm=landau_w[i],
                        phase=landau_phase[i],inverse_real=real(lp.inverse_values[i]),inverse_imag=imag(lp.inverse_values[i])))
                end
            end
            for weight in (:current,:gbu)
                # External Matsubara k0 and kinetic lambda are related by the
                # flavor shift. Both coordinate representations use g(k0;0).
                split = result.passed ? split_bu_shell(result,unitary,q,thermo.T;weight=weight) : nothing
                landau = isempty(landau_w) ? 0.0 :
                    q^2/(2π^2)*bu_phase_integral(landau_w,landau_phase,thermo.T;weight=weight)
                historical_mu = bounds.chemical_potential_shift_inv_fm
                legacy_split = result.passed && minimum([first(w);[r.omega_inv_fm for r in result.roots]]) > historical_mu ?
                    split_bu_shell(result,unitary,q,thermo.T;μ=historical_mu,weight=weight) : nothing
                legacy_landau = isempty(landau_w) ? 0.0 : (first(landau_w) > historical_mu ?
                    q^2/(2π^2)*bu_phase_integral(landau_w,landau_phase,thermo.T;μ=historical_mu,weight=weight) : missing)
                push!(shell_rows,(channel=String(channel),q_inv_fm=q,weight=String(weight),
                    count_gap=result.count,root_status=String(result.status),
                    unitary_threshold_inv_fm=thr,landau_upper_inv_fm=landau_hi,
                    bound_shell_inv_fm2=split===nothing ? NaN : split.bound_shell_inv_fm2,
                    unitary_shell_inv_fm2=split===nothing ? NaN : split.continuum_shell_inv_fm2,
                    landau_window_shell_inv_fm2=landau,
                    partial_sum_inv_fm2=split===nothing ? NaN : split.total_shell_inv_fm2+landau,
                    legacy_bose_bound_plus_unitary_inv_fm2=legacy_split===nothing ? NaN : legacy_split.total_shell_inv_fm2,
                    legacy_bose_landau_inv_fm2=legacy_landau,
                    legacy_bose_partial_sum_inv_fm2=legacy_split===nothing ? missing : legacy_split.total_shell_inv_fm2+legacy_landau,
                    threshold_phase_over_pi=gate.threshold_phase/π,
                    levinson_residual=gate.levinson_residual,
                    conditional_levinson_passed=result.passed && gate.passed,
                    tail_span=profile.tail_span,tail_before=profile.high_energy_phase_before_anchor,
                    full_state_count_certified=false,production_authorized=false))
            end
            println("[gap-spectrum] $(channel) q=$(q) roots=$(result.count) $(result.status) Levinson=$(gate.passed) residual=$(gate.levinson_residual)")
        end
    end
    all(p -> _hash(p)==SOURCE_HASHES[relpath(p,PROJECT_ROOT)],SOURCE_PATHS) ||
        error("source changed during audit; refusing to publish mismatched provenance")
    mkpath(output)
    for p in SOURCE_PATHS
        destination = joinpath(output,"source_snapshot",relpath(p,PROJECT_ROOT))
        mkpath(dirname(destination))
        cp(p,destination)
    end
    _write_root_rows(joinpath(output,"gap_roots.csv"),root_rows)
    CSV.write(joinpath(output,"split_shells.csv"),shell_rows)
    CSV.write(joinpath(output,"continuum_profiles.csv"),profile_rows)
    isempty(event_rows) || CSV.write(joinpath(output,"root_events.csv"),event_rows)
    record = Dict("status"=>"diagnostic_only_not_production","solver_called"=>false,
        "source_manifest"=>manifest_path,"source_manifest_sha256"=>_hash(manifest_path),
        "source_csv_sha256"=>_hash(csv_path),"source_rows"=>length(source),
        "git_head"=>readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`),
        "source_hashes"=>SOURCE_HASHES,
        "source_snapshot"=>"source_snapshot/ (bytes as evaluated; never staged)",
        "background"=>manifest.background,"masses_inv_fm"=>masses,"mu_inv_fm"=>mu,"thermo"=>thermo,
        "A_nodes"=>a_nodes,"A_pmax_inv_fm"=>a_max,"A_values_inv_fm2"=>A,
        "A_setting_provenance"=>"explicit rebuild; old manifest omitted A settings; inverse parity required",
        "inverse_parity_max_abs"=>parity_max,"q_values_inv_fm"=>qvalues,"root_nodes"=>nodes,
        "omega_nodes_per_segment"=>omega_nodes,"omega_min_inv_fm"=>omega_min,"omega_max_inv_fm"=>omega_max,
        "gap_scope"=>"conservative uncut envelope only; cutoff-created gaps not certified",
        "bose_coordinate"=>"external_k0: g(k0;mu=0); lambda=k0+mu1-mu2",
        "legacy_bose_column"=>"g(k0-mu_M) sensitivity only, no coordinate closure claim",
        "events"=>event_rows,"root_rows"=>length(root_rows),"shell_rows"=>length(shell_rows),
        "Mott_status"=>"not_tested_single_background","joint_convergence_status"=>"deferred_until_physical_closure")
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,record)
    end
    println("[gap-spectrum] output=$(output) inverse_parity=$(parity_max)")
    return record
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    if get(ENV,"CHARGED_GAP_MOTT_ONLY","false") == "true"
        compare_mott(ENV["CHARGED_GAP_BEFORE_DIR"],ENV["CHARGED_GAP_AFTER_DIR"],ENV["CHARGED_GAP_OUTPUT_DIR"])
    else
        main()
    end
end
end
