"""Hash-verified decision for the retained-background infinite-thermal GBU candidate.

This reducer performs no new equilibrium solve and grants no production promotion.
"""
module CausalGBUInfiniteReadiness
using CSV,JSON3,SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
hashfile(path)=bytes2hex(sha256(read(path)))

function verify_bundle(directory)
    m=JSON3.read(read(joinpath(directory,"manifest.json"),String))
    for (path,hash) in pairs(m.source_hashes)
        hashfile(joinpath(directory,"source_snapshot",String(path)))==String(hash) || error("snapshot mismatch: $(path)")
    end
    for (path,hash) in pairs(m.output_hashes)
        hashfile(joinpath(directory,String(path)))==String(hash) || error("output mismatch: $(path)")
    end
    if hasproperty(m,:background)
        for (path,hash) in pairs(m.background.input_hashes)
            hashfile(joinpath(String(m.background.input_directory),String(path)))==String(hash) || error("background mismatch")
        end
    end
    return m,(directory=abspath(directory),manifest_sha256=hashfile(joinpath(directory,"manifest.json")),
        source_files=length(m.source_hashes),output_files=length(m.output_hashes))
end

"""Weak-limit acceptance: fixed IR first, decreasing eta, then IR sensitivity."""
function eta_gate(rows;relative_target=0.01)
    isempty(rows) && return false
    cases=unique((r.channel,r.q_inv_fm) for r in rows)
    for (ch,q) in cases
        selected=filter(r->r.channel==ch && r.q_inv_fm==q,rows)
        irs=sort(unique(r.lower_inv_fm for r in selected))
        etas=sort(unique(r.eta_inv_fm for r in selected);rev=true)
        length(irs)>=2 && length(etas)>=3 && first(irs)>0 && last(etas)>0 || return false
        length(selected)==length(irs)*length(etas) || return false
        for ir in irs
            s=sort(filter(r->r.lower_inv_fm==ir,selected);by=r->r.eta_inv_fm,rev=true)
            length(unique(r.eta_inv_fm for r in s))==length(etas) || return false
            all(r->isfinite(r.relative_difference) && r.relative_difference>=0,s) || return false
            all(s[j+1].relative_difference<=s[j].relative_difference for j in 1:length(s)-1) || return false
            last(s).relative_difference<relative_target || return false
        end
        finest=filter(r->r.eta_inv_fm==last(etas),selected)
        pv=first(finest).pv_density
        isfinite(pv) && pv!=0 && all(r->r.pv_density==pv && isfinite(r.density),finest) || return false
        (maximum(r.density for r in finest)-minimum(r.density for r in finest))/abs(pv)<relative_target || return false
    end
    return true
end

function main()
    base=joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    names=["direction_b_infinite_oracle_20260907","direction_b_infinite_acceptance_v2_20260907",
        "direction_b_infinite_mott_refined_20260907","direction_b_infinite_gate_recovery_20260907"]
    bundles=[verify_bundle(joinpath(base,n)) for n in names]
    manifests=first.(bundles)
    readrows(j,file)=collect(CSV.File(joinpath(base,names[j],file)))
    probes=readrows(1,"probes.csv"); localrows=readrows(1,"shells.csv")
    acceptance=readrows(2,"acceptance.csv"); etas=readrows(2,"eta_limits.csv")
    decision=readrows(4,"decision.csv"); shells=readrows(4,"shells.csv")
    # Exact matrix completeness matters: all([]) must not certify an interrupted run.
    checks=[
        (gate="infinite_cut_radial_profile",passed=length(probes)==20 && length(localrows)==4 &&
            manifests[1].failures==0 && all(r.passed for r in probes) && all(r.relative_change<0.01 for r in localrows)),
        (gate="representative_near_axis_and_roots",passed=length(acceptance)==4 && manifests[2].passed &&
            all(r.numerics_passed && r.contour_passed && r.root_disks_passed for r in acceptance)),
        (gate="PV_eta_weighted_limit",passed=length(etas)==24 && eta_gate(etas)),
        (gate="normal_temperature_Mott",passed=manifests[3].rows==8 && manifests[3].all_passed),
        (gate="four_channel_weighted_q_integral",passed=length(decision)==4 &&
            Set(r.channel for r in decision)==Set(["pi_plus","pi_minus","K_plus","K_minus"]) &&
            length(shells)==416 && count(r->r.topology_evaluated,shells)==288 &&
            all(r.topology_passed for r in shells) && manifests[4].failures==0 && all(r.passed for r in decision))]
    output=get(ENV,"GBU_INFINITY_READINESS_OUTPUT",joinpath(base,"direction_b_infinite_readiness_20260907"))
    ispath(output) && error("refusing to overwrite $(output)")
    mkpath(output)
    CSV.write(joinpath(output,"checks.csv"),checks)
    source=relpath(@__FILE__,ROOT); snapshot=joinpath(output,"source_snapshot",source)
    mkpath(dirname(snapshot));cp(@__FILE__,snapshot)
    passed=all(r.passed for r in checks)
    record=Dict("status"=>passed ? "fixed_background_research_method_feasible" : "not_accepted",
        "all_passed"=>passed,"evidence"=>last.(bundles),"source_hashes"=>Dict(source=>hashfile(@__FILE__)),
        "output_hashes"=>Dict("checks.csv"=>hashfile(joinpath(output,"checks.csv"))),
        "thermal_target"=>"infinity","background_scope"=>"T170_MeV_muB240_MeV_quark_only_BQS_plus_retained_Mott_brackets",
        "L10_comparison_scope"=>"normal_branch_window_comparison_not_full_hard_endpoint_root_error_bound",
        "root_scope"=>"normal_analytic_gaps_and_sampled_cut_topology_not_interval_arithmetic_global_proof",
        "q_tail_scope"=>"successive_band_numerical_convergence_not_analytic_bound",
        "production_authorized"=>false,"default_provider_changed"=>false,"full_freezeout_accepted"=>false,
        "full_feedback_claimed"=>false,"solver_called"=>false)
    open(joinpath(output,"manifest.json"),"w") do io;JSON3.write(io,record);end
    println("[infinite-readiness] passed=$(passed) checks=$(count(r->r.passed,checks))/$(length(checks)) output=$(output)")
    passed || error("infinite thermal research readiness not accepted")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
