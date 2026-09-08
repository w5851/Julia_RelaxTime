"""Targeted normal-Mott refinement: one new background, no density scan."""
module CausalGBUMottRefinement
include("audit_causal_gbu_continuation_gate.jl")
include("audit_causal_gbu_validation.jl")
const C=CausalGBUContinuationGateAudit
const V=CausalGBUValidation
const R=C.R
using CSV

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    previous=joinpath(base,"direction_b_continuation_v2_20260907")
    input=joinpath(base,"method_v1_mott");bgpath=joinpath(input,"backgrounds.csv")
    previous_path=joinpath(previous,"normal_mott.csv")
    ih=R.hashfile(bgpath);ph=R.hashfile(previous_path)
    output=get(ENV,"GBU_MOTT_REFINEMENT_OUTPUT",joinpath(base,"direction_b_mott_refinement_20260907"))
    R.method_contract();hashes=R.start_output(output)
    old=collect(CSV.File(previous_path));rows=NamedTuple[]
    model=Main.Models.create_model(:PNJL)
    # Bounded refinement of the sole retained bracket wider than 0.25 MeV.
    kp=only(filter(r->r.channel=="K_plus" && r.q_inv_fm==1.,old))
    Tmid=(kp.T_low_MeV+kp.T_high_MeV)/2
    new=V.background(model,Tmid,240.,24)
    CSV.write(joinpath(output,"backgrounds.csv"),[new.row])
    k=C.kernel(new.bg,:K_plus,1.,10.)
    midpoint=C.direct_value(k,complex(k.threshold))
    fmid=real(1-4new.bg.coupling[:K_plus]*midpoint.value)
    CSV.write(joinpath(output,"new_threshold.csv"),[(channel="K_plus",q_inv_fm=1.,T_MeV=Tmid,
        threshold_inverse=fmid,pi_node_change_inv_fm2=midpoint.node_change,
        numerics_passed=midpoint.passed,production_authorized=false)])
    for row in old
        ch,q=Symbol(row.channel),row.q_inv_fm
        lo,hi=row.T_low_MeV,row.T_high_MeV
        passed=row.normal_mott_numerics_passed
        root_change=row.before_root_node_change
        if ch==:K_plus && q==1.
            if fmid<0
                lo=Tmid
            else
                hi=Tmid
            end
            passed=midpoint.passed && hi-lo<=0.25 &&
                row.before_count-row.after_count==1 && abs(row.phase_drop-1)<0.005 &&
                row.before_phase_node_passed && row.after_phase_node_passed &&
                row.before_gap_checks_passed && row.after_gap_checks_passed && root_change<1e-7
        elseif ch==:K_minus && q==1.
            bg=R.saved_background(input,row.T_pre_MeV,24)
            pre=C.normal_gap(bg,ch,q,10.;loop_nodes=(256,512))
            root_change=pre.root_change
            passed=pre.passed && pre.count-row.after_count==1 &&
                row.after_phase_node_passed && row.after_gap_checks_passed && hi-lo<=0.25 && abs(row.phase_drop-1)<0.005
        end
        margins=map((lo,hi,row.T_pre_MeV,row.T_post_MeV)) do T
            bg=T==Tmid ? new.bg : R.saved_background(input,T,24)
            i,j=R.charged_rpa_spec(ch).pair
            min(2bg.m[i],2bg.m[j],hypot(q,bg.m[i]+bg.m[j])-abs(bg.mu[i]-bg.mu[j]))
        end
        push!(rows,(channel=String(ch),q_inv_fm=q,T_low_MeV=lo,T_high_MeV=hi,
            before_count=row.before_count,after_count=row.after_count,phase_drop=row.phase_drop,
            root_node_change_inv_fm=root_change,reference_gap_lower_bound_inv_fm=minimum(margins),
            normal_mott_numerics_passed=passed,reference_band_gap_open=minimum(margins)>0,
            full_levinson_accepted=false,production_authorized=false))
        CSV.write(joinpath(output,"normal_mott_refined.csv"),rows)
        println("[mott-refine] $(ch) q=$(q) normal=$(passed)");flush(stdout)
    end
    R.hashfile(bgpath)==ih && R.hashfile(previous_path)==ph || error("retained input drift")
    R.finish_output(output,hashes,Dict("status"=>all(r.normal_mott_numerics_passed for r in rows) ? "normal_mott_refined" : "normal_mott_refinement_failed",
        "retained_input_directory"=>input,"retained_backgrounds_sha256"=>ih,
        "previous_directory"=>previous,"previous_normal_mott_sha256"=>ph,
        "rows"=>length(rows),"normal_mott_numerics_passed"=>all(r.normal_mott_numerics_passed for r in rows),
        "new_backgrounds"=>1,"solver_called"=>true,"meson_density_computed"=>false,
        "full_levinson_accepted"=>false,"full_UHP_certified"=>false,
        "limitations"=>["Normal low-energy roots only, excluding unresolved auxiliary exterior counts",
            "One p24 BQS solve refines a retained Mott bracket; no freezeout production",
            "Existing successful normal-gap checks reused with hashes; Kminus q1 loop orders256/512",
            "Adiabatic reference proof coverage is separate from full stationarity"] ))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
