"""Transfer retained normal-Mott brackets to the actual infinite thermal target."""
module CausalGBUInfiniteMottAudit
include("causal_gbu_infinite_yield.jl")
const Y=CausalGBUInfiniteYield
const P=Y.P
const I=Y.I
const R=Y.R
using CSV,JSON3

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    refined=joinpath(base,"direction_b_mott_refinement_20260907")
    previous=joinpath(base,"direction_b_continuation_v2_20260907","normal_mott.csv")
    input=joinpath(base,"method_v1_mott")
    output=get(ENV,"GBU_INFINITY_MOTT_OUTPUT",joinpath(base,"direction_b_infinite_mott_20260907"))
    bracketfile=joinpath(refined,"normal_mott_refined.csv")
    input_hashes=Dict(p=>R.hashfile(p) for p in (previous,bracketfile,joinpath(input,"backgrounds.csv"),joinpath(refined,"backgrounds.csv")))
    before=collect(CSV.File(previous));brackets=collect(CSV.File(bracketfile))
    hashes=R.start_output(output);rows=NamedTuple[];backgrounds=Any[]
    for r in brackets
        ch=Symbol(r.channel);q=r.q_inv_fm
        old=only(filter(x->x.channel==r.channel && x.q_inv_fm==q,before))
        probes=NamedTuple[]
        for (label,T) in (("low",r.T_low_MeV),("high",r.T_high_MeV),("pre",old.T_pre_MeV),("post",old.T_post_MeV))
            bg=R.saved_background(T==196.015625 ? refined : input,T,24)
            push!(backgrounds,bg)
            k=I.kernel(bg,ch,q;cut_nodes=192,split_inv_fm=36.)
            coarse=I.polarization(k,k.threshold;nodes=64)
            fine=I.polarization(k,k.threshold;nodes=128)
            f=real(1-4bg.coupling[ch]*fine)
            count=-1;phase=NaN
            if label in ("pre","post")
                p=P.profile(k;mesh=1024)
                count=length(Y.gap_roots(p)[2])
                phase=Y.cut_phase(p,k.threshold-k.shift+1e-11)
            end
            push!(probes,(label=label,T=T,f=f,change=abs(fine-coarse),count=count,phase=phase))
        end
        low,high,pre,post=probes
        passed=low.f<0<high.f && high.T-low.T<=0.25 && pre.count==1 && post.count==0 &&
            abs((pre.phase-post.phase)/pi-1)<0.005 && maximum(x.change for x in probes)<1e-6
        push!(rows,(channel=String(ch),q_inv_fm=q,T_low_MeV=low.T,T_high_MeV=high.T,
            low_threshold_inverse=low.f,high_threshold_inverse=high.f,before_count=pre.count,after_count=post.count,
            phase_drop_over_pi=(pre.phase-post.phase)/pi,node_change_inv_fm2=maximum(x.change for x in probes),
            passed=passed,production_authorized=false))
        CSV.write(joinpath(output,"normal_mott.csv"),rows)
        println("[infinite-mott] $(ch) q=$(q) passed=$(passed)");flush(stdout)
    end
    all(R.hashfile(p)==h for (p,h) in input_hashes) || error("Mott inputs drifted")
    R.finish_output(output,hashes,Dict("status"=>all(r.passed for r in rows) ? "normal_mott_infinite_passed" : "normal_mott_infinite_failed",
        "rows"=>length(rows),"input_hashes_absolute"=>input_hashes,"backgrounds"=>backgrounds,
        "all_passed"=>all(r.passed for r in rows),"solver_called"=>false,
        "thermal_target"=>"infinity","full_freezeout_authorized"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
