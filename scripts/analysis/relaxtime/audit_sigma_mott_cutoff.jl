#!/usr/bin/env julia
# quote-remote only. Attribute finite-grid changes to the moving hard s cutoff.
include(joinpath(@__DIR__, "run_sigma_mott_spectrum_v2.jl"))

function sigma_cutoff_audit()
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"s6.toml"))
    runtime=SP_OTHER.build_runtime(profile)
    old=collect(CSV.File(joinpath(INPUT_DIR,"previous","observables.csv")))
    rate_rows=collect(CSV.File(joinpath(INPUT_DIR,"previous","rates.csv")))
    for name in ("left_s6_v1","right_s6_v1")
        append!(old,collect(CSV.File(joinpath(dirname(OUTPUT_DIR),name,"observables.csv"))))
        append!(rate_rows,collect(CSV.File(joinpath(dirname(OUTPUT_DIR),name,"rates.csv"))))
    end
    plan=collect(CSV.File(joinpath(INPUT_DIR,"transport_plan.csv")))
    result,events=NamedTuple[],NamedTuple[]
    for (h,suffix) in ((2e-4,"200"),(1e-4,"100"),(5e-5,"050"),(2e-5,"020"))
        ids=["sig_m"*suffix,"sig_c000","sig_p"*suffix]
        x=[Float64(only(r for r in plan if String(r.point_id)==id).xi) for id in ids]
        rows=[only(r for r in old if abs(Float64(r.xi)-xi)<1e-12) for xi in x]
        states=[frozen_state(row,profile,32,16) for row in rows]
        for process in (:uu_to_uu,:ud_to_ud)
            problems=[AS._prepare_thermal_cache_problem(process,st.qp,st.tp,st.kc;runtime.tau_kwargs...) for st in states]
            plans=[p.plan for p in problems]
            groups=AS.angular_groups(plans[2].cos_relative,plans[2].phiw)
            bounds=[p.context.bounds for p in problems]
            correction=zeros(3);changed=zeros(Int,3);tested=0
            # Unfolded i,j enumeration is common to all backgrounds. It is
            # algebraically the same sum as the optional equal-mass folding.
            for i in eachindex(plans[2].p),j in eachindex(plans[2].p),group in groups
                angular=nothing
                cut_cosines=ntuple(3) do n
                    p=plans[n]
                    (p.mi^2+p.mj^2+2p.ei[i]*p.ej[j]-last(bounds[n]))/(2p.p[i]*p.p[j])
                end
                low,high=extrema(cut_cosines)
                # Only angular nodes between the three actual cutoff values
                # can change membership. No full thermal measure is built.
                possible=findall(c->low-1e-13<=c<=high+1e-13,group.cosines)
                isempty(possible) && continue
                for k in possible
                    invariant=ntuple(3) do n
                        p=plans[n]
                        p.mi^2+p.mj^2+2*(p.ei[i]*p.ej[j]-p.p[i]*p.p[j]*group.cosines[k])
                    end
                    accepted=ntuple(n->first(bounds[n])<invariant[n]<last(bounds[n]),3)
                    all(==(accepted[2]),accepted) && continue
                    # Only the upper physical cutoff is examined. Threshold
                    # rounding events and distribution underflow are excluded.
                    all(s>first(b)+1e-6 for (s,b) in zip(invariant,bounds)) || continue
                    all(abs(s-last(b))<0.02 for (s,b) in zip(invariant,bounds)) || continue
                    if angular===nothing
                        angular=[sum(p.cw[a]*p.cw[b]*p.fi[i,a]*p.fj[j,b] for (a,b) in group.pairs) for p in plans]
                    end
                    tested+=1
                    for side in (1,3)
                        accepted[side]==accepted[2] && continue
                        p,st,s=plans[side],states[side],invariant[side]
                        relative=AS._relative_velocity_for_mode(:cm,s,p.mi,p.mj,p.ei[i],p.ej[j])
                        weight=p.pw[i]*p.pw[j]*p.p[i]^2*p.p[j]^2*p.prefactor*angular[side]*group.weights[k]*relative
                        weight==0 && continue
                        SIG_PC.reset_cache!()
                        sigma=SIG_CS.total_cross_section(process,s,st.qp,st.tp,st.kc;n_points=24)
                        delta=weight*sigma*(Int(accepted[2])-Int(accepted[side]))
                        correction[side]+=delta;changed[side]+=1
                        push!(events,(h=h,channel=String(process),point_id=ids[side],i=i,j=j,
                            s_inv_fm2=s,s_upper_inv_fm2=last(bounds[side]),accepted_here=accepted[side],
                            accepted_center=accepted[2],thermal_weight=weight,direct_sigma_fm2=sigma,
                            fixed_membership_delta_rate_fm2=delta))
                    end
                end
            end
            rates=[Float64(only(r for r in rate_rows if String(r.point_id)==String(row.point_id) && Symbol(r.channel)==process).rate)
                   for row in rows]
            adjusted=rates.+correction
            jump(v)=(log1p((v[3]-v[2])/v[2])-log1p((v[2]-v[1])/v[1]))/h
            push!(result,(h=h,channel=String(process),raw_jump=jump(rates),fixed_membership_jump=jump(adjusted),
                support_membership_jump=jump(rates)-jump(adjusted),left_delta_rate_fm2=correction[1],
                right_delta_rate_fm2=correction[3],left_changed_events=changed[1],right_changed_events=changed[3],
                candidate_events=tested))
            CSV.write(joinpath(OUTPUT_DIR,"cutoff_effect.csv"),result)
            CSV.write(joinpath(OUTPUT_DIR,"cutoff_events.csv"),events)
            println("CUTOFF_DONE h=",h," ",process," raw=",jump(rates)," fixed=",jump(adjusted));flush(stdout)
        end
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",cases=length(result),
        interpretation="Fixed event membership at the center is a diagnostic counterfactual, not a changed cutoff or a production rate. It isolates staircase artifacts of a moving sharp integration boundary at fixed momentum/angular quadrature.",
        production_kernel_changed=false,calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        script_sha256=digest(@__FILE__),files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    sigma_cutoff_audit()
end
