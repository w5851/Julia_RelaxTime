"""Automatically recover failed threshold-limit rows and reduce the full gate.

Permitted old-source changes are the analytic phase-limit guard and adaptive
gap contour windows/resolution. The physical kernel and yield formula are fixed.
Previously successful rows retain their original provenance; no CSV is edited.
"""
module CausalGBUInfiniteGateRecovery
include("causal_gbu_infinite_qgate.jl")
const Q=CausalGBUInfiniteQGate
const Y=Q.Y
const P=Q.P
const I=Q.I
const R=Q.R
using CSV,JSON3

function evaluate(bg,ch,q,weight,label,n,mesh)
    k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=max(36.,q+28bg.T+2.))
    p=P.profile(k;mesh=mesh);s=Y.shell(p;nodes=96)
    f=q>=20 ? 0. : Y.shell(P.profile(Y.finite_kernel(bg,ch,q,10.;cut_nodes=96);mesh=mesh);nodes=96).density
    topo=label=="bulk" ? Q.topology(p) : (passed=true,cut_crossings=-1,minimum_crossing_inverse=NaN)
    return (channel=String(ch),band=label,order=n,q_inv_fm=q,quadrature_weight=weight,
        density=s.density,bound=s.bound,landau=s.landau,pair=s.pair,finite_L10_density=f,
        root_count=s.root_count,topology_evaluated=label=="bulk",topology_passed=topo.passed,
        cut_crossings=topo.cut_crossings,crossing_inverse=topo.minimum_crossing_inverse,
        omega_tail_bound=s.omega_tail_conditional_bound,production_authorized=false,
        provenance="recomputed_with_analytic_limit_and_adaptive_contour")
end

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    previous=get(ENV,"GBU_INFINITY_GATE_PREVIOUS",joinpath(base,"direction_b_infinite_production_gate_20260907"))
    output=get(ENV,"GBU_INFINITY_RECOVERY_OUTPUT",joinpath(base,"direction_b_infinite_gate_recovery_20260907"))
    manifest=JSON3.read(read(joinpath(previous,"manifest.json"),String))
    allowed=normpath.(["scripts/analysis/relaxtime/causal_gbu_infinite_yield.jl",
        "scripts/analysis/relaxtime/causal_gbu_infinite_qgate.jl"])
    all(normpath(String(path)) in allowed || R.hashfile(joinpath(R.ROOT,String(path)))==String(hash)
        for (path,hash) in pairs(manifest.source_hashes)) || error("physics/source drift beyond the limit guard")
    for (path,hash) in pairs(manifest.output_hashes)
        R.hashfile(joinpath(previous,String(path)))==String(hash) || error("previous output drift")
    end
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    all(String(manifest.background.input_hashes[Symbol(path)])==hash for (path,hash) in bg.input_hashes) || error("background drift")
    old=NamedTuple.(collect(CSV.File(joinpath(previous,"shells.csv"))))
    mott=collect(CSV.File(joinpath(previous,"mott_momenta.csv")))
    key(r)=(String(r.channel),String(r.band),r.order,r.q_inv_fm)
    available=Dict(key(r)=>r for r in old)
    channels=unique(String(r.channel) for r in mott)
    orders=Int.(manifest.q_orders);mesh=Int(manifest.mesh)
    hashes=R.start_output(output)
    rows,totals,failures=NamedTuple[],NamedTuple[],NamedTuple[]; reused=0;recomputed=0
    for channel in channels
        ch=Symbol(channel)
        edges=sort!(unique!(vcat([0.,2bg.vacuum,8.],[r.q_inv_fm for r in mott if r.channel==channel])))
        bands=vcat([(0.,8.,n,"bulk") for n in orders],
            [(8.,12.,8,"tail1"),(12.,16.,8,"tail2"),(16.,24.,8,"tail3"),(24.,32.,8,"tail4")])
        for (lo,hi,n,label) in bands
            pieces=label=="bulk" ? collect(zip(edges[1:end-1],edges[2:end])) : [(lo,hi)]
            count=0;start=length(rows)
            for (a,b) in pieces
                qs,ws=R.gauleg(a,b,n)
                for (q,w) in zip(qs,ws)
                    count+=1;k=(channel,label,n,q)
                    try
                        if haskey(available,k) && available[k].topology_passed
                            push!(rows,merge(available[k],(provenance="retained_hash_bound_pre_limit_guard",)));reused+=1
                        else
                            push!(rows,evaluate(bg,ch,q,w,label,n,mesh));recomputed+=1
                        end
                    catch err
                        err isa InterruptException && rethrow()
                        push!(failures,(channel=channel,band=label,order=n,q_inv_fm=q,reason=sprint(showerror,err)))
                    end
                end
            end
            band=rows[start+1:end];complete=length(band)==count
            sumfield(f)=complete ? sum(r.quadrature_weight*getproperty(r,f) for r in band) : NaN
            push!(totals,(channel=channel,band=label,order=n,complete=complete,
                density=sumfield(:density),bound=sumfield(:bound),landau=sumfield(:landau),pair=sumfield(:pair),
                finite_L10_density=sumfield(:finite_L10_density),omega_tail_bound=sumfield(:omega_tail_bound),
                topology_passed=complete && all(r.topology_passed for r in band),production_authorized=false))
            isempty(rows) || CSV.write(joinpath(output,"shells.csv"),rows)
            CSV.write(joinpath(output,"integrals.csv"),totals)
            isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
            println("[gate-recovery] $(channel) $(label) n=$(n) reused=$(reused) recomputed=$(recomputed) failed=$(length(failures))");flush(stdout)
        end
    end
    decisions=NamedTuple[]
    for channel in channels
        t=filter(r->r.channel==channel,totals)
        lo=only(filter(r->r.band=="bulk" && r.order==minimum(orders),t))
        hi=only(filter(r->r.band=="bulk" && r.order==maximum(orders),t))
        tails=[only(filter(r->r.band=="tail$j",t)) for j in 1:4]
        n=hi.density+sum(r.density for r in tails)
        chosen=filter(r->r.channel==channel && (r.band!="bulk" || r.order==maximum(orders)),rows)
        weighted_L_error=sum(r.quadrature_weight*abs(r.finite_L10_density-r.density) for r in chosen)
        relative=abs(hi.density-lo.density)/abs(n)
        component=max(abs(hi.bound-lo.bound),abs(hi.pair-lo.pair),abs(hi.landau-lo.landau))/abs(n)
        last_fraction=abs(last(tails).density)/abs(n)
        decreasing=all(abs(tails[j+1].density)<abs(tails[j].density) for j in 1:3)
        omega_bound=hi.omega_tail_bound+sum(r.omega_tail_bound for r in tails)
        complete=all(r.complete && r.topology_passed for r in t)
        passed=complete && n>0 && relative<0.005 && component<0.01 && weighted_L_error/abs(n)<0.005 &&
            decreasing && last_fraction<1e-6 && omega_bound/abs(n)<1e-6
        push!(decisions,(channel=channel,density_inv_fm3=n,q_relative_change=relative,
            largest_component_change_over_total=component,L10_weighted_absolute_error=weighted_L_error,
            L10_relative_error=weighted_L_error/abs(n),q_above8_fraction=sum(abs(r.density) for r in tails)/abs(n),
            final_q_band_fraction=last_fraction,omega_tail_bound=omega_bound,tail_bands_decrease=decreasing,
            passed=passed,production_authorized=false))
    end
    CSV.write(joinpath(output,"decision.csv"),decisions)
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    R.finish_output(output,hashes,Dict("status"=>isempty(failures) && all(r.passed for r in decisions) ? "fixed_background_numerically_feasible" : "production_gate_not_passed",
        "background"=>bg,"previous_directory"=>abspath(previous),"previous_manifest_sha256"=>R.hashfile(joinpath(previous,"manifest.json")),
        "rows"=>length(rows),"reused_rows"=>reused,"recomputed_rows"=>recomputed,"failures"=>length(failures),
        "all_passed"=>all(r.passed for r in decisions),"allowed_source_changes"=>allowed,
        "reuse_reason"=>"analytic threshold guard and adaptive contour window/resolution only; physical kernel and density formula unchanged",
        "solver_called"=>false,"q_tail_assessment"=>"successive_band_numerical_convergence_not_interval_bound",
        "thermal_target"=>"infinity","full_feedback_claimed"=>false,"full_freezeout_authorized"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
