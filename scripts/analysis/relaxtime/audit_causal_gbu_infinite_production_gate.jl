"""Mott-split, weighted numerical acceptance before research-production wiring."""
module CausalGBUInfiniteProductionGate
include("causal_gbu_infinite_qgate.jl")
const Q=CausalGBUInfiniteQGate
const Y=Q.Y
const P=Q.P
const I=Q.I
const R=Q.R
using CSV,JSON3

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    tag=get(ENV,"GBU_INFINITY_BACKGROUND","T170")
    bg=tag=="T170" ? R.frozen_background(joinpath(base,"negative_density_phase_fig2_like")) :
        R.saved_background(joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3"),parse(Float64,tag),48)
    output=get(ENV,"GBU_INFINITY_GATE_OUTPUT",joinpath(base,"direction_b_infinite_production_gate_20260907"))
    channels=Symbol.(split(get(ENV,"GBU_INFINITY_CHANNELS","pi_plus,pi_minus,K_plus,K_minus"),','))
    mesh=parse(Int,get(ENV,"GBU_INFINITY_DENSITY_MESH","512"))
    orders=parse.(Int,split(get(ENV,"GBU_INFINITY_Q_ORDERS","8,16"),','))
    check_topology=get(ENV,"GBU_INFINITY_TOPOLOGY","true")=="true"
    hashes=R.start_output(output)
    rows,totals,mott,failures=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    for ch in channels
        m=Q.mott_momenta(bg,ch)
        m.passed || error("Mott quadrature not converged")
        for r in m.checks
            push!(mott,merge((channel=String(ch),),r))
        end
        isempty(mott) || CSV.write(joinpath(output,"mott_momenta.csv"),mott)
        edges=sort!(unique!(vcat([0.,2bg.vacuum,8.],m.roots)))
        bands=vcat([(0.,8.,n,"bulk") for n in orders],
            [(8.,12.,8,"tail1"),(12.,16.,8,"tail2"),(16.,24.,8,"tail3"),(24.,32.,8,"tail4")])
        for (lo,hi,n,label) in bands
            pieces=label=="bulk" ? collect(zip(edges[1:end-1],edges[2:end])) : [(lo,hi)]
            count=0; start=length(rows)
            for (a,b) in pieces
                qs,ws=R.gauleg(a,b,n)
                for (q,v) in zip(qs,ws)
                    count+=1
                    try
                        k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=max(36.,q+28bg.T+2.))
                        p=P.profile(k;mesh=mesh)
                        s=Y.shell(p;nodes=96)
                        finite=q>=20 ? 0. : Y.shell(P.profile(Y.finite_kernel(bg,ch,q,10.;cut_nodes=96);mesh=mesh);nodes=96).density
                        topo=check_topology && label=="bulk" ? Q.topology(p) :
                            (passed=true,cut_crossings=-1,minimum_crossing_inverse=NaN,positive_count=-1,negative_count=-1)
                        push!(rows,(channel=String(ch),band=label,order=n,q_inv_fm=q,quadrature_weight=v,
                            density=s.density,bound=s.bound,landau=s.landau,pair=s.pair,finite_L10_density=finite,
                            root_count=s.root_count,topology_evaluated=check_topology && label=="bulk",topology_passed=topo.passed,
                            cut_crossings=topo.cut_crossings,crossing_inverse=topo.minimum_crossing_inverse,
                            omega_tail_bound=s.omega_tail_conditional_bound,production_authorized=false))
                    catch err
                        err isa InterruptException && rethrow()
                        push!(failures,(channel=String(ch),band=label,order=n,q_inv_fm=q,reason=sprint(showerror,err)))
                    end
                    isempty(rows) || CSV.write(joinpath(output,"shells.csv"),rows)
                    isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
                    println("[infinite-gate] $(ch) $(label) n=$(n) q=$(q) rows=$(length(rows)) failed=$(length(failures))");flush(stdout)
                end
            end
            band=rows[start+1:end];complete=length(band)==count
            sumfield(f)=complete ? sum(r.quadrature_weight*getproperty(r,f) for r in band) : NaN
            push!(totals,(channel=String(ch),band=label,order=n,complete=complete,
                density=sumfield(:density),bound=sumfield(:bound),landau=sumfield(:landau),pair=sumfield(:pair),
                finite_L10_density=sumfield(:finite_L10_density),omega_tail_bound=sumfield(:omega_tail_bound),
                topology_passed=complete && all(r.topology_passed for r in band),production_authorized=false))
            CSV.write(joinpath(output,"integrals.csv"),totals)
        end
    end
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    R.finish_output(output,hashes,Dict("status"=>isempty(failures) ? "weighted_gate_evaluated" : "gate_evaluation_failed",
        "background"=>bg,"rows"=>length(rows),"failures"=>length(failures),"q_orders"=>orders,"mesh"=>mesh,
        "thermal_target"=>"infinity","finite_Lth_numerical_comparison"=>10.,"solver_called"=>false,
        "outer_topology_evaluated"=>check_topology,"q_beyond_32_analytic_bound_certified"=>false,
        "relative_density_target"=>0.01,"full_production_authorized"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
