"""Automatic q integration and numerical Lth error budget; diagnostic outputs only."""
module CausalGBUInfiniteDensityAudit
include("causal_gbu_infinite_yield.jl")
const Y=CausalGBUInfiniteYield
const P=Y.P
const I=Y.I
const R=Y.R
using CSV,JSON3

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    tag=get(ENV,"GBU_INFINITY_BACKGROUND","T170")
    bg=tag=="T170" ? R.frozen_background(joinpath(base,"negative_density_phase_fig2_like")) :
        R.saved_background(joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3"),parse(Float64,tag),48)
    output=get(ENV,"GBU_INFINITY_DENSITY_OUTPUT",joinpath(base,"direction_b_infinite_density_20260907"))
    channels=Symbol.(split(get(ENV,"GBU_INFINITY_CHANNELS","pi_plus,pi_minus,K_plus,K_minus"),','))
    mesh=parse(Int,get(ENV,"GBU_INFINITY_DENSITY_MESH","512"))
    orders=parse.(Int,split(get(ENV,"GBU_INFINITY_Q_ORDERS","8,16"),','))
    hashes=R.start_output(output)
    rows,totals,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    for ch in channels
        for (a,b,n,label) in vcat([(0.,8.,n,"bulk") for n in orders],[(8.,12.,8,"tail1"),(12.,16.,8,"tail2")])
            qs,ws=R.gauleg(a,b,n)
            start=length(rows)
            for (q,v) in zip(qs,ws)
                try
                    k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=36.)
                    p=P.profile(k;mesh=mesh)
                    s=Y.shell(p;nodes=96)
                    f=P.profile(Y.finite_kernel(bg,ch,q,10.;cut_nodes=96);mesh=mesh)
                    fs=Y.shell(f;nodes=96)
                    push!(rows,(channel=String(ch),band=label,order=n,q_inv_fm=q,quadrature_weight=v,
                        density=s.density,bound=s.bound,landau=s.landau,pair=s.pair,
                        finite_L10_density=fs.density,root_count=s.root_count,
                        relative_L10_difference=abs(s.density-fs.density)/max(abs(s.density),1e-12),
                        omega_tail_bound=s.omega_tail_conditional_bound,production_authorized=false))
                catch err
                    err isa InterruptException && rethrow()
                    push!(failures,(channel=String(ch),band=label,order=n,q_inv_fm=q,reason=sprint(showerror,err)))
                end
                isempty(rows) || CSV.write(joinpath(output,"shells.csv"),rows)
                isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
                println("[infinite-density] $(ch) $(label) n=$(n) q=$(q) rows=$(length(rows)) failed=$(length(failures))");flush(stdout)
            end
            band=rows[start+1:end]
            complete=length(band)==n
            sumfield(field)=complete ? sum(r.quadrature_weight*getproperty(r,field) for r in band) : NaN
            push!(totals,(channel=String(ch),band=label,order=n,complete=complete,
                density=sumfield(:density),bound=sumfield(:bound),landau=sumfield(:landau),pair=sumfield(:pair),
                finite_L10_density=sumfield(:finite_L10_density),omega_tail_bound=sumfield(:omega_tail_bound),
                production_authorized=false))
            CSV.write(joinpath(output,"integrals.csv"),totals)
        end
    end
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    R.finish_output(output,hashes,Dict("status"=>isempty(failures) ? "weighted_integrals_computed" : "density_evaluation_failed",
        "background"=>bg,"rows"=>length(rows),"failures"=>length(failures),"q_orders"=>orders,
        "mesh"=>mesh,"thermal_target"=>"infinity","finite_Lth_numerical_comparison"=>10.,
        "solver_called"=>false,"q_beyond_16_analytic_bound_certified"=>false,
        "relative_density_target"=>0.01,"full_production_authorized"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
