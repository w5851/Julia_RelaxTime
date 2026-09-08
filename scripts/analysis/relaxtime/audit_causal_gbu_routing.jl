"""Relative-momentum regulator probe for bound-state dispersion, not an alternative density provider."""
module CausalGBURouting
include("causal_gbu_research_utils.jl")
using .CausalGBUResearch
using CSV
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
const R=CausalGBUResearch
const NC=Main.Constants_PNJL.N_color

function centered_grid(q,m1,mu1,m2,mu2,T,phi,phibar,vc,tc;np=128,nx=64)
    q>=0 && 0<vc<=tc && T>0 && min(m1,m2)>0 || throw(ArgumentError("invalid centered regulator inputs"))
    poles,weights=Float64[],Float64[]
    contact=0.0
    for (vacuum,L) in ((true,vc),(false,tc))
        ps,pws=R.gauleg(0.0,L,np)
        xs,xws=q==0 ? ([0.0],[2.0]) : R.gauleg(-1.0,1.0,nx)
        for (p,pw) in zip(ps,pws), (x,xw) in zip(xs,xws)
            e1=sqrt(p^2+q^2/4+p*q*x+m1^2)
            e2=sqrt(p^2+q^2/4-p*q*x+m2^2)
            n1=vacuum ? (1.0,0.0) : (-antiquark_distribution(e1,mu1,T,phi,phibar),quark_distribution(e1,mu1,T,phi,phibar))
            n2=vacuum ? (1.0,0.0) : (-antiquark_distribution(e2,mu2,T,phi,phibar),quark_distribution(e2,mu2,T,phi,phibar))
            measure=pw*xw*p^2/(e1*e2)
            contact+=measure*2*(e1*(n2[2]-n2[1])+e2*(n1[2]-n1[1]))
            for (si,s) in enumerate((-1,1)), (ti,t) in enumerate((-1,1))
                u=s*e1-t*e2
                residue=NC/(8pi^2)*measure*s*t*(n2[ti]-n1[si])*(u^2-q^2-(m1-m2)^2)
                residue==0 && continue
                push!(poles,u); push!(weights,residue)
            end
        end
    end
    return (poles=poles,weights=weights,contact=contact)
end
value(g,z)=sum((g.weights[i]/(z-g.poles[i]) for i in eachindex(g.poles));init=0.0+0im)

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_RESEARCH_OUTPUT",joinpath(base,"method_v1_routing"))
    hashes=R.start_output(output)
    rows=NamedTuple[]
    for channel in R.CHANNELS, q in (0.0,0.01,0.02,0.05,0.1,0.25,0.5,1.0,2.0,3.0), np in (96,192)
        a,b=R.charged_rpa_spec(channel).pair
        shift=bg.mu[a]-bg.mu[b]
        old=R.bubble_at(bg,channel,q,R.Settings(mesh=256,np=np))
        g=centered_grid(q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T,bg.Phi,bg.PhiBar,bg.vacuum,20.0;np=np)
        centered(z)=1-4bg.coupling[channel]*value(g,z+shift)
        left,right=max(0.0,old.landau),old.threshold
        roots=R.certify_gap_roots((w,_)->centered(w),q,[(left,right)];physical_sheet=true,real_axis=true)
        original=R.gap_audit(old,R.Settings();count_contour=false)
        h=min(0.02,(right-left)/20)
        count=R.contour_count(centered,left+1e-6,right-1e-6,-h,h)
        push!(rows,(channel=String(channel),q_inv_fm=q,momentum_nodes=np,
            two_line_count=original.count,centered_count=roots.count,centered_contour_count=count.count,
            two_line_root=isempty(original.roots) ? NaN : first(original.roots).omega_inv_fm,
            centered_root=isempty(roots.roots) ? NaN : first(roots.roots).omega_inv_fm,
            centered_passed=roots.passed && count.passed && roots.count==count.count,
            two_line_contact=old.grid.contact_inv_fm2,centered_contact=g.contact,
            interpretation="regulator_sensitivity_of_gap_poles_only",density_computed=false,production_authorized=false))
        CSV.write(joinpath(output,"routing_gap_comparison.csv"),rows)
        println("[gbu-routing] $(channel) q=$(q) p=$(np) roots=$(original.count)/$(roots.count)")
    end
    R.finish_output(output,hashes,Dict("status"=>"routing_probe_not_density_method","background"=>bg,
        "regulators"=>["both |p| and |p-q| < Lambda","relative |p-q/2| < Lambda"],
        "thermal_max_inv_fm"=>20.0,"solver_called"=>false,"density_computed"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
