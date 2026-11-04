# Numerical candidates for the controlled threshold study. No default backend
# is changed by loading this module.
if !isdefined(Main, :ThresholdUpgradeCommon)
    include(joinpath(@__DIR__, "ThresholdUpgradeCommon.jl"))
end

module ThresholdCacheCandidates

using Main.ThresholdUpgradeCommon
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: correction_cos_theta_coefficient
using Main.RelaxTime.PairLoopQuadrature: PairLoopPlan, pair_b0_pv, pair_occupation, occupation_derivative
const RT = Main.RelaxTime

export AdaptiveSigmaTable, build_adaptive_sigma, PairLoopPlan, pair_b0_pv,
       pair_occupation, paired_polarizations, threshold_estimates,
       paired_cross_section, functional_rate, composite_sigma

struct AdaptiveSigmaTable{R}
    nodes::Vector{Float64}
    values::Vector{Float64}
    sigma::Vector{Float64}
    estimated_error::Float64
    weighted_integral::Float64
    evaluations::Int
    converged::Bool
end

function (table::AdaptiveSigmaTable{R})(x::Float64) where {R}
    nodes=table.nodes
    first(nodes) <= x <= last(nodes) || throw(DomainError(x,"outside adaptive table"))
    i=searchsortedfirst(nodes,x)
    i==1 && return table.sigma[1]
    nodes[i]==x && return table.sigma[i]
    a,b=nodes[i-1],nodes[i]
    t=(x-a)/(b-a)
    value=(1-t)*table.values[i-1]+t*table.values[i]
    return R===:x_sigma ? value/x : value
end

@inline function segment_sigma(representation,a,b,fa,fb,x)
    t=(x-a)/(b-a)
    representation===:x_sigma && return ((1-t)*a*fa+t*b*fb)/x
    return (1-t)*fa+t*fb
end

"""Refine using an absolute weighted interpolation-defect estimate.

Quarter/midpoint probes use composite Simpson weights; after bisection the
parent quarter points become child midpoint probes and are reused. The
estimate is not a rigorous bound for an arbitrary unresolved peak; supply
physical feature seeds and verify with independent quadrature/thermal rates.
Returns converged=false on node-budget exhaustion; consumers must check it.
"""
function build_adaptive_sigma(f, seeds::AbstractVector{<:Real};
    representation::Symbol=:x_sigma, weight=identity, rtol=1e-3, atol=1e-10,
    max_nodes::Int=1024, known_values=Dict{Float64,Float64}())
    representation in (:raw,:x_sigma) || throw(ArgumentError("representation must be raw or x_sigma"))
    isfinite(rtol) && isfinite(atol) && rtol>=0 && atol>=0 && rtol+atol>0 ||
        throw(ArgumentError("finite nonnegative tolerances with at least one positive are required"))
    nodes=sort!(unique(Float64.(seeds)))
    length(nodes)>=2 && all(isfinite,nodes) && first(nodes)>0 ||
        throw(ArgumentError("at least two distinct finite positive x nodes are required"))
    max_nodes>=length(nodes) || throw(ArgumentError("max_nodes is smaller than seed count"))
    evaluated=copy(known_values)
    evaluations=Ref(0)
    function sample(x)
        value=get!(evaluated,x) do
            evaluations[]+=1
            Float64(f(x))
        end
        isfinite(value) && value>=0 || throw(DomainError((x,value),"nonfinite or negative cross section"))
        return value
    end
    sigma=sample.(nodes)
    function measure(a,b,fa,fb)
        wa,wb=Float64(weight(a)),Float64(weight(b))
        all(isfinite,(wa,wb)) && wa>=0 && wb>=0 ||
            throw(DomainError((wa,wb),"weight must be finite and nonnegative"))
        err=0.0; integral=(b-a)*(wa*fa+wb*fb)/12
        midpoint=(a+b)/2
        for (x,w) in (((a+midpoint)/2,4.0),(midpoint,2.0),((midpoint+b)/2,4.0))
            actual=sample(x)
            wx=Float64(weight(x))
            isfinite(wx) && wx>=0 || throw(DomainError(wx,"weight must be finite and nonnegative"))
            qweight=(b-a)*w/12*wx
            approx=segment_sigma(representation,a,b,fa,fb,x)
            err+=qweight*abs(actual-approx)
            integral+=qweight*actual
        end
        return err,integral
    end
    measures=[measure(nodes[i],nodes[i+1],sigma[i],sigma[i+1]) for i in 1:length(nodes)-1]
    errors=first.(measures); integrals=last.(measures)
    while sum(errors)>atol+rtol*abs(sum(integrals)) && length(nodes)<max_nodes
        i=argmax(errors)
        a,b=nodes[i],nodes[i+1]
        mid=(a+b)/2
        a<mid<b || break
        value=sample(mid)
        left=measure(a,mid,sigma[i],value)
        right=measure(mid,b,value,sigma[i+1])
        insert!(nodes,i+1,mid); insert!(sigma,i+1,value)
        errors[i]=left[1]; insert!(errors,i+1,right[1])
        integrals[i]=left[2]; insert!(integrals,i+1,right[2])
    end
    error=sum(errors); integral=sum(integrals)
    values=representation===:x_sigma ? nodes.*sigma : copy(sigma)
    return AdaptiveSigmaTable{representation}(nodes,values,sigma,error,integral,evaluations[],
        error<=atol+rtol*abs(integral))
end

function paired_polarizations(state,s,plans)
    qp=state.quark_params
    factor=-Main.Constants_PNJL.N_color/(8π^2)
    function one(flavor)
        m=getproperty(qp.m,flavor); a=getproperty(qp.A,flavor)
        plan=getproperty(plans,flavor); tp=state.thermo_params
        (plan.m,plan.mu,plan.T,plan.Phi,plan.Phibar,plan.xi,plan.cutoff)==
            (m,getproperty(qp.μ,flavor),tp.T,tp.Φ,tp.Φbar,tp.ξ,Main.Constants_PNJL.Λ_inv_fm) ||
            throw(ArgumentError("PairLoopPlan does not match the current physical background"))
        b=pair_b0_pv(plan,s)
        return (P=factor*(2a-s*b),S=factor*(2a+(4m*m-s)*b))
    end
    return (u=one(:u),s=one(:s))
end

function threshold_estimates(state; n_nodes=128)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    plans=(u=PairLoopPlan(qp.m.u,qp.μ.u,tp;n_nodes=n_nodes),
           s=PairLoopPlan(qp.m.s,qp.μ.s,tp;n_nodes=n_nodes))
    s0=4qp.m.u^2
    function denominators(s)
        pol=paired_polarizations(state,s,plans)
        mm=RT.EffectiveCouplings.mixing_matrix_elements(pol.u.P,pol.s.P,kc,:P)
        return (pi=1-4kc.K123_plus*pol.u.P,
                mixed=(mm.M00*mm.M88-mm.M08^2)/kc.det_K_plus)
    end
    d0=denominators(s0); h=s0*1e-8; dh=denominators(s0+h)
    return [(family=String(key),delta=real(getproperty(d0,key)),
        c=-imag(getproperty(dh,key))/sqrt(h),
        x_delta=(real(getproperty(d0,key))*sqrt(h)/imag(getproperty(dh,key)))^2)
        for key in (:pi,:mixed)]
end

"""Diagnostic-only temporary s-channel polarization injection.

The four affected keys are restored even on failure. This is deliberately not
a production backend API; use only in an isolated single-thread process.
"""
function paired_cross_section(process,s,state,plans;n_points=24)
    process in (:uubar_to_uubar,:uubar_to_ddbar) ||
        throw(ArgumentError("paired diagnostic currently covers the two retained light channels"))
    pc=RT.PolarizationCache
    qp,tp=state.quark_params,state.thermo_params
    pol=paired_polarizations(state,s,plans)
    previous=Pair[]
    try
        for (flavor,ns) in ((:u,0),(:s,2)), channel in (:P,:S)
            m=getproperty(qp.m,flavor); mu=getproperty(qp.μ,flavor); a=getproperty(qp.A,flavor)
            key=pc.PolarizationKey(channel,sqrt(s),0.0,m,m,mu,mu,
                tp.T,tp.Φ,tp.Φbar,tp.ξ,a,a,ns)
            push!(previous,key=>get(pc.POLARIZATION_CACHE,key,nothing))
            value=getproperty(getproperty(pol,flavor),channel)
            pc.POLARIZATION_CACHE[key]=(real(value),imag(value))
        end
        return CS.total_cross_section(process,s,qp,tp,state.K_coeffs;n_points=n_points)
    finally
        for (key,value) in previous
            value===nothing ? delete!(pc.POLARIZATION_CACHE,key) : (pc.POLARIZATION_CACHE[key]=value)
        end
    end
end

"""Use a verified adaptive table locally, a supplied cache elsewhere, and
direct evaluation below its first positive x node. Every caller supplies s0."""
function composite_sigma(table::AdaptiveSigmaTable,s0,tail,direct)
    table.converged || throw(ArgumentError("adaptive table has not converged"))
    function sigma(s)
        x=s-s0
        x<=0 && return 0.0
        x<first(table.nodes) && return direct(s)
        x<=last(table.nodes) && return table(x)
        return tail(s)
    end
    return sigma
end

"""Same 5D thermal quadrature as AverageScatteringRate, with a callable sigma.

This analysis consumer enables a fair comparison of representations without
altering the production cache contract. Its parity is checked against the
production consumer before reporting candidate rates.
"""
function functional_rate(sigma,process,state;p_nodes=128,angle_nodes=20,phi_nodes=36)
    qp,tp=state.quark_params,state.thermo_params
    i,j,c,d=CS.parse_particles_from_process(process)
    mi,mj,mc,md=(AS.get_mass(q,qp) for q in (i,j,c,d))
    mui,muj=AS.get_mu(i,qp),AS.get_mu(j,qp)
    T,Phi,Phibar,xi=tp.T,tp.Φ,tp.Φbar,tp.ξ
    p,pw=AS._cached_semi_infinite_momentum_grid(p_nodes,AS.DEFAULT_SEMI_INF_SCALE)
    ct,cw=AS._cached_interval_grid(-1.0,1.0,angle_nodes)
    phi,phiw=AS._cached_interval_grid(0.0,2π,phi_nodes)
    _,sbo,sup=AS._resolve_sigma_support_bounds(mi,mj,mc,md,
        Main.Constants_PNJL.Λ_inv_fm,nothing,nothing)
    rhoi=AS.number_density(i,mi,mui,T,Phi,Phibar,xi;
        p_nodes=AS.DEFAULT_SEMI_INF_NODES,angle_nodes=angle_nodes)
    rhoj=AS.number_density(j,mj,muj,T,Phi,Phibar,xi;
        p_nodes=AS.DEFAULT_SEMI_INF_NODES,angle_nodes=angle_nodes)
    (rhoi==0 || rhoj==0) && return 0.0
    total=0.0
    for (pi,wpi) in zip(p,pw), (pj,wpj) in zip(p,pw)
        Ei=AS.energy_from_p(pi,mi); Ej=AS.energy_from_p(pj,mj)
        for (ci,wci) in zip(ct,cw), (cj,wcj) in zip(ct,cw)
            si=sqrt(max(1-ci*ci,0.0)); sj=sqrt(max(1-cj*cj,0.0))
            fi=AS.distribution_with_anisotropy(i,pi,mi,mui,T,Phi,Phibar,xi,ci)
            fj=AS.distribution_with_anisotropy(j,pj,mj,muj,T,Phi,Phibar,xi,cj)
            (fi==0 || fj==0) && continue
            for (az,waz) in zip(phi,phiw)
                s=mi^2+mj^2+2*(Ei*Ej-pi*pj*(ci*cj+si*sj*cos(az)))
                (s<=sbo || s>=sup) && continue
                root=sqrt(s)
                Eicm=(s+mi^2-mj^2)/(2root); Ejcm=(s-mi^2+mj^2)/(2root)
                pcm=sqrt(max(0.0,(s-(mi+mj)^2)*(s-(mi-mj)^2)))/(2root)
                vnum=(Eicm*Ejcm+pcm*pcm)^2-(mi*mj)^2
                v=vnum>0 ? sqrt(vnum)/(Eicm*Ejcm) : 0.0
                (v==0 || v>2) && continue
                base=wpi*wpj*wci*wcj*waz*pi^2*pj^2*fi*fj*v
                total+=base*sigma(s)
            end
        end
    end
    return AS.DQ^2/(32π^5*rhoi*rhoj)*total
end

end
