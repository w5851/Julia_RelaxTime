"""Analysis-only thermal-regulator and continuous spectral checks.

rho means Im(Pi^R), not a meson number density. No clipping or regulator
replacement is performed. Quadrature errors below are estimates, not bounds.
"""
module CausalGBUThermalAdmissibility
include("causal_gbu_pnjl_stability.jl")
const H=CausalGBUPNJLStability
const R=H.R
const P=H.P

"""PNJL f=mean(N)/3; df/dx=-var(N)/3 for nonnegative weights."""
function occupation_statistics(x,phi,bar)
    isfinite(x) && 0<=phi<=1 && 0<=bar<=1 || throw(ArgumentError("invalid occupation inputs"))
    y=exp(-abs(x))
    weights=x>=0 ? (1.,3phi*y,3bar*y^2,y^3) : (y^3,3phi*y^2,3bar*y,1.)
    norm=sum(weights)
    mean=sum((i-1)*weights[i] for i in 1:4)/norm
    variance=sum(((i-1)-mean)^2*weights[i] for i in 1:4)/norm
    return (occupation=mean/3,derivative=-variance/3)
end

"""Positive-lambda q0 pair cut from its exact on-shell Jacobian.

For Lambda<p<Lthermal, W=-(n1+nbar2)<0 at finite T: the vacuum has
ended while the thermal subtraction has not. This is an open interval,
not a value at the hard-cutoff discontinuity.
"""
function q0_pair(p,a,u,b,v,T,phi,bar,vc,tc;Nc=3)
    P.validate(0.,a,u,b,v,T,phi,bar,vc,4)
    isfinite(p) && p>=0 && isfinite(tc) && tc>=vc && Nc>0 ||
        throw(ArgumentError("invalid q0 pair inputs"))
    e1,e2=hypot(p,a),hypot(p,b)
    lambda=e1+e2
    n1,n2=P.occupations(e1,u,T,phi,bar),P.occupations(e2,v,T,phi,bar)
    vacuum=p<vc ? 1. : 0.
    thermal=p<tc ? -(n1.quark+n2.anti) : 0.
    jacobian=2pi*p/lambda
    projection=Nc/(8pi^2)*(lambda^2-(a-b)^2)
    return (lambda_inv_fm=lambda,k0_inv_fm=lambda-u+v,
        vacuum_weight=vacuum,thermal_weight=thermal,
        pair_weight=vacuum+thermal,imaginary_inv_fm2=projection*jacobian*(vacuum+thermal),
        thermal_only_interval=vc<p<tc,production_authorized=false)
end

function direct_cut(bg,ch,q,lambda,thermal;Phi=bg.Phi,PhiBar=bg.PhiBar,nodes=96)
    a,b=R.charged_rpa_spec(ch).pair
    vac=P.reference_cut(lambda,q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T,bg.vacuum;
        Phi=Phi,PhiBar=PhiBar,component=:vacuum,nodes=nodes)
    th=P.reference_cut(lambda,q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T,thermal;
        Phi=Phi,PhiBar=PhiBar,component=:thermal,nodes=nodes)
    factor=3/(8pi^2)*(lambda^2-q^2-(bg.m[a]-bg.m[b])^2)
    return (vacuum_pair=factor*vac.pair,thermal_pair=factor*th.pair,
        landau=factor*(vac.landau+th.landau),imaginary=factor*(vac.imaginary+th.imaginary))
end

"""Adaptive two-order Gauss check on supplied physical panels.

The summed embedded difference is an ERROR ESTIMATE, never a uniform bound.
"""
function integrate_panels(f,edges;nodes=12,atol=1e-9,max_depth=12)
    length(edges)>=2 && all(isfinite,edges) && all(diff(edges).>0) &&
        nodes>=4 && isfinite(atol) && atol>0 && max_depth>=0 ||
        throw(ArgumentError("invalid quadrature settings"))
    function panel(a,b,depth)
        function rule(n)
            xs,ws=R.gauleg(a,b,n)
            return sum(w*f(x) for (x,w) in zip(xs,ws))
        end
        coarse,fine=rule(nodes),rule(2nodes)
        all(isfinite,(coarse,fine)) || throw(ArgumentError("nonfinite spectral integrand"))
        return (left=a,right=b,value=fine,error_estimate=abs(fine-coarse),depth=depth)
    end
    # Control the requested TOTAL absolute error. Proportional local budgets
    # spuriously fail on tiny kinematic panels despite a small aggregate error.
    parts=[panel(edges[i],edges[i+1],0) for i in 1:length(edges)-1]
    evaluations=3nodes*length(parts)
    error=sum(p.error_estimate for p in parts)
    while error>atol
        eligible=[i for i in eachindex(parts) if parts[i].depth<max_depth &&
            parts[i].left<(parts[i].left+parts[i].right)/2<parts[i].right]
        isempty(eligible) && break
        j=eligible[argmax([parts[i].error_estimate for i in eligible])]
        p=parts[j]
        p.error_estimate==0 && break
        mid=(p.left+p.right)/2
        parts[j]=panel(p.left,mid,p.depth+1)
        push!(parts,panel(mid,p.right,p.depth+1))
        evaluations+=6nodes
        error=sum(part.error_estimate for part in parts)
    end
    return (value=sum(p.value for p in parts),error_estimate=error,
        converged=error<=atol,evaluations=evaluations,panel_count=length(parts),
        depth_limited_panels=count(p->p.depth==max_depth,parts),rigorous_error_bound=false)
end

"""PV/UHP integral of a direct cut, without spectral interpolation.

Subtract rho(Re z) analytically to resolve the near-axis Cauchy kernel.
Caller supplies all support/cutoff breakpoints; probes at discontinuities
are outside the declared contract.
"""
function continuous_cauchy(rho,edges,z;nodes=12,atol=1e-9,max_depth=12)
    length(edges)>=2 && all(isfinite,edges) && all(diff(edges).>0) ||
        throw(ArgumentError("ordered finite spectral endpoints required"))
    isfinite(z) && imag(z)>=0 || throw(ArgumentError("closed UHP required"))
    a,b=first(edges),last(edges)
    c=real(z)
    interior=a<c<b
    imag(z)==0 && !interior && (c==a || c==b) && throw(ArgumentError("endpoint probe"))
    r=interior ? rho(c) : 0.
    panels=Float64[edges...]
    interior && push!(panels,c)
    if imag(z)>0
        h=imag(z)
        while h<b-a
            for x in (c-h,c+h)
                a<x<b && push!(panels,x)
            end
            h*=2
        end
    end
    sort!(unique!(panels))
    result=integrate_panels(x->(rho(x)-r)/(x-z),panels;nodes=nodes,atol=pi*atol,max_depth=max_depth)
    logarithm=imag(z)>0 ? log(complex(b)-z)-log(complex(a)-z) :
        complex(log(abs((b-c)/(a-c))),interior ? pi : 0.)
    return (value=(result.value+r*logarithm)/pi,
        error_estimate=result.error_estimate/pi,converged=result.converged,
        evaluations=result.evaluations,rigorous_error_bound=false)
end

"""Partition hints from the existing kinematic support, NOT spectral ordinates.

The direct integrand still comes from the independent Pereira cut. This
partition helper is not an independent proof of all support breakpoints.
"""
function cut_panels(g)
    q,a,b=g.q_inv_fm,g.m1_inv_fm,g.m2_inv_fm
    upper=hypot(a,g.thermal_cutoff_inv_fm)+hypot(b,g.thermal_cutoff_inv_fm)
    edges=Float64[-upper,0.,upper,g.mu1_inv_fm-g.mu2_inv_fm]
    for e in (hypot(q,a+b),hypot(q,a-b))
        e<upper && append!(edges,[-e,e])
    end
    for cutoff in (g.vacuum_cutoff_inv_fm,g.thermal_cutoff_inv_fm)
        q<2cutoff || continue
        s=Main.RelaxTime.OneLoopIntegrals._spectral_two_line_support(q,a,b,cutoff)
        for e in (s.pair...,s.landau...)
            -upper<e<upper && append!(edges,[-e,e])
        end
        for (p,r) in ((cutoff,cutoff),(cutoff,abs(cutoff-q)),(abs(cutoff-q),cutoff)),
                si in (-1,1),ti in (-1,1)
            e=si*hypot(a,p)-ti*hypot(b,r)
            -upper<e<upper && push!(edges,e)
        end
    end
    filter!(x->-upper<=x<=upper,edges)
    return sort!(unique!(edges))
end

"""Exact algebraic RPA sign transfer; contact terms cannot fix a cut sign."""
function rpa_spectral(pi_value,K)
    isfinite(pi_value) && isfinite(K) && K>0 || throw(ArgumentError("invalid RPA inputs"))
    inverse=1-4K*pi_value
    abs2(inverse)>0 || throw(ArgumentError("pole needs separate distributional treatment"))
    d=2K/inverse
    predicted=8K^2*imag(pi_value)/abs2(inverse)
    return (imaginary=imag(d),identity_value=predicted,inverse_abs=abs(inverse),
        production_authorized=false)
end

"""Rouche transfer requires a genuine UNIFORM error bound, not sampled drift."""
function count_transfer_margin(minimum_inverse,error;uniform_bound_certified=false)
    all(isfinite,(minimum_inverse,error)) && min(minimum_inverse,error)>=0 ||
        throw(ArgumentError("invalid count transfer inputs"))
    return (margin=minimum_inverse-error,
        count_transfer_certified=uniform_bound_certified && error<minimum_inverse,
        numerical_margin_positive=error<minimum_inverse)
end
end
