"""Direct-cut, discontinuity-aware quadrature for conditional Mott diagnostics."""
module CausalGBUContinuationGate
include("causal_gbu_endpoint_closure.jl")
const E=CausalGBUEndpointClosure
const A=E.A
const R=E.R

# A local generic affine rule: do not add methods to the shared gauleg function.
function wide_rule(f,a,b,n)
    ts,ws=R.gauleg(-1.,1.,n)
    return sum(BigFloat(w)*(b-a)/2*f((a+b)/2+(b-a)/2*BigFloat(t)) for (t,w) in zip(ts,ws))
end

function wide_intervals(A,B,C,lo,hi)
    lo<hi || return Tuple{BigFloat,BigFloat}[]
    es=BigFloat[lo,hi]
    if A==0
        B!=0 && lo< -C/B <hi && push!(es,-C/B)
    elseif B^2-4A*C>=0
        Q=-(B+copysign(sqrt(B^2-4A*C),B))/2
        for x in (Q==0 ? (-B/(2A),) : (Q/A,C/Q))
            lo<x<hi && push!(es,x)
        end
    end
    sort!(unique!(es))
    return [(es[i],es[i+1]) for i in 1:length(es)-1 if
        (x=(es[i]+es[i+1])/2; (A*x+B)*x+C<0)]
end

"""Same Pereira cut with wide coordinates for sub-ulp panel ordinates.

This preserves panels, including jumps; no clamping or deletion is used.
Gauss weights remain Float64, so this is NOT arbitrary-precision certification.
"""
function wide_cut(bg,ch,q,lambda,thermal;nodes=64)
    i,j=R.charged_rpa_spec(ch).pair
    a,u,b,v,T,phi,bar,Q,z=BigFloat.((bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T,bg.Phi,bg.PhiBar,q,lambda))
    total=BigFloat(0)
    for (component,L0) in ((:vacuum,bg.vacuum),(:thermal,thermal))
        L=BigFloat(L0)
        Q<2L || continue
        for eta in (-1,1)
            S=-eta*z
            if Q>0
                interval=S>0 ? E.P.epsilon_interval(S/2,Q,a,b,L) : nothing
                if interval!==nothing
                    total-=eta*BigFloat(pi)/(2Q)*wide_rule(interval...,nodes) do epsilon
                        n1,n2,c=E.P.weights((S-epsilon)/2,(S+epsilon)/2,u,v,T,phi,bar,component)
                        eta==1 ? c-n1.anti-n2.quark : c-n1.quark-n2.anti
                    end
                end
                lo=max(a+S/2,b-S/2);hi=min(hypot(a,L)+S/2,hypot(b,L)-S/2)
                d=a^2-b^2
                for (l,h) in wide_intervals(4*(S^2-Q^2),4S*d,Q^4+2Q^2*(a^2+b^2)-Q^2*S^2+d^2,lo,hi)
                    total+=eta*BigFloat(pi)/Q*wide_rule(l,h,nodes) do en
                        n1,n2,_=E.P.weights(en-S/2,en+S/2,u,v,T,phi,bar,component)
                        eta==1 ? n1.quark-n2.quark : n1.anti-n2.anti
                    end
                end
            elseif z!=0
                if S>a+b
                    e1=(S^2+a^2-b^2)/(2S);e2=S-e1
                    if a<e1<hypot(a,L) && b<e2<hypot(b,L)
                        p=sqrt((S^2-(a+b)^2)*(S^2-(a-b)^2))/(2S)
                        n1,n2,c=E.P.weights(e1,e2,u,v,T,phi,bar,component)
                        total-=eta*2BigFloat(pi)*p/abs(z)*(eta==1 ? c-n1.anti-n2.quark : c-n1.quark-n2.anti)
                    end
                end
                e1=(b^2-a^2-S^2)/(2S);e2=e1+S
                if a<e1<hypot(a,L) && b<e2<hypot(b,L)
                    n1,n2,_=E.P.weights(e1,e2,u,v,T,phi,bar,component)
                    total+=eta*2BigFloat(pi)*sqrt(e1^2-a^2)/abs(z)*(eta==1 ? n1.quark-n2.quark : n1.anti-n2.anti)
                end
            end
        end
    end
    return 3/(8BigFloat(pi)^2)*(z^2-Q^2-(a-b)^2)*total
end

"""Cosine-map each physical panel to resolve integrable square-root edges.

Panel discontinuities are retained; endpoints are never sampled. Comparison
between orders is a numerical error estimate, not a uniform analytic bound.
"""
function mapped_integral(f,edges;nodes=48)
    nodes>=8 && length(edges)>=2 && all(isfinite,edges) && all(diff(edges).>0) ||
        throw(ArgumentError("ordered panels and nodes>=8 required"))
    ts,ws=R.gauleg(0.,1.,nodes)
    total=0.0im
    for j in 1:length(edges)-1
        a,b=edges[j],edges[j+1]
        for (t,w) in zip(ts,ws)
            # sin^2 avoids subtractive loss at a square-root endpoint.
            x=a+(b-a)*sinpi(t/2)^2
            value=if !(a<x<b) || min(x-a,b-x)<8eps(max(abs(a),abs(b),1.))
                setprecision(128) do
                    xx=BigFloat(a)+(BigFloat(b)-BigFloat(a))*sinpi(BigFloat(t)/2)^2
                    a<xx<b || throw(ArgumentError("mapped panel has an unrepresentable wide ordinate"))
                    f(xx)
                end
            else
                f(x)
            end
            isfinite(value) || throw(ArgumentError("nonfinite mapped integrand"))
            total+=w*(b-a)*pi/2*sinpi(t)*ComplexF64(value)
        end
    end
    return total
end

"""PV/retarded boundary using the original cut, not interpolated ordinates."""
function mapped_cauchy(rho,edges,z;nodes=48)
    isfinite(z) && imag(z)>=0 || throw(ArgumentError("closed UHP required"))
    left,right=first(edges),last(edges)
    x=real(z)
    interior=left<x<right
    # At a global endpoint only a genuinely vanishing one-sided cut is allowed.
    r=interior ? rho(x) : 0.
    ps=copy(edges)
    interior && push!(ps,x)
    sort!(unique!(ps))
    integrand(t)=t isa BigFloat ? (rho(t)-(interior ? rho(BigFloat(x)) : 0))/(t-Complex{BigFloat}(z)) : (rho(t)-r)/(t-z)
    value=mapped_integral(integrand,ps;nodes=nodes)
    if r!=0
        lg=imag(z)>0 ? log(complex(right)-z)-log(complex(left)-z) :
            complex(log(abs((right-x)/(left-x))),pi)
        value+=r*lg
    end
    return value/pi
end

function source_domain(bg,ch)
    i,j=R.charged_rpa_spec(ch).pair
    margins=(u=bg.m.u-abs(bg.mu.u),d=bg.m.d-abs(bg.mu.d),s=bg.m.s-abs(bg.mu.s))
    return (pair_margin_inv_fm=min(margins[i],margins[j]),
        all_flavor_margin_inv_fm=minimum(margins),
        previous_pair_source_proof_applicable=min(margins[i],margins[j])>0,
        previous_full_source_proof_applicable=minimum(margins)>0)
end
end
