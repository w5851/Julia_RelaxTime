"""Analysis-only fast Cauchy representation, with threshold and infinite-tail terms.

Interpolation is a numerical approximation and requires comparison to the raw
infinite-domain oracle. Physical vacuum jumps are retained, not tapered away.
"""
module CausalGBUInfiniteProfile
include("causal_gbu_infinite_thermal.jl")
const I=CausalGBUInfiniteThermal
const R=I.R
const CB=Main.RelaxTime.CausalSpectralBubble

"""Integral of sqrt(u)*(1-u/H)/(u-s), from zero to H (no 1/pi)."""
function threshold_transform(s,H)
    H>0 && isfinite(H) && isfinite(s) || throw(ArgumentError("finite s and H>0 required"))
    if abs(s)>4H
        term=one(complex(s)); value=zero(complex(s))
        for n in 0:28
            value+=term/((n+1.5)*(n+2.5))
            term*=H/s
        end
        return -H^1.5/s*value
    end
    s==0 && return complex(4sqrt(H)/3)
    s==H && return complex(-2sqrt(H)/3)
    if imag(s)==0 && 0<real(s)<H
        x=real(s)
        i0=2sqrt(H)+sqrt(x)*log(abs((sqrt(H)-sqrt(x))/(sqrt(H)+sqrt(x))))
        return complex((1-x/H)*i0-2sqrt(H)/3,pi*sqrt(x)*(1-x/H))
    end
    root=sqrt(-complex(s))
    i0=2sqrt(H)-2root*atan(sqrt(H)/root)
    return (1-s/H)*i0-2sqrt(H)/3
end

function threshold_terms(k)
    i,j=R.charged_rpa_spec(k.ch).pair
    a,b=k.bg.m[i],k.bg.m[j]; Q=k.q; U=k.threshold
    E1,E2=a*U/(a+b),b*U/(a+b)
    margins=(k.bg.vacuum-a*Q/(a+b),k.bg.vacuum-b*Q/(a+b))
    minimum(abs,margins)>1e-10 || throw(ArgumentError("threshold meets vacuum projection edge"))
    vacuum=all(>(0),margins) ? 1. : 0.
    n1=I.P.occupations(E1,k.bg.mu[i],k.bg.T,k.bg.Phi,k.bg.PhiBar)
    n2=I.P.occupations(E2,k.bg.mu[j],k.bg.T,k.bg.Phi,k.bg.PhiBar)
    thermal=hasproperty(k,:thermal_limit) &&
        max(a*Q/(a+b),b*Q/(a+b))>=k.thermal_limit ? 0. : 1.
    pref=3sqrt(2)*(a*b)^1.5*sqrt(U)/(pi*(a+b)^2)
    # A compact analytic window, with its artificial endpoint continuous at zero.
    H=min(1.,k.split-U)/2
    return (U=U,H=H,positive=pref*(vacuum-thermal*(n1.quark+n2.anti)),
        negative=-pref*(vacuum-thermal*(n1.anti+n2.quark)))
end
function threshold_cut(t,x)
    u=abs(x)-t.U
    0<u<t.H || return zero(x)
    return (x>0 ? t.positive : t.negative)*sqrt(u)*(1-u/t.H)
end
function threshold_value(t,z)
    plus=threshold_transform(z-t.U,t.H)
    minus=imag(z)==0 ? conj(threshold_transform(-t.U-real(z),t.H)) :
        conj(threshold_transform(conj(-t.U-z),t.H))
    return (t.positive*plus-t.negative*minus)/pi
end

struct Profile{K,T}
    kernel::K
    threshold::T
    left::Vector{Float64}
    right::Vector{Float64}
    yleft::Vector{Float64}
    yright::Vector{Float64}
    tailx::Vector{Float64}
    tailpositive::Vector{Float64}
    tailnegative::Vector{Float64}
end

"""Cosine nodes with exact endpoints; never append a backwards rounded cell."""
function panel_nodes(a,b,n)
    a<b && n>=1 || throw(ArgumentError("ordered panel and positive order required"))
    xs=[a+(b-a)*sinpi(i/(2n))^2 for i in 0:n]
    xs[1]=a;xs[end]=b
    xs=sort!(unique!(xs))
    first(xs)==a && last(xs)==b && all(diff(xs).>0) || error("unordered panel nodes")
    return xs
end

function profile(k;mesh=256,tail_nodes=96)
    mesh>=16 && tail_nodes>=16 || throw(ArgumentError("mesh/tail_nodes >=16 required"))
    th=threshold_terms(k)
    edges=sort!(unique!(vcat(k.edges,[-th.U-th.H,th.U+th.H])))
    aa,bb,yy,zz=Float64[],Float64[],Float64[],Float64[]
    residual(x)=k.rho(x)-threshold_cut(th,x)
    for j in 1:length(edges)-1
        a,b=edges[j],edges[j+1]
        # Numerically tiny panels are retained as a single linear cell.
        n=b-a<1e-10 ? 1 : mesh
        xs=panel_nodes(a,b,n)
        ys=Float64[residual(x) for x in xs]
        if k.q==0 && b-a>1e-8
            # Interior extrapolation avoids mistaking a rounded physical edge
            # for a one-sided limit on the opposite side of the actual jump.
            h=min((b-a)/8,max((b-a)*1e-7,128eps(max(abs(a),abs(b),1.))))
            i1,i2=R.charged_rpa_spec(k.ch).pair
            edge=hypot(k.bg.m[i1],k.bg.vacuum)+hypot(k.bg.m[i2],k.bg.vacuum)
            abs(abs(a)-edge)<8eps(edge) && (ys[1]=2residual(a+h)-residual(a+2h))
            abs(abs(b)-edge)<8eps(edge) && (ys[end]=2residual(b-h)-residual(b-2h))
        end
        # Every finite-q edge is continuous. Use one shared nodal value on
        # both sides, never extrapolate two inconsistent artificial limits.
        append!(aa,xs[1:end-1]);append!(bb,xs[2:end])
        append!(yy,ys[1:end-1]);append!(zz,ys[2:end])
    end
    ts,ws=R.gauleg(0.,1.,tail_nodes); scale=2k.bg.T
    xx=k.split .+ scale.*ts./(1 .-ts)
    ww=ws.*scale./(1 .-ts).^2 ./pi
    return Profile(k,th,aa,bb,yy,zz,xx,ww.*k.rho.(xx),ww.*k.rho.(-xx))
end

function residual_cut(p,x)
    j=searchsortedlast(p.left,x)
    1<=j<=length(p.left) && x<=p.right[j] || return 0.
    return p.yleft[j]+(p.yright[j]-p.yleft[j])*(x-p.left[j])/(p.right[j]-p.left[j])
end
cut(p,x)=residual_cut(p,x)+threshold_cut(p.threshold,x)

function polarization(p::Profile,z)
    # The interpolant itself is Float64. Wide panel integration coordinates
    # must not promote every cell/log operation to BigFloat.
    z=ComplexF64(z)
    isfinite(z) && imag(z)>=0 && abs(real(z))<p.kernel.split ||
        throw(ArgumentError("closed UHP inside split window required"))
    value=0.0im
    w=real(z); r=imag(z)==0 ? residual_cut(p,w) : 0.
    if imag(z)==0 && p.kernel.q==0
        i,j=R.charged_rpa_spec(p.kernel.ch).pair
        bg=p.kernel.bg
        edge=hypot(bg.m[i],bg.vacuum)+hypot(bg.m[j],bg.vacuum)
        abs(abs(w)-edge)>8eps(edge) || throw(ArgumentError("PV at vacuum jump is undefined"))
    end
    for j in eachindex(p.left)
        a,b=p.left[j],p.right[j]; ya,yb=p.yleft[j],p.yright[j]
        if imag(z)>0
            t=(b-a)/(a-z); lg=log1p(t)
            value+=ya*lg+(yb-ya)*CB._cauchy_cell_moment(t,lg)
        elseif w==a || w==b
            abs((w==a ? ya : yb)-r)<1e-9 ||
                throw(ArgumentError("PV jump at $(w): endpoint=$(w==a ? ya : yb), shared=$(r)"))
            value+=yb-ya
        else
            t=(b-a)/(a-w)
            lg=t > -1 ? log1p(t) : log(abs((b-w)/(a-w)))
            value+=(ya-r)*lg+(yb-ya)*CB._cauchy_cell_moment(t,lg)
        end
    end
    if imag(z)==0
        value+=r*log(abs((last(p.right)-w)/(first(p.left)-w)))+im*pi*r
    end
    tail=sum(p.tailpositive[j]/(p.tailx[j]-z)+p.tailnegative[j]/(-p.tailx[j]-z)
        for j in eachindex(p.tailx))
    return value/pi+threshold_value(p.threshold,z)+tail
end
inverse(p,z)=1-4p.kernel.bg.coupling[p.kernel.ch]*polarization(p,z)
end
