"""Direction-B prerequisites: static derivatives and a subtraction toy.

Neither the toy nor the static identities define a production response,
a new equilibrium, or a GBU density.
"""
module CausalGBUDirectionB
include("causal_gbu_thermal_admissibility.jl")
const A=CausalGBUThermalAdmissibility
const R=A.R
const P=A.P
using ForwardDiff

"""PNJL color-traced log determinant; Nc=3 is in the polynomial."""
function log_partition(x,phi,bar)
    y=exp(-abs(x))
    return x>=0 ? log1p(3phi*y+3bar*y^2+y^3) :
        -3x+log1p(3bar*y+3phi*y^2+y^3)
end

"""One-flavor kinetic Omega; T, Phi and static cutoffs are held fixed."""
function static_flavor(m,mu,T,phi,bar,vc,tc;nodes=64)
    all(isfinite,(m,mu,T,phi,bar,vc,tc)) && min(m,T,vc)>0 && tc>=vc &&
        0<=phi<=1 && 0<=bar<=1 && nodes>=8 || throw(ArgumentError("invalid static inputs"))
    ps,ws=R.gauleg(0.,vc,nodes)
    vacuum=-3/pi^2*sum(w*p^2*hypot(p,m) for (p,w) in zip(ps,ws))
    vcond=-3/pi^2*sum(w*p^2*m/hypot(p,m) for (p,w) in zip(ps,ws))
    ps,ws=R.gauleg(0.,tc,nodes)
    thermal=zero(m+mu)
    condensate=vcond
    density=zero(m+mu)
    for (p,w) in zip(ps,ws)
        E=hypot(p,m)
        x,y=(E-mu)/T,(E+mu)/T
        nq,na=P.occupation(x,phi,bar),P.occupation(y,bar,phi)
        thermal-=T/pi^2*w*p^2*(log_partition(x,phi,bar)+log_partition(y,bar,phi))
        condensate+=3/pi^2*w*p^2*m/E*(nq+na)
        density+=3/pi^2*w*p^2*(nq-na)
    end
    return (omega_inv_fm4=vacuum+thermal,condensate_inv_fm3=condensate,
        density_inv_fm3=density,production_authorized=false)
end

"""Twice-subtracted SYNTHETIC rho(t)=c*t tail, t>S (not project Pi).

Delta(s)=s^2/pi * integral_S^infty rho(t)/[t^2(t-s)] dt.
Keeps the constant and linear terms at s=0 but may introduce RPA poles.
"""
function matched_linear_tail(s,c,S)
    isfinite(s) && isfinite(c) && isfinite(S) && c>0 && S>0 ||
        throw(ArgumentError("invalid tail parameters"))
    if isreal(s)
        s==S && throw(ArgumentError("hard tail endpoint"))
        return real(s)>S ? -c*s/pi*complex(log(s/S-1),-pi) :
            -c*s/pi*log1p(-s/S)
    end
    return -c*s/pi*log1p(-s/S)
end

"""Locate the SYNTHETIC RPA pole at z=iQ; no actual background is used."""
function toy_uhp_pole(c,S,K)
    all(isfinite,(c,S,K)) && min(c,S,K)>0 || throw(ArgumentError("invalid toy RPA"))
    f(Q)=1-4K*matched_linear_tail(-Q^2,c,S)
    left,right=0.,sqrt(S)
    for _ in 1:100
        f(right)<0 && break
        right*=2
    end
    f(right)<0 || error("toy root bracket failed")
    bracket=(left,right)
    for _ in 1:100
        mid=(left+right)/2
        f(mid)>0 ? (left=mid) : (right=mid)
    end
    q=(left+right)/2
    return (imaginary_frequency_inv_fm=q,residual=abs(f(q)),
        bracket_inv_fm=bracket,is_project_pole=false,production_authorized=false)
end
end
