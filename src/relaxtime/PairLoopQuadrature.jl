module PairLoopQuadrature

using ..GaussLegendre: gauleg
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: correction_cos_theta_coefficient

export PairLoopPlan, pair_b0_pv, pair_occupation

"""Reusable fixed-background quadrature for the equal-mass k=0 loop.

Only the numerical representation changes: the historical positive-energy
cut, finite three-momentum cutoff, and linear-xi loop occupations are retained.
No complex-energy extension, finite-k kernel or unequal-mass prescription.
"""
struct PairLoopPlan
    m::Float64
    mu::Float64
    T::Float64
    Phi::Float64
    Phibar::Float64
    xi::Float64
    cutoff::Float64
    p2::Vector{Float64}
    weights::Vector{Float64}
    occupation::Vector{Float64}
    F0::Float64
    threshold_integral::Float64
end

function pair_occupation(y,m,mu,T,Phi,Phibar,xi)
    E=sqrt(m*m+y)
    p=sqrt(y)
    fq=quark_distribution(E,mu,T,Phi,Phibar)
    fa=antiquark_distribution(E,mu,T,Phi,Phibar)
    if xi!=0
        fq+=correction_cos_theta_coefficient(:quark,p,m,mu,T,Phi,Phibar,xi)/3
        fa+=correction_cos_theta_coefficient(:antiquark,p,m,mu,T,Phi,Phibar,xi)/3
    end
    return (1-fq-fa)/E
end
pair_occupation(plan::PairLoopPlan,y) =
    pair_occupation(y,plan.m,plan.mu,plan.T,plan.Phi,plan.Phibar,plan.xi)

function PairLoopPlan(m::Float64,mu::Float64,tp;
    n_nodes::Int=64,cutoff::Float64=Main.Constants_PNJL.Λ_inv_fm)
    all(isfinite,(m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ,cutoff)) &&
        m>0 && tp.T>0 && cutoff>0 && n_nodes>=8 ||
        throw(ArgumentError("finite inputs, positive mass/T/cutoff and n_nodes>=8 required"))
    p,w=gauleg(0.0,cutoff,n_nodes)
    p2=p.^2
    f=[pair_occupation(y,m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ) for y in p2]
    f0=pair_occupation(0.0,m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ)
    return PairLoopPlan(m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ,cutoff,p2,w,f,f0,sum(w.*f))
end

function occupation_derivative(plan,y)
    h=cbrt(eps(Float64))*max(plan.m^2,y,1.0)
    lo=max(0.0,y-h); hi=y+h
    return (pair_occupation(plan,hi)-pair_occupation(plan,lo))/(hi-lo)
end

function pair_b0_pv(plan::PairLoopPlan,s::Float64)
    isfinite(s) && s>=0 || throw(ArgumentError("finite nonnegative invariant s required"))
    q2=(s-4plan.m^2)/4
    L=plan.cutoff
    q2==L^2 && throw(DomainError(s,"upper sharp-cutoff pole has a logarithmic divergence"))
    if q2==0
        return complex(2plan.threshold_integral,0.0)
    elseif q2<0
        kappa=sqrt(-q2)
        # Endpoint subtraction captures the nonanalytic -pi*kappa*F(0) term.
        regular=sum(w*(f-plan.F0)/(y+kappa^2)
            for (y,w,f) in zip(plan.p2,plan.weights,plan.occupation))
        return complex(2plan.threshold_integral-2kappa^2*regular-
                       2kappa*plan.F0*atan(L/kappa),0.0)
    else
        q=sqrt(q2)
        fq=pair_occupation(plan,q2)
        regular=0.0; derivative=NaN
        for (y,w,f) in zip(plan.p2,plan.weights,plan.occupation)
            d=y-q2
            divided = if abs(d)<=sqrt(eps(Float64))*max(plan.m^2,y,q2,1.0)
                isnan(derivative) && (derivative=occupation_derivative(plan,q2))
                derivative
            else
                (f-fq)/d
            end
            regular+=w*divided
        end
        logarithm=q<L ? log1p(-2q/(L+q)) : log1p(-2L/(q+L))
        re=2plan.threshold_integral+2q2*regular+q*fq*logarithm
        im=q<L ? π*q*fq : 0.0
        return complex(re,im)
    end
end

end # module PairLoopQuadrature
