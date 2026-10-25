# Fixed-background quadrature of the three positive RTA transport kernels.
if !isdefined(Main,:ThresholdUpgradeCommon)
    include(joinpath(@__DIR__,"ThresholdUpgradeCommon.jl"))
end

module TerminalTransportQuadrature

using LinearAlgebra
using Main.ThresholdUpgradeCommon: AS
const TC=Main.RelaxTime.TransportCoefficients
const SPECIES=(:u,:d,:s,:ubar,:dbar,:sbar)

export terminal_grams, gram_relative_change, transport_from_grams, SPECIES

"""Gram matrix of (p^2,1,Ekin), with positive weight p^2/Ekin^2*f*(1-f).

The polar integral is divided by two. Hence G[1,1] is the shear moment,
G[1,2] the conductivity moment, and every bulk bracket is a quadratic form
in this same matrix. No background derivatives are needed to bound the
quadrature change for *all* bulk brackets at these fixed masses.
"""
function terminal_grams(qp,tp,p,pw,ct,cw)
    return NamedTuple{SPECIES}(map(SPECIES) do species
        m=AS.get_mass(species,qp);mu=AS.get_mu(species,qp)
        G=zeros(3,3)
        for (momentum,w) in zip(p,pw)
            p2=momentum^2;E=sqrt(p2+m*m)
            angular=0.0
            for (c,wc) in zip(ct,cw)
                f=AS.distribution_with_anisotropy(species,momentum,m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ,c)
                angular+=wc*f*(1-f)/2
            end
            base=w*p2/(E*E)*angular
            basis=(p2,1.0,E)
            for i in 1:3,j in 1:3
                G[i,j]+=base*basis[i]*basis[j]
            end
        end
        G
    end)
end

"""Maximum relative quadratic-form change, over every nonzero coefficient vector.

This is a numerical generalized-eigenvalue bound relative to the supplied
reference matrix, not a claim that the reference equals the exact integral.
Diagonal scaling limits avoidable conditioning from unequal units of the
three basis functions.
"""
function gram_relative_change(candidate,reference)
    scale=sqrt.(diag(reference))
    all(s->isfinite(s) && s>0,scale) || throw(DomainError(scale,"positive Gram diagonal required"))
    normalizer=scale*scale'
    R=Symmetric(reference./normalizer)
    D=Symmetric((candidate-reference)./normalizer)
    return maximum(abs,eigvals(D,R))
end

function bulk_vector(qp,tp,coeffs,species)
    flavor=species in (:u,:ubar) ? :u : species in (:d,:dbar) ? :d : :s
    i=flavor===:u ? 1 : flavor===:d ? 2 : 3
    m=getproperty(qp.m,flavor);mu=getproperty(qp.μ,flavor)
    m==coeffs.masses[i] || throw(ArgumentError("bulk coefficients and Gram masses must match"))
    v=coeffs.v_n_sq;muprime=coeffs.dμB_dT_sigma
    a=1-3v
    b=3v*(tp.T*m*(coeffs.dM_dT[i]+coeffs.dM_dμB[i]*muprime)-m*m)
    c=3v*mu-v*tp.T*muprime
    species in (:ubar,:dbar,:sbar) && (c=-c)
    return [a,b,c]
end

function transport_from_grams(grams,qp,tp,tau;bulk_coeffs=nothing)
    charges=TC.default_charges();g=TC.degeneracy_default()
    eta=g/(30pi^2*tp.T)*sum(getproperty(tau,s)*getproperty(grams,s)[1,1] for s in SPECIES)
    sigma=g/(6pi^2*tp.T)*sum(TC._q2_for_species(s,charges)*getproperty(tau,s)*getproperty(grams,s)[1,2] for s in SPECIES)
    zeta=bulk_coeffs===nothing ? NaN : Main.Constants_PNJL.N_color/(9pi^2*tp.T)*sum(SPECIES) do s
        b=bulk_vector(qp,tp,bulk_coeffs,s)
        getproperty(tau,s)*dot(b,getproperty(grams,s)*b)
    end
    return (eta=eta,sigma=sigma,zeta=zeta)
end

end
