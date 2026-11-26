"""Independent fifth-chapter / MFIR diagnostic formulas, in fm natural units.

This module does not call Models, PNJLCore or MagneticIntegrals.  The thermal
integral is split around each Landau level's Fermi surface; the vacuum MFIR
condensate uses loggamma, independently of the production Abel-Plana kernel.
"""
module MagneticReference

using FastGaussQuadrature: gausslegendre
using SpecialFunctions: loggamma
using LinearAlgebra: norm, svdvals

const HBARC = 197.3269804
const CHARGES = (2 / 3, 1 / 3, 1 / 3)
const NODE_CACHE = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

Base.@kwdef struct Parameters
    cutoff::Float64 = 602.3 / HBARC
    G0::Float64 = 1.835 / (602.3 / HBARC)^2
    K::Float64 = 12.36 / (602.3 / HBARC)^5
    ml::Float64 = 5.5 / HBARC
    ms::Float64 = 140.7 / HBARC
    T0::Float64 = 210 / HBARC
    a0::Float64 = 3.51
    a1::Float64 = -2.47
    a2::Float64 = 15.2
    b3::Float64 = -1.75
    imc_a::Float64 = 0.0108805
    imc_b::Float64 = -1.0133e-4
    imc_c::Float64 = 0.02228
    imc_d::Float64 = 1.84558e-4
    qcd_scale::Float64 = 300 / HBARC
    colors::Int = 3
end

Base.@kwdef struct Controls
    order::Int = 32
    pz_max::Float64 = 40.0
    n_max::Int = 32
    cutoff_N::Int = 10
    thermal_width::Float64 = 6.0
    thermal_tail::Float64 = 40.0
    route::Symbol = :mfir
end

function coupling(eB, p::Parameters)
    z = eB / p.qcd_scale^2
    return p.G0 * (1 + p.imc_a * z^2 + p.imc_b * z^3) /
        (1 + p.imc_c * z^2 + p.imc_d * z^4)
end

function masses(phi, eB, p::Parameters)
    G = coupling(eB, p)
    u, d, s = phi[1:3]
    return [p.ml - 4G*u + 2p.K*d*s, p.ml - 4G*d + 2p.K*u*s,
        p.ms - 4G*s + 2p.K*u*d]
end

function mass_jacobian(phi, eB, p::Parameters)
    G = coupling(eB, p)
    u, d, s = phi[1:3]
    return [-4G 2p.K*s 2p.K*d; 2p.K*s -4G 2p.K*u; 2p.K*d 2p.K*u -4G]
end

# Euler-Maclaurin derivative of Hurwitz zeta at s=-1. The derivative of
# (s)_(2r-1) at -1 is -(2r-3)! for r>=2, so no s difference is used.
function zeta_derivative(x)
    x > 0 || throw(ArgumentError("positive zeta argument required"))
    shift = 24
    q = x + shift
    val = -sum((x+k)*log(x+k) for k in 0:shift-1) +
        q^2 * (log(q)/2 - 1/4) - q*log(q)/2 + (1+log(q))/12
    bernoulli = (-1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
    for (i, B) in enumerate(bernoulli)
        r = i+1
        val -= B*q^(2-2r)/((2r)*(2r-1)*(2r-2))
    end
    return val
end

function integrate(f, edges, order)
    nodes, weights = get!(NODE_CACHE, order) do
        gausslegendre(order)
    end
    acc = nothing
    for k in 2:length(edges)
        lo, hi = edges[k-1], edges[k]
        hi > lo || continue
        mid, halfwidth = (hi+lo)/2, (hi-lo)/2
        for j in eachindex(nodes)
            v = f(mid+halfwidth*nodes[j]) .* (halfwidth*weights[j])
            acc = acc === nothing ? v : acc .+ v
        end
    end
    return acc
end

# logZ, occupation f, dlogZ/dPhi, dlogZ/dPhibar, independently stabilized
# by rewriting the polynomial in exp(-abs((mu-E)/T)).
@inline function occupation(delta, Phi, Phibar)
    r = exp(-abs(delta))
    if delta <= 0
        z = 1+3Phi*r+3Phibar*r^2+r^3
        z > 0 || throw(DomainError(z, "Polyakov polynomial is nonpositive"))
        return log(z), (Phi*r+2Phibar*r^2+r^3)/z, 3r/z, 3r^2/z
    end
    z = r^3+3Phi*r^2+3Phibar*r+1
    z > 0 || throw(DomainError(z, "Polyakov polynomial is nonpositive"))
    return 3delta+log(z), (Phi*r^2+2Phibar*r+1)/z, 3r^2/z, 3r/z
end

function thermal_level(mass, mu, T, Phi, Phibar, qB, n, ctl::Controls, colors)
    mn2 = mass^2 + 2n*qB
    mn = sqrt(mn2)
    # Every omitted level has Boltzmann weight < exp(-thermal_tail). Keeping
    # n_max fixed avoids changing the equation's Landau budget along PALC.
    maximum_energy = abs(mu)+ctl.thermal_tail*T
    mn >= maximum_energy && return zeros(5)
    # Subtract the exactly integrable T=0 filled sea. Otherwise very small
    # constituent masses make M/E sharply varying near pz=0, even when the
    # Fermi surface itself is resolved. Integrate only the thermal correction.
    kF = min(sqrt(max(mu^2-mn2,0)),ctl.pz_max)
    filled = if kF > 0
        EF = sqrt(mn2+kF^2)
        angle = asinh(kF/mn)
        [ -colors*(abs(mu)*kF-(kF*EF+mn2*angle)/2),
            colors*mass*angle, colors*sign(mu)*kF, 0.0, 0.0]
    else
        zeros(5)
    end
    energy_edges = (max(abs(mu)-ctl.thermal_tail*T, 0),
        max(abs(mu)-12T, 0), max(abs(mu)-ctl.thermal_width*T, 0), abs(mu),
        abs(mu)+ctl.thermal_width*T, abs(mu)+12T, maximum_energy)
    edges = sort!(unique([min(sqrt(max(E^2-mn2,0)),ctl.pz_max)
        for E in energy_edges]))
    values = integrate(edges, ctl.order) do pz
        E = sqrt(mn2+pz^2)
        lp, fp, dp, dbp = occupation((mu-E)/T, Phi, Phibar)
        lm, fm, dbm, dm = occupation((-mu-E)/T, Phibar, Phi)
        occupied = E < abs(mu) ? 1.0 : 0.0
        [ -T*(lp+lm)+colors*(abs(mu)-E)*occupied,
            colors*mass/E*(fp+fm-occupied), colors*(fp-fm-sign(mu)*occupied),
            -T*(dp+dm), -T*(dbp+dbm)]
    end
    return values === nothing ? filled : filled+values
end

function vacuum(mass, qB, p::Parameters, ctl::Controls)
    M = abs(mass)
    M > 0 || throw(ArgumentError("nonzero constituent mass required"))
    L = p.cutoff
    if ctl.route === :mfir
        x = M^2/(2qB)
        base = -p.colors/(8pi^2) * (L*sqrt(L^2+M^2)*(2L^2+M^2) -
            M^4*asinh(L/M))
        mag = -p.colors*qB^2/(2pi^2) * (zeta_derivative(x) -
            (x^2-x)*log(x)/2+x^2/4)
        loop = -p.colors*mass/(2pi^2) * (L*sqrt(L^2+M^2)-M^2*asinh(L/M)) -
            p.colors*mass*qB/(2pi^2) * (loggamma(x)-log(2pi)/2+x-(x-1/2)*log(x))
        return base+mag, loop
    end
    ctl.route === :landau_legacy || throw(ArgumentError("unknown route"))
    total = zeros(2)
    for n in 0:ctl.n_max
        transverse2 = 2n*qB
        boundaries = [sqrt(max((r*L)^2-transverse2,0)) for r in (0.8,1.0,1.2,2.0,4.0)]
        edges = sort!(unique(vcat(0.0, min.(boundaries,ctl.pz_max),ctl.pz_max)))
        value = integrate(edges, ctl.order) do pz
            E = sqrt(M^2+transverse2+pz^2)
            cutoff = 1/(1+((transverse2+pz^2)/L^2)^ctl.cutoff_N)
            [cutoff*E, cutoff*mass/E]
        end
        total .-= p.colors*qB/(2pi^2)*(n==0 ? 1 : 2).*value
    end
    return total[1], total[2]
end

function polyakov(Phi, Phibar, T, p::Parameters)
    a = p.a0+p.a1*p.T0/T+p.a2*(p.T0/T)^2
    b = p.b3*(p.T0/T)^3
    H = 1-6Phi*Phibar+4*(Phi^3+Phibar^3)-3*(Phi*Phibar)^2
    H > 0 || throw(DomainError(H, "nonpositive Haar polynomial"))
    du = T^4*(-a*Phibar/2+b*(-6Phibar+12Phi^2-6Phi*Phibar^2)/H)
    db = T^4*(-a*Phi/2+b*(-6Phi+12Phibar^2-6Phibar*Phi^2)/H)
    return T^4*(-a*Phi*Phibar/2+b*log(H)), du, db
end

function evaluate(x, T, mu, eB, p::Parameters=Parameters(), ctl::Controls=Controls())
    phi = x[1:3]
    Phi, Phibar = x[4:5]
    M = masses(phi,eB,p)
    G = coupling(eB,p)
    chi = 2G*sum(abs2,phi)-4p.K*prod(phi)
    poly, du, db = polyakov(Phi,Phibar,T,p)
    vac, therm = 0.0, 0.0
    loop, rho = zeros(3), zeros(3)
    level_rho = Vector{Float64}[]
    for f in 1:3
        qB = CHARGES[f]*eB
        vf, lf = vacuum(M[f],qB,p,ctl)
        vac += vf
        loop[f] = lf
        densities = Float64[]
        for n in 0:ctl.n_max
            values = thermal_level(M[f],mu[f],T,Phi,Phibar,qB,n,ctl,p.colors)
            values .*= qB/(2pi^2)*(n==0 ? 1 : 2)
            therm += values[1]
            loop[f] += values[2]
            rho[f] += values[3]
            du += values[4]
            db += values[5]
            push!(densities,values[3])
        end
        push!(level_rho,densities)
    end
    gap = loop-phi
    J = mass_jacobian(phi,eB,p)
    gradient = vcat(J'*gap,du,db)
    return (;omega=chi+poly+vac+therm,chi,poly,vac,therm,masses=M,
        condensate=loop,condensate_gap=gap,rho,level_rho,
        residual=vcat(gap,du,db),stationarity=gradient,
        mass_jacobian_sigma_min=minimum(svdvals(J)))
end

end # module
