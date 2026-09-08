"""Ordered P/S bubble from one regulated two-line spectral sum.

This diagnostic has its own real and imaginary parts and matching contact
term. It never imports the retained PV/log continuation or an external A.
Finite quadrature atoms are not physical poles. Only upper-half-plane values
and values in a conservative open analytic gap may be evaluated directly.
"""
module CausalSpectralBubble

using ..GaussLegendre: gauleg
using ..OneLoopIntegrals: B0_spectral_cut, _spectral_two_line_support
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.Constants_PNJL: N_color, Λ_inv_fm

export SpectralBubbleGrid, build_spectral_bubble, spectral_bubble, spectral_bubble_cut
export PiecewiseSpectralFunction, cauchy_transform, build_bubble_dispersion

"""Compact, piecewise-linear Im Pi profile; endpoints must vanish.

Its Cauchy transform is exact for the interpolant. Mesh convergence remains
separate from the eta limit. No density positivity or phase fold is imposed.
"""
struct PiecewiseSpectralFunction
    energy::Vector{Float64}
    imaginary::Vector{Float64}
    function PiecewiseSpectralFunction(energy,imaginary)
        x,y = Float64.(energy),Float64.(imaginary)
        length(x) == length(y) && length(x) >= 3 || throw(ArgumentError("spectral arrays need equal length >=3"))
        all(isfinite,x) && all(isfinite,y) && all(diff(x) .> 0) ||
            throw(ArgumentError("spectral nodes must be finite and strictly increasing"))
        first(y) == last(y) == 0 || throw(ArgumentError("compact spectral endpoints must be zero"))
        new(x,y)
    end
end

# Exact linear-cell integral written without slope * distance cancellation.
# H(t)=1-log(1+t)/t = t/2-t^2/3+...; the series handles remote/narrow cells.
function _cauchy_cell_moment(t,logarithm)
    if abs(t)<0.01
        term=t
        value=zero(t)
        for n in 1:8
            value+=term/(n+1)
            term*=-t
        end
        return value
    end
    return 1-logarithm/t
end

function cauchy_transform(profile::PiecewiseSpectralFunction,z::Number)
    z = ComplexF64(z)
    isfinite(real(z)) && isfinite(imag(z)) && imag(z) >= 0 ||
        throw(ArgumentError("Cauchy energy must be finite in the closed upper half-plane"))
    x,y = profile.energy,profile.imaginary
    if imag(z) > 0
        value = 0.0im
        for i in 1:length(x)-1
            a,b = x[i],x[i+1]
            t=(b-a)/(a-z)
            logarithm=log1p(t)
            value += y[i]*logarithm+(y[i+1]-y[i])*_cauchy_cell_moment(t,logarithm)
        end
        return value/π
    end
    w = real(z)
    j = searchsortedlast(x,w)
    rho = 1 <= j < length(x) ? y[j]+(y[j+1]-y[j])*(w-x[j])/(x[j+1]-x[j]) : 0.0
    value = 0.0
    for i in 1:length(x)-1
        a,b = x[i],x[i+1]
        # At a knot the subtraction coefficient is exactly zero analytically.
        if w != a && w != b
            t=(b-a)/(a-w)
            # For a point inside the cell t<-1: take the real PV logarithm.
            logarithm=t > -1 ? log1p(t) : log(abs((b-w)/(a-w)))
            value += (y[i]-rho)*logarithm+(y[i+1]-y[i])*_cauchy_cell_moment(t,logarithm)
        else
            value += y[i+1]-y[i]
        end
    end
    if rho != 0
        value += rho*log(abs((last(x)-w)/(first(x)-w)))
    end
    return complex(value/π,rho)
end

struct SpectralBubbleGrid
    poles_inv_fm::Vector{Float64}
    b0_residues_inv_fm::Vector{Float64}
    pi_residues_inv_fm3::Vector{Float64}
    contact_inv_fm2::Float64
    q_inv_fm::Float64
    m1_inv_fm::Float64
    m2_inv_fm::Float64
    mu1_inv_fm::Float64
    mu2_inv_fm::Float64
    T_inv_fm::Float64
    Phi::Float64
    PhiBar::Float64
    mass_term_inv_fm2::Float64
    vacuum_cutoff_inv_fm::Float64
    thermal_cutoff_inv_fm::Float64
    momentum_nodes::Int
    angle_nodes::Int
end

function _append_component!(poles, bweights, piweights, q, m1, mu1, m2, mu2,
                            T, phi, phibar, mass_term, cutoff, np, nx, component)
    q >= 2cutoff && return 0.0
    # Split at the overlap boundary of the two cutoff spheres.
    edges = sort!(unique([max(0.0,q-cutoff),abs(cutoff-q),cutoff]))
    contact = 0.0
    for j in 1:length(edges)-1
        left,right = edges[j],edges[j+1]
        left < right || continue
        ps, pws = gauleg(left,right,np)
        for (p,pw) in zip(ps,pws)
            E1 = hypot(p,m1)
            xmin = q == 0 ? -1.0 : max(-1.0,(p^2+q^2-cutoff^2)/(2p*q))
            xs,xws = q == 0 ? ([0.0],[2.0]) : gauleg(xmin,1.0,nx)
            for (x,xw) in zip(xs,xws)
                E2 = sqrt(p^2+q^2-2p*q*x+m2^2)
                n1 = component === :vacuum ? (1.0,0.0) :
                    (-antiquark_distribution(E1,mu1,T,phi,phibar),quark_distribution(E1,mu1,T,phi,phibar))
                n2 = component === :vacuum ? (1.0,0.0) :
                    (-antiquark_distribution(E2,mu2,T,phi,phibar),quark_distribution(E2,mu2,T,phi,phibar))
                measure = pw*xw*p^2/(E1*E2)
                # Sum R_st*u analytically. It is the matching A1+A2 at q=0.
                contact += measure*2*(E1*(n2[2]-n2[1])+E2*(n1[2]-n1[1]))
                for (si,s) in enumerate((-1,1)), (ti,t) in enumerate((-1,1))
                    u = s*E1-t*E2
                    residue = measure*s*t*(n2[ti]-n1[si])
                    residue == 0 && continue
                    push!(poles,u)
                    push!(bweights,residue)
                    # On-shell Dirac trace, not the external z-dependent prefactor.
                    trace_weight = u^2-q^2-mass_term
                    push!(piweights,N_color/(8π^2)*trace_weight*residue)
                end
            end
        end
    end
    return contact
end

"""Build a solver-free ordered bubble grid (all energies in fm^-1).

Both vacuum line momenta are below `vacuum_cutoff_inv_fm`. Both thermal
line momenta are below `thermal_cutoff_inv_fm` (default: same cutoff).
Taking a larger thermal cutoff is an explicit vacuum/thermal split, not
a change to PNJLCore. Its infinite thermal-cutoff limit is not certified.
P uses (m1-m2)^2, S uses (m1+m2)^2. Phi/PhiBar must be in [0,1].
"""
function build_spectral_bubble(q::Real,m1::Real,mu1::Real,m2::Real,mu2::Real,T::Real;
        channel::Symbol=:P, Phi::Real=0.0, PhiBar::Real=0.0,
        vacuum_cutoff_inv_fm::Real=Λ_inv_fm,
        thermal_cutoff_inv_fm::Real=vacuum_cutoff_inv_fm,
        momentum_nodes::Integer=64,angle_nodes::Integer=32)
    q,a,u,b,v,temp,phi,phibar,vc,tc = Float64.((q,m1,mu1,m2,mu2,T,Phi,PhiBar,
        vacuum_cutoff_inv_fm,thermal_cutoff_inv_fm))
    all(isfinite,(q,a,u,b,v,temp,phi,phibar,vc,tc)) || throw(ArgumentError("bubble inputs must be finite"))
    q >= 0 && a > 0 && b > 0 && temp > 0 && vc > 0 && tc >= vc ||
        throw(ArgumentError("requires q>=0, positive masses/T/cutoffs and thermal_cutoff>=vacuum_cutoff"))
    0 <= phi <= 1 && 0 <= phibar <= 1 || throw(ArgumentError("Phi/PhiBar must be in [0,1]"))
    channel in (:P,:S) || throw(ArgumentError("channel must be :P or :S"))
    momentum_nodes >= 4 && angle_nodes >= 4 || throw(ArgumentError("quadrature nodes must be at least 4"))
    mass_term = channel === :P ? (a-b)^2 : (a+b)^2
    poles,bweights,piweights = Float64[],Float64[],Float64[]
    contact = 0.0
    for (component,cutoff) in ((:vacuum,vc),(:thermal,tc))
        contact += _append_component!(poles,bweights,piweights,q,a,u,b,v,temp,phi,phibar,
            mass_term,cutoff,Int(momentum_nodes),Int(angle_nodes),component)
    end
    return SpectralBubbleGrid(poles,bweights,piweights,contact,q,a,b,u,v,temp,phi,phibar,
        mass_term,vc,tc,Int(momentum_nodes),Int(angle_nodes))
end

"""Evaluate the full Pi and B0 at external k0+i*eta.

eta=0 is allowed only strictly between the uncut Landau and unitary
envelopes (either frequency sign). No naive real-axis quadrature on cuts,
static-limit choice, second sheet, or production certification is implied.
"""
function spectral_bubble(grid::SpectralBubbleGrid,k0::Real;eta_inv_fm::Real=1e-2)
    w,eta = Float64(k0),Float64(eta_inv_fm)
    isfinite(w) && isfinite(eta) && eta >= 0 || throw(ArgumentError("k0 must be finite and eta nonnegative"))
    lambda = w+grid.mu1_inv_fm-grid.mu2_inv_fm
    if eta == 0
        lo = hypot(grid.q_inv_fm,grid.m1_inv_fm-grid.m2_inv_fm)
        hi = hypot(grid.q_inv_fm,grid.m1_inv_fm+grid.m2_inv_fm)
        lo < abs(lambda) < hi || throw(ArgumentError("eta=0 evaluation requires an open analytic gap"))
    end
    z = complex(lambda,eta)
    b0,value,derivative = 0.0im,0.0im,0.0im
    for i in eachindex(grid.poles_inv_fm)
        denominator = z-grid.poles_inv_fm[i]
        b0 += grid.b0_residues_inv_fm[i]/denominator
        value += grid.pi_residues_inv_fm3[i]/denominator
        derivative -= grid.pi_residues_inv_fm3[i]/denominator^2
    end
    reconstructed = N_color/(8π^2)*((z^2-grid.q_inv_fm^2-grid.mass_term_inv_fm2)*b0-grid.contact_inv_fm2)
    all(isfinite,(real(b0),imag(b0),real(value),imag(value),real(derivative),imag(derivative))) ||
        throw(ArgumentError("spectral bubble result is not finite"))
    return (value=value,B0=b0,derivative=derivative,contact_inv_fm2=grid.contact_inv_fm2,
        reconstructed_value=reconstructed,contact_identity_residual=abs(value-reconstructed),
        lambda_inv_fm=lambda,eta_inv_fm=eta,analytic_scope=eta>0 ? :upper_half_plane : :open_gap,
        regulator=:two_line_vacuum_thermal,physical_cut_certified=false,production_authorized=false)
end

"""Exact on-shell cut of the SAME regulated bubble; does not supply a PV real part."""
function spectral_bubble_cut(grid::SpectralBubbleGrid,k0::Real;energy_nodes::Integer=64)
    lambda = Float64(k0)+grid.mu1_inv_fm-grid.mu2_inv_fm
    return _spectral_bubble_cut_lambda(grid,lambda;energy_nodes=energy_nodes)
end

function _spectral_bubble_cut_lambda(grid::SpectralBubbleGrid,lambda;energy_nodes::Integer=64)
    pair,landau = 0.0,0.0
    unresolved = false
    for (component,cutoff) in ((:vacuum,grid.vacuum_cutoff_inv_fm),(:thermal,grid.thermal_cutoff_inv_fm))
        cut = B0_spectral_cut(lambda,grid.q_inv_fm,grid.m1_inv_fm,grid.mu1_inv_fm,
            grid.m2_inv_fm,grid.mu2_inv_fm,grid.T_inv_fm;Φ=grid.Phi,Φbar=grid.PhiBar,
            pmax_inv_fm=cutoff,energy_nodes=energy_nodes,component=component)
        pair += cut.pair
        landau += cut.landau
        unresolved |= cut.static_degeneracy_unresolved
    end
    factor = N_color/(8π^2)*(lambda^2-grid.q_inv_fm^2-grid.mass_term_inv_fm2)
    return (imaginary=factor*(pair+landau),pair=factor*pair,landau=factor*landau,
        B0_imaginary=pair+landau,static_degeneracy_unresolved=unresolved,
        production_authorized=false)
end

"""Build a real-axis spectral interpolant, using internal lambda coordinates.

Includes kinematic/cutoff knots and the exact zero-external-frequency point.
Cosine-spaced subinterval nodes resolve threshold square roots. At q=0 a
hard-cutoff jump is approximated by the last mesh cell; refine the mesh to
test its integrated effect. This is not an error-controlled PV solver.
"""
function build_bubble_dispersion(grid::SpectralBubbleGrid;
        segment_nodes::Integer=64,energy_nodes::Integer=64)
    segment_nodes >= 8 || throw(ArgumentError("segment_nodes must be at least 8"))
    q,a,b = grid.q_inv_fm,grid.m1_inv_fm,grid.m2_inv_fm
    shift = grid.mu1_inv_fm-grid.mu2_inv_fm
    upper = hypot(a,grid.thermal_cutoff_inv_fm)+hypot(b,grid.thermal_cutoff_inv_fm)
    edges = [-upper,upper,0.0]
    abs(shift) < upper && push!(edges,shift)
    for e in (hypot(q,a+b),hypot(q,a-b))
        e < upper && append!(edges,[-e,e])
    end
    for cutoff in (grid.vacuum_cutoff_inv_fm,grid.thermal_cutoff_inv_fm)
        q < 2cutoff || continue
        support = _spectral_two_line_support(q,a,b,cutoff)
        for e in (support.pair...,support.landau...)
            -upper < e < upper && append!(edges,[-e,e])
        end
        for (p,r) in ((cutoff,cutoff),(cutoff,abs(cutoff-q)),(abs(cutoff-q),cutoff))
            for s in (-1,1),t in (-1,1)
                e = s*hypot(a,p)-t*hypot(b,r)
                -upper < e < upper && push!(edges,e)
            end
        end
    end
    sort!(unique!(edges))
    x = Float64[]
    for j in 1:length(edges)-1
        left,right = edges[j],edges[j+1]
        append!(x,[left+(right-left)*(1-cospi(i/segment_nodes))/2 for i in 0:segment_nodes-1])
    end
    push!(x,upper)
    sort!(unique!(x))
    y = [_spectral_bubble_cut_lambda(grid,l;energy_nodes=energy_nodes).imaginary for l in x]
    # Compact-support endpoints have zero measure; no interior values are clipped.
    y[1] = y[end] = 0.0
    return PiecewiseSpectralFunction(x,y)
end

end # module CausalSpectralBubble
