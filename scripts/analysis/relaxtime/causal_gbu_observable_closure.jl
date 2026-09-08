"""B2/B3 algebra and counterexamples; not a density provider or 2PI solver."""
module CausalGBUObservableClosure
include("causal_gbu_thermal_admissibility.jl")
const A=CausalGBUThermalAdmissibility
const R=A.R
const P=A.P

"""Stable W(delta), including the small-phase cubic cancellation."""
function gbu_weight(delta)
    isfinite(delta) || throw(ArgumentError("finite phase required"))
    if abs(delta)<0.01
        return delta^3*(2/3+delta^2*(-2/15+delta^2*(4/315-2*delta^2/2835)))
    end
    return delta-sin(2delta)/2
end

"""Scalar RPA identity in project units: D0=2K and meson selfenergy=2Pi.

No positivity assumption on ImPi. The supplied phase is locally congruent
to -arg(F); this function never determines bound-state counts or a branch.
Only the REAL-AXIS value has the intended GBU observable interpretation.
"""
function scalar_gbu_parts(pi_value,pi_derivative,K;phase=nothing)
    all(isfinite,(pi_value,pi_derivative,K)) && K>0 ||
        throw(ArgumentError("finite Pi, Pi derivative and positive K required"))
    f=1-4K*pi_value
    abs(f)>0 || throw(ArgumentError("pole requires distributional treatment"))
    delta=phase===nothing ? -angle(f) : phase
    isfinite(delta) && abs(cis(delta)-conj(f)/abs(f))<1e-10 ||
        throw(ArgumentError("phase inconsistent with inverse"))
    d=2K/f
    dd=8K^2*pi_derivative/f^2
    sigma=2pi_value
    dsigma=2pi_derivative
    correction=imag(sigma)*real(d)
    correction_derivative=imag(dsigma)*real(d)+imag(sigma)*real(dd)
    delta_derivative=imag(4K*pi_derivative/f)
    return (phase=delta,selfenergy_correction=correction,
        optical_error=abs(imag(d)-imag(sigma)*abs2(d)),
        weight=gbu_weight(delta),selfenergy_weight=delta-correction,
        weight_error=abs(gbu_weight(delta)-(delta-correction)),
        derivative=2sin(delta)^2*delta_derivative,
        selfenergy_derivative=delta_derivative-correction_derivative,
        derivative_error=abs(2sin(delta)^2*delta_derivative-
            (delta_derivative-correction_derivative)),
        full_product_imaginary=imag(sigma*d),
        born_subtracted_phase=delta-imag(4K*pi_value),
        independent_counting_required=true,production_authorized=false)
end

"""Three determinant pieces share a source, but their RPA logarithms do not add."""
function subtraction_logs(vacuum_small,medium_full,reference_large,K)
    all(isfinite,(vacuum_small,medium_full,reference_large,K)) && K>0 ||
        throw(ArgumentError("invalid subtraction inputs"))
    pieces=(vacuum_small,medium_full,reference_large)
    fs=complex.(1 .-4K.*collect(pieces))
    f=complex(1-4K*(vacuum_small+medium_full-reference_large))
    all(x->abs(x)>0,fs) && abs(f)>0 || throw(ArgumentError("singular logarithm"))
    return (combined=log(f),separately_resummed=log(fs[1])+log(fs[2])-log(fs[3]),
        difference=log(f)-log(fs[1])-log(fs[2])+log(fs[3]),
        physical_equivalence_claimed=false)
end

bose(w,T,nu)=w>nu && T>0 ? inv(expm1((w-nu)/T)) :
    throw(ArgumentError("Bose support requires w>nu and T>0"))

"""A tagged, fixed-profile partial pressure per q shell (q measure omitted).

nu is a BOOKKEEPING derivative at zero, not a fitted fugacity or a changed
quark chemical potential. Using this definition does not prove a stationary
2PI functional, positivity, or the relation to final experimental yields.
"""
function tagged_pressure(phase,edges,T,nu=0.;nodes=64)
    isfinite(T) && T>0 && length(edges)>=2 && all(isfinite,edges) &&
        all(diff(edges).>0) && first(edges)>nu && nodes>=8 ||
        throw(ArgumentError("invalid tagged profile window"))
    total=zero(nu+T)
    for j in 1:length(edges)-1
        xs,ws=R.gauleg(edges[j],edges[j+1],nodes)
        total+=sum(v*bose(w,T,nu)*gbu_weight(phase(w))/pi for (w,v) in zip(xs,ws))
    end
    return total
end

"""Separate explicit Bose and profile-dependence terms, on a fixed window."""
function partial_density_terms(phase,phase_derivative,phase_mu_derivative,edges,T;nodes=64)
    isfinite(T) && T>0 && length(edges)>=2 && all(isfinite,edges) &&
        all(diff(edges).>0) && first(edges)>0 && nodes>=8 ||
        throw(ArgumentError("invalid partial-density window"))
    tag,spectral,response=0.,0.,0.
    for j in 1:length(edges)-1
        xs,ws=R.gauleg(edges[j],edges[j+1],nodes)
        for (w,v) in zip(xs,ws)
            delta=phase(w); g=bose(w,T,0.)
            tag+=v*g*(1+g)/T*gbu_weight(delta)/pi
            spectral+=v*g*2sin(delta)^2*phase_derivative(w)/pi
            response+=v*g*2sin(delta)^2*phase_mu_derivative(w)/pi
        end
    end
    a,b=first(edges),last(edges)
    boundary=(bose(b,T,0.)*gbu_weight(phase(b))-bose(a,T,0.)*gbu_weight(phase(a)))/pi
    return (tagged_density=tag,spectral_density=spectral,
        endpoint_term=boundary,integration_by_parts_error=abs(tag-(spectral-boundary)),
        profile_response=response,total_chemical_derivative=tag+response,
        fixed_profile_is_total_derivative=false,production_authorized=false)
end

"""If F_R>=a>0, |F_I|<=b, then |W|<=2 atan(b/a)^3/3.

Uniformity and the near-zero (no winding) phase branch are the CALLER's
obligation. F alone does not fix that branch. The shell integral bound uses
integral_beta*g*(1+g) = g(lower); it does not discard the tail.
"""
function tail_bound(real_lower,imaginary_upper,omega_lower,T,q;
                    uniform_bounds_certified=false,near_zero_branch_certified=false)
    all(isfinite,(real_lower,imaginary_upper,omega_lower,T,q)) &&
        min(real_lower,omega_lower,T)>0 && min(imaginary_upper,q)>=0 ||
        throw(ArgumentError("invalid tail bound inputs"))
    epsilon=atan(imaginary_upper/real_lower)
    return (phase_bound=epsilon,weight_bound=2epsilon^3/3,
        shell_bound_inv_fm2=q^2/(2pi^3)*(2epsilon^3/3)*bose(omega_lower,T,0.),
        rigorous_bound_certified=uniform_bounds_certified && near_zero_branch_certified,
        production_authorized=false)
end

"""Exterior RPA zero exclusion for Pi(z)=integral rho(x)/(x-z)/pi.

For support in [-S,S], integral |rho| <= M, |z|>=R>S, the triangle
inequality gives |4K Pi|<=4K M/(pi*(R-S)). No positivity is assumed.
Certification requires a genuine total-variation upper bound from caller;
quadrature estimates alone must not set mass_bound_certified=true.
"""
function exterior_exclusion(support_radius,total_variation,K,radius;mass_bound_certified=false)
    all(isfinite,(support_radius,total_variation,K,radius)) &&
        min(support_radius,total_variation)>=0 && K>0 && radius>support_radius ||
        throw(ArgumentError("finite support, total variation, K>0 and R>S required"))
    bubble_bound=total_variation/(pi*(radius-support_radius))
    margin=1-4K*bubble_bound
    return (bubble_bound_inv_fm2=bubble_bound,inverse_lower_bound=margin,
        exterior_zero_free_certified=mass_bound_certified && margin>0,
        production_authorized=false)
end

"""Exactly soluble signed auxiliary determinant F=(z^2-E^2)/(z^2-P^2).

E,P>0 implies no UHP zeros/poles; E>P has negative positive-frequency Pi
spectral weight and negative signed excess count. Not a physical meson.
"""
toy_inverse(z,E2,P2)=(z^2-E2)/(z^2-P2)

function toy_matsubara(E,P,T;terms=512)
    all(isfinite,(E,P,T)) && min(E,P,T)>0 && terms>=16 ||
        throw(ArgumentError("stable toy requires positive energies/T"))
    a,b=(E/(2pi*T))^2,(P/(2pi*T))^2
    raw=T/2*log(E^2/P^2)+T*sum(log1p(a/n^2)-log1p(b/n^2) for n in 1:terms)
    # Analytic leading tail, with a bound on the omitted log expansion.
    zeta2tail=pi^2/6-sum(inv(Float64(n)^2) for n in 1:terms)
    value=raw+T*(a-b)*zeta2tail
    remainder_bound=T*(a^2+b^2)/(6terms^3)
    exact=(E-P)/2+T*(log(-expm1(-E/T))-log(-expm1(-P/T)))
    return (matsubara=value,exact=exact,remainder_bound=remainder_bound,
        thermal=exact-(E-P)/2,signed_excess_count=bose(E,T,0.)-bose(P,T,0.),
        positive_frequency_pi_weight_sign=sign(P^2-E^2),is_project_result=false)
end
end
