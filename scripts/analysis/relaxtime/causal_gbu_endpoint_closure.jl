"""Analysis-only endpoint obstruction to finite-cutoff all-root counting."""
module CausalGBUEndpointClosure
include("causal_gbu_observable_closure.jl")
const O=CausalGBUObservableClosure
const A=O.A
const R=O.R
const P=O.P

"""Exact equal-ball intersection volume (fm^-3)."""
function lens_volume(q,L)
    all(isfinite,(q,L)) && q>=0 && L>0 || throw(ArgumentError("q>=0, L>0 required"))
    return q>=2L ? 0. : pi*(4L+q)*(2L-q)^2/12
end

"""Analytic bound on integral |rho|/pi, without spectral interpolation.

Projected P/S Dirac traces obey T_st>=0, sum_st T_st=4. Full/vacuum
occupation differences have modulus <=1. Triangle inequality for
vacuum_Lambda+full_Lth-vacuum_Lth gives the bound. This real-arithmetic
theorem does not claim directed-rounding certification of floating values.
"""
function projected_exterior(q,a,b,vc,tc,K;Nc=3)
    all(isfinite,(q,a,b,vc,tc,K,Nc)) && q>=0 && min(a,b,vc,K,Nc)>0 && tc>=vc ||
        throw(ArgumentError("invalid projected bound parameters"))
    mass=Nc*(lens_volume(q,vc)+2lens_volume(q,tc))/(2pi^3)
    support=hypot(a,tc)+hypot(b,tc)
    radius=support+max(1.,8K*mass)
    return (support_radius_inv_fm=support,total_variation_over_pi_bound_inv_fm3=mass,
        exterior_radius_inv_fm=radius,inverse_deviation_bound=4K*mass/(radius-support),
        analytic_exterior_exclusion=true,directed_rounding_certified=false)
end

"""q0 pair spectral endpoint approached from INSIDE a finite thermal cutoff.

For tc>vc its weight is -(n1+nbar2), strictly negative at finite T before
floating underflow. Assigning the isolated endpoint zero does not change it.
"""
function thermal_endpoint(a,u,b,v,T,phi,bar,vc,tc;Nc=3)
    P.validate(0.,a,u,b,v,T,phi,bar,vc,4)
    isfinite(tc) && tc>vc && Nc>0 || throw(ArgumentError("strict thermal extension required"))
    e1,e2=hypot(tc,a),hypot(tc,b)
    n1,n2=P.occupations(e1,u,T,phi,bar),P.occupations(e2,v,T,phi,bar)
    upper=e1+e2
    rho=-Nc/(8pi^2)*(upper^2-(a-b)^2)*(2pi*tc/upper)*(n1.quark+n2.anti)
    return (lambda_endpoint_inv_fm=upper,k0_endpoint_inv_fm=upper-u+v,
        rho_inside_inv_fm2=rho,occupation_sum=n1.quark+n2.anti,
        nonzero_in_float64=rho<0,positive_temperature_negative_endpoint=true)
end

"""c=rho(S-)<0 and finite Pi_reg(S)=h imply F(S+0)=-Inf, F(Inf)=1.

This proves existence, not uniqueness. The log-distance is an asymptotic
estimate, not a root residual. Neither a UHP pole nor a physical bound
state is inferred from this intermediate-value argument.
"""
function endpoint_obstruction(c,S,width,h,K,T,shift)
    all(isfinite,(c,S,width,h,K,T,shift)) && c<0 && S>width>0 && K>0 && T>0 && S>shift ||
        throw(ArgumentError("negative endpoint, positive width/K/T and Bose support required"))
    logdistance=log(width)+pi*(1-4K*h)/(4K*c)
    return (at_least_one_positive_exterior_zero=true,
        asymptotic_log_distance_inv_fm=logdistance,
        asymptotic_log10_distance_inv_fm=logdistance/log(10.),
        sub_float64_resolution=logdistance<log(eps(S)),
        upper_bose_weight_per_zero=inv(expm1((S-shift)/T)),
        bose_bound_is_total_density_bound=false,uniqueness_proved=false,
        distance_is_asymptotic=true,UHP_instability_implied=false,
        production_authorized=false)
end

"""Soluble signed band: rho(x)=-c*sign(x) on a<|x|<S, c>0."""
function toy_inverse(z,a,S,c,K)
    0<a<S && c>0 && K>0 || throw(ArgumentError("invalid signed-band toy"))
    z=complex(z)
    return 1+4K*c/pi*(log(S-z)+log(S+z)-log(a-z)-log(a+z))
end

function toy_exterior_root(a,S,c,K)
    0<a<S && c>0 && K>0 || throw(ArgumentError("invalid signed-band toy"))
    t=exp(-pi/(4K*c))
    root=sqrt((S^2-t*a^2)/(1-t))
    slope=4K*c/pi*(2root/(root^2-S^2)-2root/(root^2-a^2))
    return (root=root,residue=-2K/slope,phase_jump=pi,
        is_positive_physical_bound_state=false,is_synthetic=true)
end

"""Dependency gate: running a check is not the same as passing its physics."""
function stage_status(;analytic_structure,physical_counting,integral_convergence)
    return (step1=analytic_structure ? "passed" : "failed",
        step2=!analytic_structure ? "blocked_by_step1" : physical_counting ? "passed" : "failed",
        step3=!(analytic_structure && physical_counting) ? "blocked_by_counting" :
            integral_convergence ? "passed" : "failed",
        research_production_ready=analytic_structure && physical_counting && integral_convergence,
        production_authorized=false)
end
end
