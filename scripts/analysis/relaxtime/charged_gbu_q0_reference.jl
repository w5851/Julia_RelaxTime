"""Analysis-only timelike q=0 extrapolation of the current infinite-thermal GBU kernel.

This is the retained internal-lambda reference approximation, not a finite-q
retarded loop or an eta continuation. Bose weights always use external omega.
The spacelike/negative-lambda sector is outside this approximation (phase zero).
The default density is a finite-window derivative integral. A nonzero boosted
static imaginary part is recorded, while the historical zero-limit prescription
is available only as an explicit strict comparison.
"""
module ChargedGBUQ0Reference

const Workflow = Main.ChargedGBUResearchWorkflow
const P = Workflow.P
const Y = Workflow.Y
const ENDPOINT_POLICIES = ("finite_window", "strict_zero_limit")

function validate_prescription(policy, omega_lower_inv_fm)
    policy in ENDPOINT_POLICIES || throw(ArgumentError("unsupported q0 endpoint_policy: $(policy)"))
    isfinite(omega_lower_inv_fm) && omega_lower_inv_fm > 0 ||
        throw(ArgumentError("q0 omega_lower_inv_fm must be positive and finite"))
    return nothing
end

function reference_coordinate(omega_inv_fm, q_inv_fm, shift_inv_fm)
    all(isfinite, (omega_inv_fm, q_inv_fm, shift_inv_fm)) && q_inv_fm >= 0 ||
        throw(ArgumentError("finite coordinates and nonnegative q_inv_fm required"))
    lambda = omega_inv_fm + shift_inv_fm
    lambda >= q_inv_fm || return nothing
    # Factorization avoids cancellation at the light-cone onset.
    return sqrt((lambda - q_inv_fm) * (lambda + q_inv_fm))
end

boosted_energy(lambda_inv_fm, q_inv_fm, shift_inv_fm) =
    hypot(lambda_inv_fm, q_inv_fm) - shift_inv_fm

struct Reference{T,R}
    profile::T
    rootsets::R
end

function reference(profile)
    profile.kernel.q == 0 || throw(ArgumentError("reference requires a q=0 kernel"))
    origin = P.inverse(profile, 0.)
    real(origin) > 0 && abs(imag(origin)) < 1e-8 ||
        throw(ArgumentError("q0 reference onset is not positive real"))
    return Reference(profile, Y.gap_roots(profile))
end

function cut_phase(ref::Reference, omega_inv_fm, q_inv_fm)
    p = ref.profile
    lambda0 = reference_coordinate(omega_inv_fm, q_inv_fm, p.kernel.shift)
    lambda0 === nothing && return 0.
    f = P.inverse(p, lambda0)
    abs(f) > 1e-12 || throw(ArgumentError("unresolved zero in q0 cut integration"))
    return -angle(f)
end

function cut_integral(ref::Reference, lo, hi, q_inv_fm, T; nodes=64)
    lo < hi || return 0.
    p = ref.profile; shift = p.kernel.shift
    source_edges = vcat(filter(>=(0), p.kernel.edges), [0., p.threshold.U + p.threshold.H])
    mapped_edges = [boosted_energy(x, q_inv_fm, shift) for x in source_edges]
    edges = sort!(unique!(vcat([lo, hi], filter(x -> lo < x < hi, mapped_edges))))
    # Integration remains in external omega: no hidden coordinate Jacobian.
    return real(Y.I.C.mapped_integral(edges; nodes=nodes) do omega
        g = inv(expm1(omega / T))
        g * (1 + g) / T * Y.O.gbu_weight(cut_phase(ref, omega, q_inv_fm)) / pi
    end)
end

"""W(delta)-W(anchor), without subtracting nearly equal small cubic weights."""
function weight_difference(delta, anchor)
    all(isfinite, (delta, anchor)) || throw(ArgumentError("finite phases required"))
    h = delta - anchor
    h_minus_sin = abs(h) < .1 ? h^3 * (1/6 + h^2 * (-1/120 +
        h^2 * (1/5040 + h^2 * (-1/362880 + h^2/39916800)))) : h - sin(h)
    return 2h * sin((delta + anchor)/2)^2 + cos(delta + anchor) * h_minus_sin
end

"""Integrate g W'/pi on [lo,hi] with both boundary terms, without differentiating phase.

Subtracting the constant W(lo) before integration avoids cancellation between
the infrared bulk and its lower boundary. This is exactly integration by parts
on the finite window, not a claim about the limit lo -> 0.
"""
function finite_window_integral(phase, lo, hi, T; nodes=64, edges=Float64[],
        phase_lo=nothing, phase_hi=nothing)
    all(isfinite, (lo, hi, T)) && 0 < lo < hi && T > 0 && nodes >= 8 ||
        throw(ArgumentError("finite ordered positive frequency window, T>0 and nodes>=8 required"))
    a = phase_lo === nothing ? phase(lo) : phase_lo
    b = phase_hi === nothing ? phase(hi) : phase_hi
    # Resolve each logarithmic decade explicitly; the lower bound is never an
    # implicit consequence of the first Gauss node.
    panels = max(1, ceil(Int, log(hi/lo)/log(4.)))
    log_edges = exp.(range(log(lo), log(hi); length=panels+1))
    log_edges[1] = lo; log_edges[end] = hi
    breaks = sort!(unique!(vcat([lo, hi], filter(x -> lo < x < hi, vcat(edges, log_edges)))))
    value = real(Y.I.C.mapped_integral(breaks; nodes=nodes) do omega
        g = inv(expm1(omega/T))
        g*(1+g)/T * weight_difference(phase(omega), a)/pi
    end) + inv(expm1(hi/T))*weight_difference(b, a)/pi
    lower_boundary = inv(expm1(lo/T))*Y.O.gbu_weight(a)/pi
    upper_boundary = inv(expm1(hi/T))*Y.O.gbu_weight(b)/pi
    return (derivative=value, lower_boundary=lower_boundary, upper_boundary=upper_boundary)
end

function finite_cut_integral(ref::Reference, lo, hi, q_inv_fm, T; nodes=64,
        phase_lo=nothing, phase_hi=nothing)
    lo < hi || return (derivative=0., lower_boundary=0., upper_boundary=0.)
    p = ref.profile; shift = p.kernel.shift
    source_edges = vcat(filter(>=(0), p.kernel.edges), [0., p.threshold.U + p.threshold.H])
    edges = [boosted_energy(x, q_inv_fm, shift) for x in source_edges]
    return finite_window_integral(w -> cut_phase(ref, w, q_inv_fm), lo, hi, T;
        nodes=nodes, edges=edges, phase_lo=phase_lo, phase_hi=phase_hi)
end

"""Boost roots and integrate the signed GBU derivative with explicit endpoint policy."""
function shell(ref::Reference, q_inv_fm; nodes=64, upper=nothing,
        endpoint_policy="finite_window", omega_lower_inv_fm=1e-5)
    validate_prescription(endpoint_policy, omega_lower_inv_fm)
    isfinite(q_inv_fm) && q_inv_fm >= 0 || throw(ArgumentError("nonnegative finite q_inv_fm required"))
    p = ref.profile; k = p.kernel; bg = k.bg; T = bg.T; shift = k.shift
    upper = upper === nothing ? max(24., q_inv_fm + 24T) : Float64(upper)
    unitary = hypot(k.threshold, q_inv_fm)
    landau_edge = hypot(k.landau, q_inv_fm)
    lower = max(0., landau_edge - shift); threshold = unitary - shift
    abs(shift) < unitary && threshold < upper < hypot(k.split, q_inv_fm) - abs(shift) ||
        throw(ArgumentError("Bose-safe q0 reference geometry required"))
    omega_lower_inv_fm < upper || throw(ArgumentError("omega_lower_inv_fm must be below upper"))
    static_coordinate = reference_coordinate(0., q_inv_fm, shift)
    f0 = static_coordinate === nothing ? complex(1.) : P.inverse(p, static_coordinate)
    isfinite(f0) && real(f0) > 0 ||
        throw(ArgumentError("q0 extrapolated static inverse is not finite positive real-part"))
    endpoint_warning = abs(imag(f0)) >= 1e-8
    endpoint_policy == "strict_zero_limit" && endpoint_warning &&
        throw(ArgumentError("static instability or unresolved Bose endpoint"))
    negative, positive = ref.rootsets
    length(negative) <= 1 && length(positive) <= 1 ||
        throw(ArgumentError("q0 reference requires reviewed extra-root topology"))
    roots = [boosted_energy(x, q_inv_fm, shift) for x in positive]
    all(>(0), roots) || throw(ArgumentError("nonpositive-energy root"))
    n = length(roots)
    landau_end = lower > 1e-8 ? cut_phase(ref, lower - 1e-8, q_inv_fm) : 0.
    # The coordinate map multiplies the sqrt-cut coefficient by a positive
    # sqrt(unitary/k.threshold); only its sign enters this independent limit.
    pair_start = Y.Limits.threshold_limit(P.inverse(p, k.threshold), p.threshold.positive;
        inverse_error_budget=4bg.coupling[k.ch] * 1e-6).phase
    abs(landau_end) < 0.005 && abs(pair_start - n * pi) < 0.005 ||
        throw(ArgumentError("independent q0 root count and cut phase limits disagree"))
    highphase = cut_phase(ref, upper, q_inv_fm)
    abs(highphase) < 0.005 || throw(ArgumentError("UV phase not near zero"))
    landau, pair, bound, lower_boundary, upper_boundary = if endpoint_policy == "finite_window"
        l = finite_cut_integral(ref, omega_lower_inv_fm, lower, q_inv_fm, T;
            nodes=nodes, phase_hi=0.)
        pair_lower = max(omega_lower_inv_fm, threshold)
        u = finite_cut_integral(ref, pair_lower, upper, q_inv_fm, T;
            nodes=nodes, phase_lo=pair_lower == threshold ? pair_start : nothing, phase_hi=highphase)
        b = sum(x -> inv(expm1(x/T)), filter(x -> omega_lower_inv_fm < x < upper, roots); init=0.)
        (l.derivative, u.derivative, b, l.lower_boundary+u.lower_boundary,
            l.upper_boundary+u.upper_boundary)
    else
        l = cut_integral(ref, 0., lower, q_inv_fm, T; nodes=nodes)
        u = cut_integral(ref, threshold, upper, q_inv_fm, T; nodes=nodes)
        b = sum(x -> inv(expm1(x/T)), roots; init=0.)
        boundary = n * inv(expm1(threshold/T))
        (l, u-boundary, b, boundary, 0.)
    end
    measure = q_inv_fm^2 / (2pi^2)
    return (density=measure * (bound + landau + pair),
        bound=measure * bound, landau=measure * landau, pair=measure * pair,
        root_count=n, negative_root_count=length(negative), roots=roots,
        landau_phase=landau_end, threshold_phase=pair_start, high_phase=highphase,
        omega_tail_conditional_bound=measure * inv(expm1(upper / T)),
        static_inverse=real(f0), static_inverse_imag=imag(f0), static_phase=-angle(f0),
        static_gbu_weight=Y.O.gbu_weight(-angle(f0)), endpoint_warning=endpoint_warning,
        warning_code=endpoint_warning ? "q0_extrapolated_static_imaginary" : "",
        endpoint_policy=endpoint_policy,
        omega_lower_inv_fm=endpoint_policy == "finite_window" ? Float64(omega_lower_inv_fm) : 0.,
        omega_upper_inv_fm=upper, lower_boundary=measure*lower_boundary,
        upper_boundary=measure*upper_boundary,
        production_authorized=false, route="q0_lambda_reference")
end

end
