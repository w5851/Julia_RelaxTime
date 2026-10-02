"""Analysis-only timelike q=0 extrapolation of the current infinite-thermal GBU kernel.

This is the retained internal-lambda reference approximation, not a finite-q
retarded loop or an eta continuation. Bose weights always use external omega.
The spacelike/negative-lambda sector is outside this approximation (phase zero).
"""
module ChargedGBUQ0Reference

const Workflow = Main.ChargedGBUResearchWorkflow
const P = Workflow.P
const Y = Workflow.Y

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

"""Boost independently counted roots; retain signed continuum and its boundary term."""
function shell(ref::Reference, q_inv_fm; nodes=64, upper=nothing)
    isfinite(q_inv_fm) && q_inv_fm >= 0 || throw(ArgumentError("nonnegative finite q_inv_fm required"))
    p = ref.profile; k = p.kernel; bg = k.bg; T = bg.T; shift = k.shift
    upper = upper === nothing ? max(24., q_inv_fm + 24T) : Float64(upper)
    unitary = hypot(k.threshold, q_inv_fm)
    landau_edge = hypot(k.landau, q_inv_fm)
    lower = max(0., landau_edge - shift); threshold = unitary - shift
    abs(shift) < unitary && threshold < upper < hypot(k.split, q_inv_fm) - abs(shift) ||
        throw(ArgumentError("Bose-safe q0 reference geometry required"))
    static_coordinate = reference_coordinate(0., q_inv_fm, shift)
    f0 = static_coordinate === nothing ? complex(1.) : P.inverse(p, static_coordinate)
    real(f0) > 0 && abs(imag(f0)) < 1e-8 ||
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
    landau = cut_integral(ref, 0., lower, q_inv_fm, T; nodes=nodes)
    pair = cut_integral(ref, threshold, upper, q_inv_fm, T; nodes=nodes)
    bound = sum(x -> inv(expm1(x / T)), roots; init=0.)
    boundary = n * inv(expm1(threshold / T))
    measure = q_inv_fm^2 / (2pi^2)
    return (density=measure * (bound + landau + pair - boundary),
        bound=measure * bound, landau=measure * landau, pair=measure * (pair - boundary),
        root_count=n, negative_root_count=length(negative), roots=roots,
        landau_phase=landau_end, threshold_phase=pair_start, high_phase=highphase,
        omega_tail_conditional_bound=measure * inv(expm1(upper / T)),
        static_inverse=real(f0), production_authorized=false, route="q0_lambda_reference")
end

end
