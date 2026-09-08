"""
    BUPhaseGates

Pure contracts for the real-axis Beth-Uhlenbeck energy measure, high-energy
phase anchoring, subthreshold root counting, and Levinson/Mott acceptance
gates. The module does not evaluate a propagator or change any production
default.
"""
module BUPhaseGates

export STRICT_SINGLE_CHARGE_OMEGA_MEASURE, LEGACY_POSITIVE_ENERGY_OMEGA_MEASURE
export bu_omega_measure, bu_omega_measure_factor
export anchor_phase_high_energy, count_subthreshold_roots
export count_bound_states, continue_bound_state_counts
export levinson_phase_gate, mott_phase_gate
export bose_support_gate, convergence_gate, four_density_algorithm_labels
export joint_convergence_gate

const STRICT_SINGLE_CHARGE_OMEGA_MEASURE = :single_charge_domega_over_pi
const LEGACY_POSITIVE_ENERGY_OMEGA_MEASURE = :legacy_domega_over_2pi

"""Return the canonical positive-energy BU measure label."""
@inline function bu_omega_measure(measure::Symbol)::Symbol
    if measure === STRICT_SINGLE_CHARGE_OMEGA_MEASURE ||
       measure === :domega_over_pi || measure === :strict
        return STRICT_SINGLE_CHARGE_OMEGA_MEASURE
    elseif measure === LEGACY_POSITIVE_ENERGY_OMEGA_MEASURE ||
           measure === :domega_over_2pi || measure === :legacy
        return LEGACY_POSITIVE_ENERGY_OMEGA_MEASURE
    end
    throw(ArgumentError(
        "unknown BU omega measure $(measure); use :single_charge_domega_over_pi or :legacy_domega_over_2pi",
    ))
end

"""Return the multiplicative factor for the selected `domega` measure."""
@inline function bu_omega_measure_factor(measure::Symbol)::Float64
    canonical = bu_omega_measure(measure)
    return canonical === STRICT_SINGLE_CHARGE_OMEGA_MEASURE ? inv(pi) : inv(2pi)
end

function _validated_phase_profile(
    omega::AbstractVector{<:Real},
    phase::AbstractVector{<:Real},
)
    length(omega) == length(phase) ||
        throw(ArgumentError("omega and phase must have the same length"))
    length(omega) >= 2 || throw(ArgumentError("phase profile must contain at least two points"))
    omega_values = Float64.(omega)
    phase_values = Float64.(phase)
    all(isfinite, omega_values) || throw(ArgumentError("omega values must be finite"))
    all(isfinite, phase_values) || throw(ArgumentError("phase values must be finite"))
    all(diff(omega_values) .> 0.0) || throw(ArgumentError("omega values must be strictly increasing"))
    return omega_values, phase_values
end

function _unwrap_from_high_energy(phase::Vector{Float64}, branch_tol::Float64)
    branch_tol >= 0.0 || throw(ArgumentError("branch_tol must be nonnegative"))
    branch_tol < pi || throw(ArgumentError("branch_tol must be smaller than pi"))
    reversed = reverse(phase)
    unwrapped_reversed = similar(reversed)
    unwrapped_reversed[1] = reversed[1]
    shift = 0.0
    for i in 2:length(reversed)
        delta = reversed[i] - reversed[i - 1]
        if delta > pi - branch_tol
            shift -= 2pi
        elseif delta < -pi + branch_tol
            shift += 2pi
        end
        unwrapped_reversed[i] = reversed[i] + shift
    end
    return reverse(unwrapped_reversed)
end

"""
    anchor_phase_high_energy(omega, phase; target=0, branch_tol=0, tail_points=4)

Unwrap a real-axis phase from the largest energy toward lower energy and shift
the continuous branch so that its high-energy endpoint equals `target`.
Returns the raw, unwrapped, and anchored profiles plus tail diagnostics.
"""
function anchor_phase_high_energy(
    omega::AbstractVector{<:Real},
    phase::AbstractVector{<:Real};
    target::Real=0.0,
    branch_tol::Real=0.0,
    tail_points::Integer=4,
)
    omega_values, phase_values = _validated_phase_profile(omega, phase)
    target_value = Float64(target)
    isfinite(target_value) || throw(ArgumentError("target must be finite"))
    tail_points >= 2 || throw(ArgumentError("tail_points must be at least 2"))
    unwrapped = _unwrap_from_high_energy(phase_values, Float64(branch_tol))
    applied_shift = target_value - unwrapped[end]
    anchored = unwrapped .+ applied_shift
    first_tail = max(1, length(anchored) - Int(tail_points) + 1)
    tail = @view anchored[first_tail:end]
    adjacent_jumps = abs.(diff(anchored))
    return (
        omega=omega_values,
        raw_phase=phase_values,
        unwrapped_phase=unwrapped,
        anchored_phase=anchored,
        target=target_value,
        high_energy_phase_before_anchor=unwrapped[end],
        high_energy_phase_after_anchor=anchored[end],
        applied_shift=applied_shift,
        tail_points=length(tail),
        tail_span=maximum(tail) - minimum(tail),
        max_adjacent_jump=isempty(adjacent_jumps) ? 0.0 : maximum(adjacent_jumps),
    )
end

"""
    count_subthreshold_roots(omega, inverse_values, threshold; kwargs...)

Count simple real-axis roots below `threshold` by sign changes of the real part
of an inverse propagator. A non-negligible subthreshold imaginary part makes
the result diagnostic-only and returns `status=:complex_subthreshold`.
"""
function count_subthreshold_roots(
    omega::AbstractVector{<:Real},
    inverse_values::AbstractVector{<:Number},
    threshold::Real;
    zero_tolerance::Real=1.0e-10,
    imag_tolerance::Real=1.0e-8,
)
    length(omega) == length(inverse_values) ||
        throw(ArgumentError("omega and inverse_values must have the same length"))
    threshold_value = Float64(threshold)
    isfinite(threshold_value) || throw(ArgumentError("threshold must be finite"))
    zero_tolerance >= 0.0 || throw(ArgumentError("zero_tolerance must be nonnegative"))
    imag_tolerance >= 0.0 || throw(ArgumentError("imag_tolerance must be nonnegative"))
    omega_values = Float64.(omega)
    all(isfinite, omega_values) || throw(ArgumentError("omega values must be finite"))
    all(diff(omega_values) .> 0.0) || throw(ArgumentError("omega values must be strictly increasing"))
    all(z -> isfinite(real(z)) && isfinite(imag(z)), inverse_values) ||
        throw(ArgumentError("inverse_values must be finite"))

    indices = findall(<(threshold_value), omega_values)
    length(indices) >= 2 || return (
        count=0,
        brackets=Tuple{Float64,Float64}[],
        max_abs_imag=NaN,
        threshold=threshold_value,
        status=:insufficient_subthreshold_grid,
        passed=false,
    )
    re_values = Float64[real(inverse_values[i]) for i in indices]
    im_values = Float64[imag(inverse_values[i]) for i in indices]
    sub_omega = omega_values[indices]
    max_abs_imag = maximum(abs, im_values)

    brackets = Tuple{Float64,Float64}[]
    previous_index = nothing
    previous_sign = 0
    for i in eachindex(re_values)
        abs(re_values[i]) <= zero_tolerance && continue
        current_sign = signbit(re_values[i]) ? -1 : 1
        if previous_index !== nothing && current_sign != previous_sign
            push!(brackets, (sub_omega[previous_index], sub_omega[i]))
        end
        previous_index = i
        previous_sign = current_sign
    end
    status = max_abs_imag <= imag_tolerance ? :ok : :complex_subthreshold
    return (
        count=length(brackets),
        brackets=brackets,
        max_abs_imag=max_abs_imag,
        threshold=threshold_value,
        status=status,
        passed=status === :ok,
    )
end

"""
    count_bound_states(inverse_fn, q, threshold; kwargs...)

Independently sample a retarded inverse propagator below the two-particle
threshold and count simple real-axis zero brackets.  A finite imaginary part
does not get silently discarded: the result is marked `:complex_subthreshold`
and is diagnostic-only.  This routine is intentionally separate from phase
unwrapping/Levinson evaluation.
"""
function count_bound_states(
    inverse_fn,
    q::Real,
    threshold::Real;
    omega_min::Real=0.0,
    omega_nodes::Integer=64,
    zero_tolerance::Real=1.0e-10,
    imag_tolerance::Real=1.0e-8,
)
    q_value = Float64(q)
    threshold_value = Float64(threshold)
    lower = Float64(omega_min)
    all(isfinite, (q_value, threshold_value, lower)) ||
        throw(ArgumentError("q, threshold, and omega_min must be finite"))
    q_value >= 0.0 || throw(ArgumentError("q must be nonnegative"))
    threshold_value > lower || throw(ArgumentError("threshold must exceed omega_min"))
    omega_nodes >= 4 || throw(ArgumentError("omega_nodes must be at least 4"))
    omega = collect(range(lower, threshold_value; length=Int(omega_nodes) + 1))[1:end-1]
    values = ComplexF64[inverse_fn(Float64(w), q_value) for w in omega]
    roots = count_subthreshold_roots(
        omega,
        values,
        threshold_value;
        zero_tolerance=zero_tolerance,
        imag_tolerance=imag_tolerance,
    )
    return merge(roots, (q=q_value, omega_nodes=Int(omega_nodes),
                         omega_min=lower, independent=true,
                         counting_method=:subthreshold_sign_brackets))
end

"""
    continue_bound_state_counts(inverse_fn, q_values, threshold_fn; kwargs...)

Run the independent bound-state counter along a momentum continuation.  No
count is inferred from the previous point; `continuation_delta` is reported so
that a branch jump can be reviewed explicitly.
"""
function continue_bound_state_counts(
    inverse_fn,
    q_values::AbstractVector{<:Real},
    threshold_fn;
    kwargs...
)
    isempty(q_values) && throw(ArgumentError("q_values must not be empty"))
    previous_count = nothing
    rows = NamedTuple[]
    for q in q_values
        q_value = Float64(q)
        threshold = threshold_fn isa Function ? threshold_fn(q_value) : threshold_fn
        result = count_bound_states(inverse_fn, q_value, Float64(threshold); kwargs...)
        delta = previous_count === nothing ? missing : result.count - previous_count
        push!(rows, merge(result, (
            continuation_previous_count=previous_count === nothing ? missing : previous_count,
            continuation_delta=delta,
            continuation_index=length(rows) + 1,
        )))
        previous_count = result.count
    end
    return rows
end

function _linear_interpolate(x::Vector{Float64}, y::Vector{Float64}, point::Float64)
    x[1] <= point <= x[end] || throw(ArgumentError("interpolation point is outside the profile"))
    index = searchsortedlast(x, point)
    index == length(x) && return y[end]
    x[index] == point && return y[index]
    fraction = (point - x[index]) / (x[index + 1] - x[index])
    return muladd(fraction, y[index + 1] - y[index], y[index])
end

"""
    levinson_phase_gate(omega, phase, threshold; bound_state_count, ...)

Anchor the phase at high energy and test
`delta(threshold) - delta(infinity) = pi * bound_state_count`. The gate also
requires a stable high-energy tail before accepting the profile.
"""
function levinson_phase_gate(
    omega::AbstractVector{<:Real},
    phase::AbstractVector{<:Real},
    threshold::Real;
    bound_state_count::Integer,
    target::Real=0.0,
    branch_tol::Real=0.0,
    tail_points::Integer=4,
    phase_tolerance::Real=0.05pi,
    tail_tolerance::Real=0.02pi,
)
    bound_state_count >= 0 || throw(ArgumentError("bound_state_count must be nonnegative"))
    phase_tolerance >= 0.0 || throw(ArgumentError("phase_tolerance must be nonnegative"))
    tail_tolerance >= 0.0 || throw(ArgumentError("tail_tolerance must be nonnegative"))
    profile = anchor_phase_high_energy(
        omega,
        phase;
        target=target,
        branch_tol=branch_tol,
        tail_points=tail_points,
    )
    threshold_value = Float64(threshold)
    threshold_phase = _linear_interpolate(profile.omega, profile.anchored_phase, threshold_value)
    expected_phase = profile.target + pi * Int(bound_state_count)
    residual = threshold_phase - expected_phase
    phase_passed = abs(residual) <= phase_tolerance
    tail_passed = profile.tail_span <= tail_tolerance
    return (
        passed=phase_passed && tail_passed,
        bound_state_count=Int(bound_state_count),
        threshold=threshold_value,
        threshold_phase=threshold_phase,
        high_energy_target=profile.target,
        expected_threshold_phase=expected_phase,
        levinson_residual=residual,
        phase_tolerance=Float64(phase_tolerance),
        tail_span=profile.tail_span,
        tail_tolerance=Float64(tail_tolerance),
        phase_passed=phase_passed,
        tail_passed=tail_passed,
        profile=profile,
    )
end

"""
    mott_phase_gate(before, after; expected_bound_state_drop=1, ...)

Require a Levinson-consistent loss of bound states and the matching decrease
of the threshold phase across a proposed Mott transition.
"""
function mott_phase_gate(
    before::NamedTuple,
    after::NamedTuple;
    expected_bound_state_drop::Integer=1,
    phase_tolerance::Real=0.10pi,
)
    expected_bound_state_drop >= 0 ||
        throw(ArgumentError("expected_bound_state_drop must be nonnegative"))
    phase_tolerance >= 0.0 || throw(ArgumentError("phase_tolerance must be nonnegative"))
    required = (:passed, :bound_state_count, :threshold_phase)
    all(name -> hasproperty(before, name) && hasproperty(after, name), required) ||
        throw(ArgumentError("before and after must be levinson_phase_gate results"))
    count_drop = Int(before.bound_state_count) - Int(after.bound_state_count)
    phase_drop = Float64(before.threshold_phase) - Float64(after.threshold_phase)
    expected_phase_drop = pi * Int(expected_bound_state_drop)
    phase_residual = phase_drop - expected_phase_drop
    count_passed = count_drop == expected_bound_state_drop
    phase_passed = abs(phase_residual) <= phase_tolerance
    endpoints_passed = Bool(before.passed) && Bool(after.passed)
    return (
        passed=endpoints_passed && count_passed && phase_passed,
        before_passed=Bool(before.passed),
        after_passed=Bool(after.passed),
        bound_state_count_before=Int(before.bound_state_count),
        bound_state_count_after=Int(after.bound_state_count),
        bound_state_count_drop=count_drop,
        expected_bound_state_drop=Int(expected_bound_state_drop),
        threshold_phase_drop=phase_drop,
        expected_threshold_phase_drop=expected_phase_drop,
        phase_residual=phase_residual,
        phase_tolerance=Float64(phase_tolerance),
        count_passed=count_passed,
        phase_passed=phase_passed,
    )
end

"""
    bose_support_gate(mass, chemical_potential; omega_min, omega_max)

Check the normal-phase Bose support needed by a positive-energy integral. The
minimum physical excitation is `mass` at `q=0`; both it and the requested
integration lower bound must lie above the chemical potential.
"""
function bose_support_gate(
    mass::Real,
    chemical_potential::Real;
    omega_min::Real,
    omega_max::Real,
)
    m = Float64(mass)
    μ = Float64(chemical_potential)
    lo = Float64(omega_min)
    hi = Float64(omega_max)
    all(isfinite, (m, μ, lo, hi)) || throw(ArgumentError("mass, chemical potential, and omega bounds must be finite"))
    m >= 0.0 || throw(ArgumentError("mass must be nonnegative"))
    hi > lo || throw(ArgumentError("omega_max must exceed omega_min"))
    safe_excitation = m > μ
    safe_window = lo > μ
    passed = safe_excitation && safe_window
    status = passed ? :safe_normal_domain : :unsafe_bose_domain
    return (
        passed=passed,
        status=status,
        mass=m,
        chemical_potential=μ,
        omega_min=lo,
        omega_max=hi,
        min_E_minus_mu=m - μ,
        support_lower_bound=max(lo, μ),
        excitation_safe=safe_excitation,
        integration_window_safe=safe_window,
    )
end

"""Compare one numerical value with a reference under an explicit tolerance."""
function convergence_gate(
    reference::Real,
    candidate::Real;
    rtol::Real=1.0e-3,
    atol::Real=1.0e-10,
)
    ref = Float64(reference)
    cand = Float64(candidate)
    rel_tol = Float64(rtol)
    abs_tol = Float64(atol)
    all(isfinite, (ref, cand, rel_tol, abs_tol)) || throw(ArgumentError("convergence inputs must be finite"))
    rel_tol >= 0.0 || throw(ArgumentError("rtol must be nonnegative"))
    abs_tol >= 0.0 || throw(ArgumentError("atol must be nonnegative"))
    difference = cand - ref
    scale = max(abs(ref), abs(cand), abs_tol)
    relative_difference = abs(difference) / scale
    passed = abs(difference) <= abs_tol + rel_tol * abs(ref)
    return (
        passed=passed,
        reference=ref,
        candidate=cand,
        absolute_difference=abs(difference),
        relative_difference=relative_difference,
        rtol=rel_tol,
        atol=abs_tol,
    )
end

"""
    joint_convergence_gate(samples; value_field=:density, rtol=..., atol=...)

Apply pairwise numerical convergence checks to an ordered sequence of strict
diagnostic results.  The contract also requires explicit acceptance and,
when present, a stable endpoint flag on every sample.  Optional metadata such
as `eta_inv_fm`, `q_nodes`, `omega_nodes`, `qmax_inv_fm`, and `omega_max_inv_fm`
is copied into the returned audit record; no axis is silently treated as
converged merely because a density is finite.
"""
function joint_convergence_gate(
    samples::AbstractVector;
    value_field::Symbol=:density,
    rtol::Real=0.05,
    atol::Real=1.0e-10,
    require_accepted::Bool=true,
)
    length(samples) >= 2 || throw(ArgumentError("at least two convergence samples are required"))
    all(sample -> hasproperty(sample, value_field), samples) ||
        throw(ArgumentError("every sample must provide $(value_field)"))
    finite = all(sample -> isfinite(Float64(getproperty(sample, value_field))), samples)
    pairwise = [
        convergence_gate(
            Float64(getproperty(samples[i], value_field)),
            Float64(getproperty(samples[i + 1], value_field));
            rtol=rtol,
            atol=atol,
        ) for i in 1:(length(samples) - 1)
    ]
    accepted = all(sample -> !require_accepted ||
        (hasproperty(sample, :accepted) && Bool(getproperty(sample, :accepted))), samples)
    endpoint_stable = all(sample -> !hasproperty(sample, :tail_stable) ||
        Bool(getproperty(sample, :tail_stable)), samples)
    return (
        passed=finite && accepted && endpoint_stable && all(item -> item.passed, pairwise),
        finite=finite,
        accepted=accepted,
        endpoint_stable=endpoint_stable,
        pairwise=pairwise,
        sample_count=length(samples),
        value_field=value_field,
        rtol=Float64(rtol),
        atol=Float64(atol),
        metadata=map(samples) do sample
            (
                eta_inv_fm=hasproperty(sample, :eta_inv_fm) ? getproperty(sample, :eta_inv_fm) : missing,
                qmax_inv_fm=hasproperty(sample, :qmax_inv_fm) ? getproperty(sample, :qmax_inv_fm) : missing,
                q_nodes=hasproperty(sample, :q_nodes) ? getproperty(sample, :q_nodes) : missing,
                omega_max_inv_fm=hasproperty(sample, :omega_max_inv_fm) ? getproperty(sample, :omega_max_inv_fm) : missing,
                omega_nodes=hasproperty(sample, :omega_nodes) ? getproperty(sample, :omega_nodes) : missing,
            )
        end,
    )
end

"""Canonical four density algorithms used by the comparison and scan layer."""
@inline four_density_algorithm_labels() = (
    :stable_particle_limit,
    :reduced_strict_bw,
    :q_pole_strict_bw,
    :phase_shift_bu,
)

end # module BUPhaseGates
