raw"""
    ChargedPhaseBackend

Strict, solver-independent real-axis phase backend for charged RPA/BU.

The backend consumes an ordered retarded inverse propagator ``Δ^R(ω,q)`` and
implements the project formula contract

```math
δ(ω,q) = -\operatorname{arg} Δ^R(ω,q),
\qquad
n_M = \frac{d_M}{T}\int\frac{dq\,q^2}{2π^2}
       \int\frac{dω}{π}\,g_B(ω)\frac{∂δ}{∂ω}.
```

It is intentionally parallel to the legacy `MesonDensity` phase path.  It
does not change any production default, solve a PNJL equilibrium state, or
perform second-sheet pole continuation.
"""
module ChargedPhaseBackend

using ..GaussLegendre: gauleg
using ..ChargedRPAKernel: ChargedRPAKernelSpec, charged_rpa_inverse
using ..BUPhaseGates: STRICT_SINGLE_CHARGE_OMEGA_MEASURE,
                      bu_omega_measure,
                      bu_omega_measure_factor,
                      bose_support_gate,
                      count_subthreshold_roots,
                      levinson_phase_gate,
                      mott_phase_gate,
                      anchor_phase_high_energy,
                      convergence_gate

export StrictChargedPhaseSpec
export strict_retarded_phase
export strict_phase_profile, strict_phase_gate
export strict_charged_bu_density
export strict_charged_rpa_bu_density
export strict_mott_gate
export strict_density_convergence_gate

const _VALID_PHASE_OBJECTS = (:inverse_propagator, :propagator)

"""Numerical choices that define the strict phase convention."""
struct StrictChargedPhaseSpec
    phase_object::Symbol
    phase_sign::Int
    target::Float64
    branch_tol::Float64
    tail_points::Int
    tail_tolerance::Float64
end

function StrictChargedPhaseSpec(;
    phase_object::Symbol=:inverse_propagator,
    phase_sign::Integer=(phase_object === :inverse_propagator ? -1 : 1),
    target::Real=0.0,
    branch_tol::Real=0.0,
    tail_points::Integer=4,
    tail_tolerance::Real=0.02 * π,
)
    phase_object in _VALID_PHASE_OBJECTS || throw(ArgumentError(
        "phase_object must be :inverse_propagator or :propagator",
    ))
    phase_sign in (-1, 1) || throw(ArgumentError("phase_sign must be +1 or -1"))
    values = (Float64(target), Float64(branch_tol), Float64(tail_tolerance))
    all(isfinite, values) || throw(ArgumentError("phase specification values must be finite"))
    branch_tol >= 0.0 || throw(ArgumentError("branch_tol must be nonnegative"))
    tail_tolerance >= 0.0 || throw(ArgumentError("tail_tolerance must be nonnegative"))
    tail_points >= 2 || throw(ArgumentError("tail_points must be at least 2"))
    return StrictChargedPhaseSpec(
        phase_object,
        Int(phase_sign),
        values[1],
        values[2],
        Int(tail_points),
        values[3],
    )
end

@inline function _finite_complex(value, label::AbstractString)::ComplexF64
    value isa Number || throw(ArgumentError("$(label) must be numeric"))
    z = ComplexF64(value)
    isfinite(real(z)) && isfinite(imag(z)) || throw(ArgumentError("$(label) must be finite"))
    return z
end

@inline function _validate_grid(omega::AbstractVector{<:Real})
    values = Float64.(omega)
    length(values) >= 2 || throw(ArgumentError("omega grid must contain at least two points"))
    all(isfinite, values) || throw(ArgumentError("omega grid must be finite"))
    all(diff(values) .> 0.0) || throw(ArgumentError("omega grid must be strictly increasing"))
    return values
end

"""Evaluate the canonical principal-value retarded phase at one point."""
@inline function strict_retarded_phase(
    inverse_value::Number;
    phase_object::Symbol=:inverse_propagator,
    phase_sign::Integer=(phase_object === :inverse_propagator ? -1 : 1),
)
    phase_object in _VALID_PHASE_OBJECTS || throw(ArgumentError(
        "phase_object must be :inverse_propagator or :propagator",
    ))
    phase_sign in (-1, 1) || throw(ArgumentError("phase_sign must be +1 or -1"))
    z = _finite_complex(inverse_value, "inverse_value")
    return Float64(phase_sign) * atan(imag(z), real(z))
end

"""
    strict_phase_profile(omega, inverse_values; spec=StrictChargedPhaseSpec())

Construct the raw, high-energy-unwrapped and anchored phase profile.  The
anchor is not silently accepted: `tail_stable` is false when the finite high
energy tail is still changing by more than `spec.tail_tolerance`.
"""
function strict_phase_profile(
    omega::AbstractVector{<:Real},
    inverse_values::AbstractVector{<:Number};
    spec::StrictChargedPhaseSpec=StrictChargedPhaseSpec(),
)
    ω = _validate_grid(omega)
    length(ω) == length(inverse_values) || throw(ArgumentError(
        "omega and inverse_values must have the same length",
    ))
    values = ComplexF64[_finite_complex(value, "inverse_values[$i]") for (i, value) in enumerate(inverse_values)]
    raw = Float64[
        strict_retarded_phase(value; phase_object=spec.phase_object, phase_sign=spec.phase_sign)
        for value in values
    ]
    anchored = anchor_phase_high_energy(
        ω,
        raw;
        target=spec.target,
        branch_tol=spec.branch_tol,
        tail_points=spec.tail_points,
    )
    return (
        omega=ω,
        inverse_values=values,
        raw_phase=raw,
        unwrapped_phase=anchored.unwrapped_phase,
        anchored_phase=anchored.anchored_phase,
        phase_object=spec.phase_object,
        phase_sign=spec.phase_sign,
        target=spec.target,
        branch_tol=spec.branch_tol,
        applied_shift=anchored.applied_shift,
        high_energy_phase_before_anchor=anchored.high_energy_phase_before_anchor,
        high_energy_phase_after_anchor=anchored.high_energy_phase_after_anchor,
        tail_span=anchored.tail_span,
        tail_points=anchored.tail_points,
        max_adjacent_jump=anchored.max_adjacent_jump,
        tail_stable=anchored.tail_span <= spec.tail_tolerance,
        tail_tolerance=spec.tail_tolerance,
    )
end

"""Apply the endpoint and Levinson checks to one strict phase profile."""
function strict_phase_gate(
    profile;
    threshold::Real,
    bound_state_count::Integer,
    phase_tolerance::Real=0.05 * π,
    imag_tolerance::Real=1.0e-8,
)
    isfinite(Float64(threshold)) || throw(ArgumentError("threshold must be finite"))
    roots = count_subthreshold_roots(
        profile.omega,
        profile.inverse_values,
        threshold;
        imag_tolerance=imag_tolerance,
    )
    levinson = levinson_phase_gate(
        profile.omega,
        profile.raw_phase,
        threshold;
        bound_state_count=bound_state_count,
        target=profile.target,
        branch_tol=profile.branch_tol,
        phase_tolerance=phase_tolerance,
        tail_tolerance=profile.tail_tolerance,
        tail_points=profile.tail_points,
    )
    passed = profile.tail_stable && roots.passed && levinson.passed
    return (
        passed=passed,
        tail_stable=profile.tail_stable,
        roots=roots,
        levinson=levinson,
        threshold=Float64(threshold),
        bound_state_count=Int(bound_state_count),
    )
end

"""Run the same endpoint/Levinson contract on two profiles and then apply
the project Mott-transition gate.

The two profiles must use the same phase convention.  This helper does not
infer a bound-state count or a Mott temperature; both are explicit inputs so
that an unsuccessful physical identification remains visible in diagnostics.
"""
function strict_mott_gate(
    before_profile,
    after_profile;
    before_threshold::Real,
    after_threshold::Real,
    before_bound_state_count::Integer,
    after_bound_state_count::Integer,
    phase_tolerance::Real=0.05 * π,
    expected_bound_state_drop::Integer=1,
    transition_phase_tolerance::Real=0.10 * π,
)
    before = strict_phase_gate(
        before_profile;
        threshold=before_threshold,
        bound_state_count=before_bound_state_count,
        phase_tolerance=phase_tolerance,
    )
    after = strict_phase_gate(
        after_profile;
        threshold=after_threshold,
        bound_state_count=after_bound_state_count,
        phase_tolerance=phase_tolerance,
    )
    transition = mott_phase_gate(
        before.levinson,
        after.levinson;
        expected_bound_state_drop=expected_bound_state_drop,
        phase_tolerance=transition_phase_tolerance,
    )
    return (
        passed=before.passed && after.passed && transition.passed,
        before=before,
        after=after,
        transition=transition,
    )
end

@inline function _bose(E::Float64, μ::Float64, T::Float64)
    E > μ || throw(ArgumentError("strict BU requires omega > chemical potential"))
    x = (E - μ) / T
    x > 700.0 && return 0.0
    return inv(expm1(x))
end

function _trapz(x::Vector{Float64}, y::Vector{Float64})
    total = 0.0
    @inbounds for i in 1:(length(x) - 1)
        total += 0.5 * (x[i + 1] - x[i]) * (y[i + 1] + y[i])
    end
    return total
end

function _phase_derivative(omega::Vector{Float64}, phase::Vector{Float64})
    out = similar(phase)
    out[1] = (phase[2] - phase[1]) / (omega[2] - omega[1])
    out[end] = (phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    @inbounds for i in 2:(length(phase) - 1)
        out[i] = (phase[i + 1] - phase[i - 1]) / (omega[i + 1] - omega[i - 1])
    end
    return out
end

"""
    strict_charged_bu_density(inverse_fn, mass, T; kwargs...) -> NamedTuple

Integrate a charge-resolved strict BU density from an ordered retarded inverse
propagator callback `inverse_fn(omega, q)`.  `threshold` and
`bound_state_count` are required for an accepted Levinson gate.  The function
always returns diagnostics; when a gate fails, `accepted=false` and
`status=:gate_failed` prevent the value from being treated as production.
"""
function strict_charged_bu_density(
    inverse_fn,
    mass::Real,
    T::Real;
    μ::Real=0.0,
    degeneracy::Integer=1,
    qmax::Real=12.0,
    q_nodes::Integer=24,
    omega_min::Real=0.05,
    omega_max::Real=10.0,
    omega_nodes::Integer=96,
    threshold::Union{Nothing,Real,Function}=nothing,
    bound_state_count::Union{Nothing,Integer,Function}=nothing,
    phase_spec::StrictChargedPhaseSpec=StrictChargedPhaseSpec(),
    omega_measure::Symbol=STRICT_SINGLE_CHARGE_OMEGA_MEASURE,
    phase_tolerance::Real=0.05 * π,
    require_levinson::Bool=true,
)
    m, temperature, chemical_potential = Float64.((mass, T, μ))
    all(isfinite, (m, temperature, chemical_potential)) ||
        throw(ArgumentError("mass, T, and mu must be finite"))
    m >= 0.0 || throw(ArgumentError("mass must be nonnegative"))
    temperature > 0.0 || throw(ArgumentError("T must be positive"))
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive"))
    qmax_value, ωmin, ωmax = Float64.((qmax, omega_min, omega_max))
    qmax_value > 0.0 || throw(ArgumentError("qmax must be positive"))
    ωmax > ωmin || throw(ArgumentError("omega_max must exceed omega_min"))
    q_nodes >= 2 || throw(ArgumentError("q_nodes must be at least 2"))
    omega_nodes >= 4 || throw(ArgumentError("omega_nodes must be at least 4"))
    if require_levinson
        threshold === nothing && throw(ArgumentError("threshold is required when require_levinson=true"))
        bound_state_count === nothing && throw(ArgumentError("bound_state_count is required when require_levinson=true"))
    end
    threshold isa Real && !isfinite(Float64(threshold)) && throw(ArgumentError("threshold must be finite"))
    threshold isa Function || threshold === nothing || threshold isa Real ||
        throw(ArgumentError("threshold must be a real value, callable, or nothing"))
    bound_state_count isa Integer && bound_state_count < 0 &&
        throw(ArgumentError("bound_state_count must be nonnegative"))
    bound_state_count isa Function || bound_state_count === nothing || bound_state_count isa Integer ||
        throw(ArgumentError("bound_state_count must be an integer, callable, or nothing"))
    support = bose_support_gate(m, chemical_potential; omega_min=ωmin, omega_max=ωmax)
    measure = bu_omega_measure(omega_measure)
    factor = bu_omega_measure_factor(measure)
    if !support.passed
        return (
            density=NaN,
            accepted=false,
            status=:unsafe_bose_domain,
            message="strict charged BU normal-phase Bose support failed",
            omega_measure=measure,
            omega_measure_factor=factor,
            bose_support=support,
            qmax=qmax_value,
            q_nodes=Int(q_nodes),
            omega_min=ωmin,
            omega_max=ωmax,
            omega_nodes=Int(omega_nodes),
            q_profiles=NamedTuple[],
            failed_q_count=0,
            threshold=threshold isa Real ? Float64(threshold) : NaN,
            bound_state_count=bound_state_count isa Integer ? Int(bound_state_count) : -1,
            threshold_policy=threshold isa Function ? :q_callable : (threshold === nothing ? :none : :scalar),
            bound_state_policy=bound_state_count isa Function ? :q_callable :
                (bound_state_count === nothing ? :none : :scalar),
            density_finite=false,
            density_nonnegative=false,
        )
    end

    q_grid, q_weights = gauleg(0.0, qmax_value, Int(q_nodes))
    ω_grid, _ = gauleg(ωmin, ωmax, Int(omega_nodes))
    q_profiles = NamedTuple[]
    q_integral = 0.0
    failed_q_count = 0

    @inbounds for iq in eachindex(q_grid, q_weights)
        q = Float64(q_grid[iq])
        values = ComplexF64[_finite_complex(inverse_fn(Float64(ω), q), "inverse_fn value") for ω in ω_grid]
        profile = strict_phase_profile(ω_grid, values; spec=phase_spec)
        threshold_q = threshold isa Function ? Float64(threshold(q)) :
            (threshold === nothing ? NaN : Float64(threshold))
        bound_state_count_q = bound_state_count isa Function ? Int(bound_state_count(q)) :
            (bound_state_count === nothing ? 0 : Int(bound_state_count))
        isfinite(threshold_q) && ωmin < threshold_q < ωmax || threshold === nothing ||
            throw(ArgumentError("threshold(q) must be finite and inside the omega grid"))
        bound_state_count_q >= 0 || throw(ArgumentError("bound_state_count(q) must be nonnegative"))
        gate = threshold === nothing ? nothing : strict_phase_gate(
            profile;
            threshold=threshold_q,
            bound_state_count=bound_state_count_q,
            phase_tolerance=Float64(phase_tolerance),
        )
        accepted_q = profile.tail_stable && (!require_levinson || Bool(gate.passed))
        accepted_q || (failed_q_count += 1)

        derivative = _phase_derivative(profile.omega, profile.anchored_phase)
        # The phase profile contains the positive bound-state jump below the
        # continuum threshold and the compensating continuum fall-off.  The
        # BU spectral measure is the signed derivative `d(delta)/domega`.
        integrand = Float64[
            _bose(Float64(ω), chemical_potential, temperature) * derivative[i]
            for (i, ω) in enumerate(profile.omega)
        ]
        shell = (q^2 / (2.0 * π^2)) * factor * _trapz(profile.omega, integrand)
        q_integral += Float64(q_weights[iq]) * shell
        push!(q_profiles, (
            q=q,
            threshold=threshold_q,
            bound_state_count=bound_state_count_q,
            shell=Float64(shell),
            accepted=accepted_q,
            profile=profile,
            gate=gate,
        ))
    end

    density = Float64(degeneracy) * q_integral / temperature
    density_finite = isfinite(density)
    density_nonnegative = density >= -1.0e-12
    accepted = failed_q_count == 0 && density_finite && density_nonnegative
    status = !density_finite ? :invalid_density :
        (!density_nonnegative ? :invalid_density : (accepted ? :ok : :gate_failed))
    return (
        density=density,
        accepted=accepted,
        status=status,
        message=accepted ? "" :
            (density_nonnegative ? "one or more q shells failed strict phase gates" :
             "strict charged BU density is negative or non-finite"),
        omega_measure=measure,
        omega_measure_factor=factor,
        bose_support=support,
        qmax=qmax_value,
        q_nodes=Int(q_nodes),
        omega_min=ωmin,
        omega_max=ωmax,
        omega_nodes=Int(omega_nodes),
        threshold=threshold isa Real ? Float64(threshold) : NaN,
        bound_state_count=bound_state_count isa Integer ? Int(bound_state_count) : -1,
        threshold_policy=threshold isa Function ? :q_callable : (threshold === nothing ? :none : :scalar),
        bound_state_policy=bound_state_count isa Function ? :q_callable :
            (bound_state_count === nothing ? :none : :scalar),
        density_finite=density_finite,
        density_nonnegative=density_nonnegative,
        failed_q_count=failed_q_count,
        q_integral=q_integral,
        q_profiles=q_profiles,
    )
end

"""
    strict_charged_rpa_bu_density(spec, coupling, polarization_fn, mass, T; kwargs...)

Compose the validated charged single-channel RPA inverse
`1 - 4*K*Pi^R_ordered` with the strict phase/BU integrator.  `polarization_fn`
must return the ordered retarded complex bubble at `(omega, q)`.  This is the
non-public diagnostic adapter between `ChargedRPAProvider` and the generic
phase backend; it does not alter any legacy density entrypoint.
"""
function strict_charged_rpa_bu_density(
    spec::ChargedRPAKernelSpec,
    coupling::Real,
    polarization_fn,
    mass::Real,
    T::Real;
    kwargs...
)
    inverse_fn = (omega, q) -> charged_rpa_inverse(
        spec,
        coupling,
        polarization_fn(omega, q),
    )
    return strict_charged_bu_density(inverse_fn, mass, T; kwargs...)
end

"""Compare coarse/refined strict density outputs without hiding gate failures."""
function strict_density_convergence_gate(
    coarse,
    refined;
    rtol::Real=0.05,
    atol::Real=1.0e-10,
)
    hasproperty(coarse, :density) && hasproperty(refined, :density) ||
        throw(ArgumentError("coarse and refined results must provide density"))
    finite = isfinite(Float64(coarse.density)) && isfinite(Float64(refined.density))
    numeric = finite ? convergence_gate(Float64(coarse.density), Float64(refined.density); rtol=rtol, atol=atol) : nothing
    accepted = finite && Bool(getproperty(coarse, :accepted)) && Bool(getproperty(refined, :accepted)) && Bool(numeric.passed)
    return (
        passed=accepted,
        finite=finite,
        numeric=numeric,
        coarse_accepted=Bool(getproperty(coarse, :accepted)),
        refined_accepted=Bool(getproperty(refined, :accepted)),
        rtol=Float64(rtol),
        atol=Float64(atol),
    )
end

end # module ChargedPhaseBackend
