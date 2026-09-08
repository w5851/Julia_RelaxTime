"""
    PhaseNormalization

Pure algebraic contracts connecting a scattering phase, an on-shell
S-matrix, and the Beth-Uhlenbeck energy measure.  This module deliberately
does not infer a physical phase from a propagator: that mapping is model and
analytic-continuation dependent and remains owned by the charged phase
backend.
"""
module PhaseNormalization

using LinearAlgebra: tr

export SCATTERING_PHASE, SMATRIX_ARGUMENT
export phase_variable
export phase_to_s_matrix, s_matrix_argument, s_matrix_to_phase
export phase_measure_factor, phase_measure
export s_matrix_log_derivative, s_matrix_density_of_states
export phase_to_s_matrix_profile, s_matrix_to_phase_profile

const SCATTERING_PHASE = :delta
const SMATRIX_ARGUMENT = :s_matrix_argument

@inline function phase_variable(variable::Symbol)::Symbol
    variable === SCATTERING_PHASE && return SCATTERING_PHASE
    variable === SMATRIX_ARGUMENT && return SMATRIX_ARGUMENT
    variable === :arg_s && return SMATRIX_ARGUMENT
    throw(ArgumentError(
        "phase variable must be :delta or :s_matrix_argument",
    ))
end

@inline function _finite_real(value, label::AbstractString)::Float64
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("$(label) must be finite"))
    return result
end

@inline function _finite_complex(value, label::AbstractString)::ComplexF64
    value isa Number || throw(ArgumentError("$(label) must be numeric"))
    result = ComplexF64(value)
    isfinite(real(result)) && isfinite(imag(result)) ||
        throw(ArgumentError("$(label) must be finite"))
    return result
end

"""
    phase_to_s_matrix(phase; variable=:delta)

Map the physical scattering phase ``delta`` to the scalar elastic
S-matrix.  If `variable=:s_matrix_argument`, the input is already
`arg(S)` and is exponentiated without the factor two.
"""
@inline function phase_to_s_matrix(phase::Real; variable::Symbol=SCATTERING_PHASE)
    canonical = phase_variable(variable)
    value = _finite_real(phase, "phase")
    return cis(canonical === SCATTERING_PHASE ? 2.0 * value : value)
end

"""Return the principal scalar argument of a finite S-matrix value."""
@inline function s_matrix_argument(value::Number)::Float64
    z = _finite_complex(value, "S-matrix value")
    abs(z) > 0.0 || throw(ArgumentError("S-matrix value must be nonzero"))
    return atan(imag(z), real(z))
end

"""
    s_matrix_to_phase(value; variable=:delta)

Return the principal phase represented by a scalar S-matrix value.  The
result is a local principal value; continuous branch reconstruction belongs
to `BUPhaseGates.anchor_phase_high_energy`.
"""
@inline function s_matrix_to_phase(value::Number; variable::Symbol=SCATTERING_PHASE)
    canonical = phase_variable(variable)
    argument = s_matrix_argument(value)
    return canonical === SCATTERING_PHASE ? argument / 2.0 : argument
end

"""Return the BU derivative measure for the selected phase variable."""
@inline function phase_measure_factor(variable::Symbol=SCATTERING_PHASE)::Float64
    canonical = phase_variable(variable)
    return canonical === SCATTERING_PHASE ? inv(pi) : inv(2pi)
end

@inline function phase_measure(variable::Symbol=SCATTERING_PHASE)::Symbol
    canonical = phase_variable(variable)
    return canonical === SCATTERING_PHASE ? :d_delta_over_pi : :d_arg_s_over_2pi
end

"""
    s_matrix_log_derivative(S, dS)

Return the imaginary logarithmic derivative `Im(S⁻¹ dS)`.  For a unitary
scalar `S=exp(2im*delta)`, this is exactly `2*d_delta`.
"""
@inline function s_matrix_log_derivative(S::Number, dS::Number)::Float64
    s = _finite_complex(S, "S-matrix value")
    ds = _finite_complex(dS, "S-matrix derivative")
    abs(s) > 0.0 || throw(ArgumentError("S-matrix value must be nonzero"))
    return imag(ds / s)
end

"""
    s_matrix_density_of_states(S, dS)

Return the scalar DMB/Beth-Uhlenbeck density-of-states contribution
`Im(S⁻¹ dS)/(2π)`.  For `S=exp(2im*delta)` this equals `d_delta/π`.
"""
@inline function s_matrix_density_of_states(S::Number, dS::Number)::Float64
    return phase_measure_factor(SMATRIX_ARGUMENT) * s_matrix_log_derivative(S, dS)
end

"""Return a diagonal synthetic S-matrix from a phase profile."""
function phase_to_s_matrix_profile(
    phase::AbstractVector{<:Real};
    variable::Symbol=SCATTERING_PHASE,
)
    isempty(phase) && throw(ArgumentError("phase profile must not be empty"))
    return ComplexF64[phase_to_s_matrix(value; variable=variable) for value in phase]
end

"""Convert a scalar S-matrix profile to principal phase values."""
function s_matrix_to_phase_profile(
    values::AbstractVector{<:Number};
    variable::Symbol=SCATTERING_PHASE,
)
    isempty(values) && throw(ArgumentError("S-matrix profile must not be empty"))
    return Float64[s_matrix_to_phase(value; variable=variable) for value in values]
end

"""
    s_matrix_log_derivative(S, dS)

For a diagonal matrix S and its energy derivative dS, return the connected
trace `Im tr(S⁻¹ dS)`.  A matrix route must use the argument-of-determinant
measure `d arg(det S)/(2π)`; individual eigenphases may still be tracked
separately when the channel basis is physically fixed.
"""
function s_matrix_log_derivative(S::AbstractMatrix, dS::AbstractMatrix)::Float64
    size(S) == size(dS) || throw(ArgumentError("S and dS must have the same size"))
    size(S, 1) == size(S, 2) || throw(ArgumentError("S must be square"))
    all(isfinite, real.(S)) && all(isfinite, imag.(S)) ||
        throw(ArgumentError("S must be finite"))
    all(isfinite, real.(dS)) && all(isfinite, imag.(dS)) ||
        throw(ArgumentError("dS must be finite"))
    result = imag(tr(S \ dS))
    isfinite(result) || throw(ArgumentError("S-matrix logarithmic derivative must be finite"))
    return Float64(result)
end

@inline function s_matrix_density_of_states(S::AbstractMatrix, dS::AbstractMatrix)::Float64
    return phase_measure_factor(SMATRIX_ARGUMENT) * s_matrix_log_derivative(S, dS)
end

end # module PhaseNormalization
