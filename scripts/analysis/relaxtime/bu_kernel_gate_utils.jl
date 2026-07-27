"""
Pure numerical helpers for the finite-window Beth-Uhlenbeck identity gate.

This file belongs to the analysis layer. It deliberately does not import the
production meson-density implementation, so the integration-by-parts identity
can be unit-tested without solving a PNJL state.
"""
module BUKernelGateUtils

export bose_occupation
export finite_window_bu_identity
export nonuniform_three_point_derivative

@inline function _require_finite(name::AbstractString, value::Real)::Float64
    x = Float64(value)
    isfinite(x) || throw(ArgumentError("$(name) must be finite, got $(value)"))
    return x
end

function _finite_vector(name::AbstractString, values::AbstractVector)
    out = Float64.(values)
    all(isfinite, out) || throw(ArgumentError("$(name) must contain only finite values"))
    return out
end

"""
    bose_occupation(omega, mu, T) -> Float64

Return `1 / (exp((omega-mu)/T) - 1)` on the strict normal Bose domain.
Unlike exploratory cut prescriptions, this gate never moves an unsafe lower
bound above the pole silently.
"""
function bose_occupation(omega::Real, mu::Real, T::Real)::Float64
    omega_f = _require_finite("omega", omega)
    mu_f = _require_finite("mu", mu)
    T_f = _require_finite("T", T)
    T_f > 0.0 || throw(ArgumentError("T must be positive, got $(T)"))
    omega_f > mu_f || throw(ArgumentError(
        "Bose occupation requires omega > mu, got omega=$(omega), mu=$(mu)",
    ))
    return inv(expm1((omega_f - mu_f) / T_f))
end

"""
    nonuniform_three_point_derivative(x, y) -> Vector{Float64}

Differentiate tabulated values with the quadratic three-point formula on a
strictly increasing, potentially nonuniform grid. The two-point case falls
back to its unique secant slope.
"""
function nonuniform_three_point_derivative(x::AbstractVector, y::AbstractVector)
    x_f = _finite_vector("x", x)
    y_f = _finite_vector("y", y)
    n = length(x_f)
    length(y_f) == n || throw(ArgumentError("x/y length mismatch: $(n) != $(length(y_f))"))
    n >= 2 || throw(ArgumentError("need at least 2 points for derivative"))
    all(diff(x_f) .> 0.0) || throw(ArgumentError("x must be strictly increasing"))

    dy = similar(y_f)
    if n == 2
        slope = (y_f[2] - y_f[1]) / (x_f[2] - x_f[1])
        fill!(dy, slope)
        return dy
    end

    h1 = x_f[2] - x_f[1]
    h2 = x_f[3] - x_f[2]
    dy[1] = (
        -(2.0 * h1 + h2) * y_f[1] / (h1 * (h1 + h2)) +
        (h1 + h2) * y_f[2] / (h1 * h2) -
        h1 * y_f[3] / (h2 * (h1 + h2))
    )

    @inbounds for i in 2:(n - 1)
        h_left = x_f[i] - x_f[i - 1]
        h_right = x_f[i + 1] - x_f[i]
        dy[i] = (
            -h_right * y_f[i - 1] / (h_left * (h_left + h_right)) +
            (h_right - h_left) * y_f[i] / (h_left * h_right) +
            h_left * y_f[i + 1] / (h_right * (h_left + h_right))
        )
    end

    h1 = x_f[n - 1] - x_f[n - 2]
    h2 = x_f[n] - x_f[n - 1]
    dy[n] = (
        h2 * y_f[n - 2] / (h1 * (h1 + h2)) -
        (h1 + h2) * y_f[n - 1] / (h1 * h2) +
        (h1 + 2.0 * h2) * y_f[n] / (h2 * (h1 + h2))
    )
    return dy
end

"""
    finite_window_bu_identity(omega_nodes, omega_weights, F_values, dF_values;
                              T, mu=0, omega_min, omega_max, F_min, F_max)

Evaluate the finite-window identity

```text
integral(g dF/domega)/(2pi)
  = integral(g(1+g)F/T)/(2pi) + [gF]_{omega_min}^{omega_max}/(2pi).
```

`omega_nodes` and `omega_weights` describe one quadrature rule on the open
integration window. `F_min` and `F_max` must be evaluated on the same smooth
phase branch as `F_values`.
"""
function finite_window_bu_identity(
    omega_nodes::AbstractVector,
    omega_weights::AbstractVector,
    F_values::AbstractVector,
    dF_values::AbstractVector;
    T::Real,
    mu::Real=0.0,
    omega_min::Real,
    omega_max::Real,
    F_min::Real,
    F_max::Real,
)
    omega = _finite_vector("omega_nodes", omega_nodes)
    weights = _finite_vector("omega_weights", omega_weights)
    F = _finite_vector("F_values", F_values)
    dF = _finite_vector("dF_values", dF_values)

    n = length(omega)
    n > 0 || throw(ArgumentError("omega_nodes cannot be empty"))
    length(weights) == n || throw(ArgumentError(
        "omega_nodes/omega_weights length mismatch: $(n) != $(length(weights))",
    ))
    length(F) == n || throw(ArgumentError(
        "omega_nodes/F_values length mismatch: $(n) != $(length(F))",
    ))
    length(dF) == n || throw(ArgumentError(
        "omega_nodes/dF_values length mismatch: $(n) != $(length(dF))",
    ))

    T_f = _require_finite("T", T)
    mu_f = _require_finite("mu", mu)
    omega_min_f = _require_finite("omega_min", omega_min)
    omega_max_f = _require_finite("omega_max", omega_max)
    F_min_f = _require_finite("F_min", F_min)
    F_max_f = _require_finite("F_max", F_max)

    T_f > 0.0 || throw(ArgumentError("T must be positive, got $(T)"))
    omega_max_f > omega_min_f || throw(ArgumentError(
        "omega_max must exceed omega_min, got $(omega_max) <= $(omega_min)",
    ))
    omega_min_f > mu_f || throw(ArgumentError(
        "omega_min must exceed mu on the strict Bose domain, got omega_min=$(omega_min), mu=$(mu)",
    ))
    all(>(0.0), weights) || throw(ArgumentError("omega_weights must be positive"))
    issorted(omega) || throw(ArgumentError("omega_nodes must be sorted in ascending order"))
    all(x -> omega_min_f < x < omega_max_f, omega) || throw(ArgumentError(
        "omega_nodes must lie strictly inside (omega_min, omega_max)",
    ))

    derivative = 0.0
    bulk = 0.0
    @inbounds for i in eachindex(omega, weights, F, dF)
        g = bose_occupation(omega[i], mu_f, T_f)
        derivative += weights[i] * g * dF[i]
        bulk += weights[i] * g * (1.0 + g) * F[i] / T_f
    end
    derivative /= 2.0 * pi
    bulk /= 2.0 * pi

    g_min = bose_occupation(omega_min_f, mu_f, T_f)
    g_max = bose_occupation(omega_max_f, mu_f, T_f)
    boundary_lower = -g_min * F_min_f / (2.0 * pi)
    boundary_upper = g_max * F_max_f / (2.0 * pi)
    boundary = boundary_lower + boundary_upper
    total = bulk + boundary
    closure = derivative - total
    scale = max(abs(derivative), abs(total), eps(Float64))

    return (
        derivative=derivative,
        byparts_bulk=bulk,
        boundary_lower=boundary_lower,
        boundary_upper=boundary_upper,
        boundary=boundary,
        byparts_total=total,
        closure_abs=abs(closure),
        closure_rel=abs(closure) / scale,
        closure_signed=closure,
        bose_min=g_min,
        bose_max=g_max,
        weighted_phase_min=g_min * F_min_f,
        weighted_phase_max=g_max * F_max_f,
    )
end

end # module BUKernelGateUtils
