Base.@kwdef struct PhaseGeometryTolerances
    position_MeV::Float64 = 0.05
    density::Float64 = 0.005
    maxwell_area::Float64 = 1e-4
    response_rtol::Float64 = 0.05
end

Base.@kwdef struct PhaseGeometryError
    comparable::Bool = false
    converged::Bool = false
    position_MeV::Float64 = Inf
    density::Float64 = Inf
    maxwell_area::Float64 = Inf
    response_rtol::Float64 = Inf
    reason::String = "not_compared"
end

function _validate_phase_geometry_tolerances(tol::PhaseGeometryTolerances)
    isfinite(tol.position_MeV) && tol.position_MeV > 0 ||
        throw(ArgumentError("phase position tolerance must be finite and positive, got $(tol.position_MeV)"))
    isfinite(tol.density) && tol.density > 0 ||
        throw(ArgumentError("phase density tolerance must be finite and positive, got $(tol.density)"))
    isfinite(tol.maxwell_area) && tol.maxwell_area > 0 ||
        throw(ArgumentError("Maxwell area tolerance must be finite and positive, got $(tol.maxwell_area)"))
    isfinite(tol.response_rtol) && tol.response_rtol > 0 ||
        throw(ArgumentError("phase response relative tolerance must be finite and positive, got $(tol.response_rtol)"))
    return tol
end

@inline function _phase_max_abs_difference(left::Tuple, right::Tuple)
    length(left) == length(right) || return Inf
    isempty(left) && return 0.0
    err = 0.0
    for i in eachindex(left)
        a = Float64(left[i])
        b = Float64(right[i])
        (isfinite(a) && isfinite(b)) || return Inf
        err = max(err, abs(a - b))
    end
    return err
end

@inline function _phase_response_relative_error(actual::Real, expected::Real)
    a = Float64(actual)
    e = Float64(expected)
    (isfinite(a) && isfinite(e)) || return Inf
    scale = max(abs(a), abs(e), 1e-12)
    return abs(a - e) / scale
end

function _phase_geometry_snapshot(result)
    status = Symbol(result.status)
    reason = String(result.reason)
    if status == :valid
        sres = result.sres
        position = (
            Float64(something(result.mu_transition, NaN)),
            Float64(something(sres.mu_spinodal_hadron, NaN)),
            Float64(something(sres.mu_spinodal_quark, NaN)),
        )
        density = (
            Float64(something(result.rho_hadron, NaN)),
            Float64(something(result.rho_quark, NaN)),
            Float64(something(sres.rho_spinodal_hadron, NaN)),
            Float64(something(sres.rho_spinodal_quark, NaN)),
        )
        return (
            status=status,
            reason=reason,
            position=position,
            density=density,
            maxwell_area=Float64(something(result.area_residual, Inf)),
            derivative_sign_changes=Int(sres.derivative_sign_changes),
        )
    end
    return (
        status=status,
        reason=reason,
        position=(),
        density=(),
        maxwell_area=Inf,
        derivative_sign_changes=Int(result.sres.derivative_sign_changes),
    )
end

function _compare_phase_geometry(coarse, fine, tol::PhaseGeometryTolerances)
    _validate_phase_geometry_tolerances(tol)
    a = _phase_geometry_snapshot(coarse)
    b = _phase_geometry_snapshot(fine)

    if a.status == :invalid && b.status == :invalid && a.reason == "no_s_shape" && b.reason == "no_s_shape"
        return PhaseGeometryError(
            comparable=true,
            converged=true,
            position_MeV=0.0,
            density=0.0,
            maxwell_area=0.0,
            response_rtol=0.0,
            reason="stable_no_s_shape",
        )
    end

    if a.status != :valid || b.status != :valid
        return PhaseGeometryError(reason="classification_changed_or_unresolved:$(a.status)->$(b.status)")
    end

    position_error = _phase_max_abs_difference(a.position, b.position)
    density_error = _phase_max_abs_difference(a.density, b.density)
    area_gate_value = max(a.maxwell_area, b.maxwell_area)
    signs_stable = a.derivative_sign_changes == b.derivative_sign_changes
    converged = signs_stable &&
                position_error <= tol.position_MeV &&
                density_error <= tol.density &&
                area_gate_value <= tol.maxwell_area
    return PhaseGeometryError(
        comparable=true,
        converged=converged,
        position_MeV=position_error,
        density=density_error,
        maxwell_area=area_gate_value,
        response_rtol=0.0,
        reason=signs_stable ? (converged ? "converged" : "geometry_tolerance_exceeded") : "spinodal_topology_changed",
    )
end

function _phase_geometry_midpoint_error(left, midpoint, right, tol::PhaseGeometryTolerances)
    _validate_phase_geometry_tolerances(tol)
    a = _phase_geometry_snapshot(left)
    m = _phase_geometry_snapshot(midpoint)
    b = _phase_geometry_snapshot(right)

    if all(x -> x.status == :invalid && x.reason == "no_s_shape", (a, m, b))
        return PhaseGeometryError(
            comparable=true,
            converged=true,
            position_MeV=0.0,
            density=0.0,
            maxwell_area=0.0,
            response_rtol=0.0,
            reason="stable_no_s_shape",
        )
    end
    if any(x -> x.status != :valid, (a, m, b))
        return PhaseGeometryError(reason="midpoint_classification_changed_or_unresolved:$(a.status),$(m.status),$(b.status)")
    end

    position_error = 0.0
    for i in eachindex(m.position)
        position_error = max(position_error,
            _linear_midpoint_error(m.position[i], a.position[i], b.position[i]))
    end
    density_error = 0.0
    for i in eachindex(m.density)
        density_error = max(density_error,
            _linear_midpoint_error(m.density[i], a.density[i], b.density[i]))
    end
    signs_stable = a.derivative_sign_changes == m.derivative_sign_changes == b.derivative_sign_changes
    area_gate_value = m.maxwell_area
    converged = signs_stable &&
                position_error <= tol.position_MeV &&
                density_error <= tol.density &&
                area_gate_value <= tol.maxwell_area
    return PhaseGeometryError(
        comparable=true,
        converged=converged,
        position_MeV=position_error,
        density=density_error,
        maxwell_area=area_gate_value,
        response_rtol=0.0,
        reason=signs_stable ? (converged ? "converged" : "interpolation_tolerance_exceeded") : "spinodal_topology_changed",
    )
end

@inline function _linear_midpoint_error(actual::Real, left::Real, right::Real)
    a = Float64(actual)
    l = Float64(left)
    r = Float64(right)
    (isfinite(a) && isfinite(l) && isfinite(r)) || return Inf
    return abs(a - 0.5 * (l + r))
end

function _rows_by_key(rows, key::Symbol; digits::Int=8)
    mapped = Dict{Float64, Any}()
    for row in rows
        value = Float64(getproperty(row, key))
        isfinite(value) || continue
        mapped[round(value; digits=digits)] = row
    end
    return mapped
end

function _phase_result_midpoint_error(
        left::PhasePipelineResult,
        midpoint::PhasePipelineResult,
        right::PhasePipelineResult,
        tol::PhaseGeometryTolerances)
    _validate_phase_geometry_tolerances(tol)
    position_error = 0.0
    density_error = 0.0
    maxwell_area = 0.0
    response_error = 0.0
    comparable_count = 0

    boundary_left = _rows_by_key(left.first_order_boundary, :T_MeV)
    boundary_mid = _rows_by_key(midpoint.first_order_boundary, :T_MeV)
    boundary_right = _rows_by_key(right.first_order_boundary, :T_MeV)
    for key in intersect(intersect(keys(boundary_left), keys(boundary_mid)), keys(boundary_right))
        a, m, b = boundary_left[key], boundary_mid[key], boundary_right[key]
        position_error = max(position_error,
            _linear_midpoint_error(m.mu_transition_MeV, a.mu_transition_MeV, b.mu_transition_MeV))
        density_error = max(density_error,
            _linear_midpoint_error(m.rho_hadron, a.rho_hadron, b.rho_hadron),
            _linear_midpoint_error(m.rho_quark, a.rho_quark, b.rho_quark))
        maxwell_area = max(maxwell_area, Float64(m.area_residual))
        comparable_count += 1
    end

    spinodal_left = _rows_by_key(left.spinodal, :T_MeV)
    spinodal_mid = _rows_by_key(midpoint.spinodal, :T_MeV)
    spinodal_right = _rows_by_key(right.spinodal, :T_MeV)
    for key in intersect(intersect(keys(spinodal_left), keys(spinodal_mid)), keys(spinodal_right))
        a, m, b = spinodal_left[key], spinodal_mid[key], spinodal_right[key]
        position_error = max(position_error,
            _linear_midpoint_error(m.mu_spinodal_hadron_MeV, a.mu_spinodal_hadron_MeV, b.mu_spinodal_hadron_MeV),
            _linear_midpoint_error(m.mu_spinodal_quark_MeV, a.mu_spinodal_quark_MeV, b.mu_spinodal_quark_MeV))
        density_error = max(density_error,
            _linear_midpoint_error(m.rho_spinodal_hadron, a.rho_spinodal_hadron, b.rho_spinodal_hadron),
            _linear_midpoint_error(m.rho_spinodal_quark, a.rho_spinodal_quark, b.rho_spinodal_quark))
        comparable_count += 1
    end

    crossover_left = _rows_by_key(left.crossover_line, :mu_MeV)
    crossover_mid = _rows_by_key(midpoint.crossover_line, :mu_MeV)
    crossover_right = _rows_by_key(right.crossover_line, :mu_MeV)
    for key in intersect(intersect(keys(crossover_left), keys(crossover_mid)), keys(crossover_right))
        a, m, b = crossover_left[key], crossover_mid[key], crossover_right[key]
        position_error = max(position_error,
            _linear_midpoint_error(m.T_crossover_MeV, a.T_crossover_MeV, b.T_crossover_MeV))
        density_error = max(density_error,
            _linear_midpoint_error(m.rho, a.rho, b.rho))
        expected_response = 0.5 * (Float64(a.derivative) + Float64(b.derivative))
        response_error = max(response_error,
            _phase_response_relative_error(m.derivative, expected_response))
        comparable_count += 1
    end

    if left.cep.found && midpoint.cep.found && right.cep.found
        position_error = max(position_error,
            _linear_midpoint_error(midpoint.cep.T_cep_MeV, left.cep.T_cep_MeV, right.cep.T_cep_MeV),
            _linear_midpoint_error(midpoint.cep.mu_cep_MeV, left.cep.mu_cep_MeV, right.cep.mu_cep_MeV))
        comparable_count += 1
    end

    if comparable_count == 0
        all_empty = all(result -> isempty(result.first_order_boundary) && isempty(result.spinodal) &&
                                  isempty(result.crossover_line) && !result.cep.found,
                        (left, midpoint, right))
        return PhaseGeometryError(
            comparable=all_empty,
            converged=all_empty,
            position_MeV=all_empty ? 0.0 : Inf,
            density=all_empty ? 0.0 : Inf,
            maxwell_area=all_empty ? 0.0 : Inf,
            response_rtol=all_empty ? 0.0 : Inf,
            reason=all_empty ? "stable_no_transition_signal" : "no_common_phase_geometry",
        )
    end

    converged = position_error <= tol.position_MeV &&
                density_error <= tol.density &&
                maxwell_area <= tol.maxwell_area &&
                response_error <= tol.response_rtol
    return PhaseGeometryError(
        comparable=true,
        converged=converged,
        position_MeV=position_error,
        density=density_error,
        maxwell_area=maxwell_area,
        response_rtol=response_error,
        reason=converged ? "converged" : "interpolation_tolerance_exceeded",
    )
end
