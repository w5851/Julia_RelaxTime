module PhaseGuidedTransportScanPlan

struct PhaseGuidedPoint
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    mode::Symbol
    phase_reference_kind::Symbol
    scan_group::String
    group_label::String
    T_phase_base_MeV::Float64
    alpha_T::Float64
end

struct PhaseGuidedPlan
    points::Vector{PhaseGuidedPoint}
    total::Int
end

function _interp(x::Float64, xs::Vector{Float64}, ys::Vector{Float64})
    isempty(xs) && return NaN
    length(xs) == 1 && return ys[1]
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    for i in 1:(length(xs) - 1)
        left = xs[i]
        right = xs[i + 1]
        if left <= x <= right
            weight = (x - left) / (right - left)
            return ys[i] + weight * (ys[i + 1] - ys[i])
        end
    end
    return NaN
end

function _load_first_order_line(xi::Float64)
    data = Main.GapTransportScanPhaseEquilibrium.load_phase_boundary_data(xi)
    return data, data.xi
end

function _interpolate_first_order_temperature(muB_MeV::Float64, xi::Float64)
    data, xi_used = _load_first_order_line(xi)
    isempty(data.mu_values) && return NaN, xi_used

    mus = copy(data.mu_values)
    Ts = copy(data.T_values)
    order = sortperm(mus)
    mus = mus[order]
    Ts = Ts[order]
    return _interp(muB_MeV, mus, Ts), xi_used
end

function _interpolate_crossover_temperature(muB_MeV::Float64, xi::Float64)
    Tc_MeV, xi_used = Main.GapTransportScanPhaseEquilibrium.interpolate_crossover_temperature(xi, muB_MeV / 3.0)
    return Tc_MeV, xi_used
end

function _cep_muB(xi::Float64)
    data, _ = _load_first_order_line(xi)
    return data.mu_CEP
end

function _phase_reference_for_mode_a(muB_MeV::Float64, xi::Float64)
    T_cross, _ = _interpolate_crossover_temperature(muB_MeV, xi)
    T_first, _ = _interpolate_first_order_temperature(muB_MeV, xi)
    mu_CEP = _cep_muB(xi)

    if isfinite(mu_CEP) && muB_MeV > mu_CEP + 1e-6 && isfinite(T_first)
        return (:first_order, T_first)
    elseif isfinite(T_cross)
        return (:crossover, T_cross)
    elseif isfinite(T_first)
        return (:first_order, T_first)
    else
        return (:unknown, NaN)
    end
end

function _phase_reference_for_mode_b(T_MeV::Float64, muB_MeV::Float64, xi::Float64)
    data, _ = _load_first_order_line(xi)
    mu_CEP = data.mu_CEP
    T_CEP = data.T_CEP
    T_first, _ = _interpolate_first_order_temperature(muB_MeV, xi)
    if isfinite(T_CEP) && isfinite(mu_CEP) && abs(T_MeV - T_CEP) <= 8.0 && abs(muB_MeV - mu_CEP) <= 40.0
        return :cep_neighbor
    elseif isfinite(T_first) && T_MeV <= T_first && isfinite(mu_CEP) && muB_MeV >= mu_CEP
        return :first_order
    else
        return :crossover
    end
end

function build_plan(opts)
    points = PhaseGuidedPoint[]

    if opts.mode == :mode_a_fixed_muB_phase_scaled
        for muB_MeV in opts.muB_values
            for alpha_T in opts.alpha_T_values
                scan_group = "muB$(round(muB_MeV; digits=3))_alpha$(round(alpha_T; digits=3))"
                group_label = "muB=$(muB_MeV) MeV, alpha_T=$(alpha_T)"
                for xi in opts.xi_values
                    phase_reference_kind, T_phase_base_MeV = _phase_reference_for_mode_a(muB_MeV, xi)
                    isfinite(T_phase_base_MeV) || error("no phase reference temperature for mode a point: muB=$(muB_MeV), xi=$(xi)")
                    push!(points, PhaseGuidedPoint(
                        alpha_T * T_phase_base_MeV,
                        muB_MeV,
                        xi,
                        opts.mode,
                        phase_reference_kind,
                        scan_group,
                        group_label,
                        T_phase_base_MeV,
                        alpha_T,
                    ))
                end
            end
        end
    elseif opts.mode == :mode_b_fixed_T_sparse_muB
        for T_MeV in opts.T_values
            for muB_MeV in opts.muB_values
                scan_group = "T$(round(T_MeV; digits=3))_muB$(round(muB_MeV; digits=3))"
                group_label = "T=$(T_MeV) MeV, muB=$(muB_MeV) MeV"
                for xi in opts.xi_values
                    phase_reference_kind = _phase_reference_for_mode_b(T_MeV, muB_MeV, xi)
                    push!(points, PhaseGuidedPoint(
                        T_MeV,
                        muB_MeV,
                        xi,
                        opts.mode,
                        phase_reference_kind,
                        scan_group,
                        group_label,
                        T_MeV,
                        NaN,
                    ))
                end
            end
        end
    else
        error("unsupported mode: $(opts.mode)")
    end

    return PhaseGuidedPlan(points, length(points))
end

function execute_plan!(point_runner!::Function, opts, plan, existing)
    stats_success = 0
    stats_error = 0
    stats_skipped = 0
    done = 0

    grouped = Dict{String, Vector{PhaseGuidedPoint}}()
    group_order = String[]
    for point in plan.points
        if !haskey(grouped, point.scan_group)
            grouped[point.scan_group] = PhaseGuidedPoint[]
            push!(group_order, point.scan_group)
        end
        push!(grouped[point.scan_group], point)
    end

    for group in group_order
        previous_solution = nothing
        previous_phase = :unknown
        for point in grouped[group]
            done += 1
            key = (point.T_MeV, point.muB_MeV, point.xi)
            if opts.resume && (key in existing)
                stats_skipped += 1
                continue
            end

            point_result = point_runner!(point, previous_solution, previous_phase)
            previous_solution = point_result.next_solution
            previous_phase = point_result.next_phase
            if point_result.success
                stats_success += 1
            else
                stats_error += 1
            end
        end
    end

    return (
        success=stats_success,
        error=stats_error,
        skipped=stats_skipped,
        done=done,
        total=plan.total,
    )
end

export PhaseGuidedPoint
export PhaseGuidedPlan
export build_plan
export execute_plan!

end # module PhaseGuidedTransportScanPlan
