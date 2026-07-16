module PhaseGuidedTransportScanPlan

struct PhaseGuidedPoint
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    mode::Symbol
    phase_reference_kind::Symbol
    scan_group::String
    group_label::String
    plot_panel::String
    plot_panel_label::String
    plot_series::String
    plot_series_label::String
    T_phase_base_MeV::Float64
    alpha_T::Float64
    phase_anchor_method::Symbol
    coexistence_side::Symbol
    coexistence_certified::Bool
    coexistence_delta_xi::Float64
    coexistence_T_lower_MeV::Float64
    coexistence_T_upper_MeV::Float64
    anchor_p_num::Int
    anchor_t_num::Int
    anchor_convergence_delta_MeV::Float64
    anchor_convergence_certified::Bool
end

struct PhaseGuidedPlan
    points::Vector{PhaseGuidedPoint}
    total::Int
end

const MODE_A_REFERENCE_XI = 0.0

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

    muq_MeV = muB_MeV / 3.0
    mus = copy(data.mu_values)
    Ts = copy(data.T_values)
    order = sortperm(mus)
    mus = mus[order]
    Ts = Ts[order]
    return _interp(muq_MeV, mus, Ts), xi_used
end

function _interpolate_crossover_temperature(muB_MeV::Float64, xi::Float64)
    Tc_MeV, xi_used = Main.GapTransportScanPhaseEquilibrium.interpolate_crossover_temperature(xi, muB_MeV / 3.0)
    return Tc_MeV, xi_used
end

function _cep_muB(xi::Float64)
    data, _ = _load_first_order_line(xi)
    return data.muB_CEP
end

function _phase_reference_for_mode_a(muB_MeV::Float64, xi::Float64)
    T_cross, _ = _interpolate_crossover_temperature(muB_MeV, MODE_A_REFERENCE_XI)
    T_first, _ = _interpolate_first_order_temperature(muB_MeV, MODE_A_REFERENCE_XI)
    mu_CEP = _cep_muB(MODE_A_REFERENCE_XI)

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
    mu_CEP = data.muB_CEP
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

function _mode_a_anchor(opts, muB_MeV::Float64)
    phase_reference_kind, interpolated_T_MeV = _phase_reference_for_mode_a(muB_MeV, MODE_A_REFERENCE_XI)
    isfinite(interpolated_T_MeV) || error("no phase reference temperature for mode a point: muB=$(muB_MeV)")

    if phase_reference_kind === :first_order && opts.phase_anchor_policy === :direct_coexistence
        numerics = (p_num=opts.p_num, t_num=opts.t_num)
        anchor = Main.GapTransportScanPhaseEquilibrium.direct_coexistence_anchor(
            muB_MeV,
            interpolated_T_MeV,
            numerics;
            xi=MODE_A_REFERENCE_XI,
        )
        certification = if any(isapprox(xi, 0.0; atol=1e-12, rtol=0.0) for xi in opts.xi_values) &&
                           any(isapprox(alpha, 1.0; atol=1e-12, rtol=0.0) for alpha in opts.alpha_T_values)
            Main.GapTransportScanPhaseEquilibrium.certify_coexistence_side_points(
                anchor,
                muB_MeV,
                numerics,
            )
        else
            nothing
        end
        return (
            phase_reference_kind=phase_reference_kind,
            T_phase_base_MeV=anchor.T_mid_MeV,
            phase_anchor_method=anchor.method,
            anchor=anchor,
            certification=certification,
        )
    end

    return (
        phase_reference_kind=phase_reference_kind,
        T_phase_base_MeV=interpolated_T_MeV,
        phase_anchor_method=:reference_interpolation,
        anchor=nothing,
        certification=nothing,
    )
end

function _mode_a_xi_values(opts, anchor_info, alpha_T::Float64)
    certification = anchor_info.certification
    exact_coexistence_slice = anchor_info.phase_reference_kind === :first_order &&
        anchor_info.phase_anchor_method === :direct_two_branch_equal_omega_bisection &&
        isapprox(alpha_T, 1.0; atol=1e-12, rtol=0.0) &&
        certification !== nothing && certification.certified
    exact_coexistence_slice || return opts.xi_values

    values = [xi for xi in opts.xi_values if !isapprox(xi, 0.0; atol=1e-12, rtol=0.0)]
    append!(values, (certification.minus_xi, certification.plus_xi))
    return unique(sort(Float64.(values)))
end

function build_plan(opts)
    points = PhaseGuidedPoint[]

    if opts.mode == :mode_a_fixed_muB_phase_scaled
        for muB_MeV in opts.muB_values
            anchor_info = _mode_a_anchor(opts, muB_MeV)
            for alpha_T in opts.alpha_T_values
                scan_group = "muB$(round(muB_MeV; digits=3))_alpha$(round(alpha_T; digits=3))"
                group_label = "muB=$(muB_MeV) MeV, alpha_T=$(alpha_T)"
                plot_panel = "muB$(round(muB_MeV; digits=3))"
                plot_panel_label = "muB=$(muB_MeV) MeV"
                plot_series = "alpha$(round(alpha_T; digits=3))"
                T_phase_base_MeV = anchor_info.T_phase_base_MeV
                phase_reference_kind = anchor_info.phase_reference_kind
                for xi in _mode_a_xi_values(opts, anchor_info, alpha_T)
                    plot_series_label = "alpha_T=$(alpha_T), T=$(round(alpha_T * T_phase_base_MeV; digits=3)) MeV"
                    certification = anchor_info.certification
                    coexistence_side = if certification !== nothing && isapprox(xi, certification.minus_xi; atol=1e-12, rtol=0.0)
                        :quark_side
                    elseif certification !== nothing && isapprox(xi, certification.plus_xi; atol=1e-12, rtol=0.0)
                        :hadron_side
                    else
                        :none
                    end
                    anchor = anchor_info.anchor
                    push!(points, PhaseGuidedPoint(
                        alpha_T * T_phase_base_MeV,
                        muB_MeV,
                        xi,
                        opts.mode,
                        phase_reference_kind,
                        scan_group,
                        group_label,
                        plot_panel,
                        plot_panel_label,
                        plot_series,
                        plot_series_label,
                        T_phase_base_MeV,
                        alpha_T,
                        anchor_info.phase_anchor_method,
                        coexistence_side,
                        certification !== nothing && certification.certified,
                        certification === nothing ? NaN : certification.delta_xi,
                        anchor === nothing ? NaN : anchor.T_lower_MeV,
                        anchor === nothing ? NaN : anchor.T_upper_MeV,
                        opts.p_num,
                        opts.t_num,
                        certification === nothing ? NaN : certification.anchor_convergence_delta_MeV,
                        certification !== nothing && certification.anchor_convergence_certified,
                    ))
                end
            end
        end
    elseif opts.mode == :mode_b_fixed_T_sparse_muB
        for T_MeV in opts.T_values
            for muB_MeV in opts.muB_values
                scan_group = "T$(round(T_MeV; digits=3))_muB$(round(muB_MeV; digits=3))"
                group_label = "T=$(T_MeV) MeV, muB=$(muB_MeV) MeV"
                plot_panel = "T$(round(T_MeV; digits=3))"
                plot_panel_label = "T=$(T_MeV) MeV"
                plot_series = "muB$(round(muB_MeV; digits=3))"
                plot_series_label = "muB=$(muB_MeV) MeV"
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
                        plot_panel,
                        plot_panel_label,
                        plot_series,
                        plot_series_label,
                        T_MeV,
                        NaN,
                        :reference_classification,
                        :none,
                        false,
                        NaN,
                        NaN,
                        NaN,
                        opts.p_num,
                        opts.t_num,
                        NaN,
                        false,
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
