struct CEPResult
    found::Bool
    # Phase solvers work in the quark chemical potential mu_q. Keep the
    # historical field name for API compatibility; public artifacts should also
    # expose explicit muq_cep_MeV and muB_cep_MeV names.
    T_cep_MeV::Float64
    mu_cep_MeV::Float64
    uncertainty_T_MeV::Float64
    T_bracket_low_MeV::Float64
    T_bracket_high_MeV::Float64
    bracket_width_T_MeV::Float64
    result_status::Symbol
    T_last_first_order_MeV::Float64
    mu_last_first_order_MeV::Float64
    T_first_monotone_MeV::Float64
    ambiguity_width_T_MeV::Float64
    temperature_resolution_target_MeV::Float64
    eval_count::Int
    unknown_count::Int
    reason::Union{Nothing, String}
    method::Symbol
end

const CEP_RESULT_STATUSES = (:resolved, :ambiguous, :not_found)

"""Construct a CEP result while preserving the historical keyword contract.

When `result_status` is omitted, the old `found` flag is used to infer a
legacy resolved/not-found result. New callers should pass the explicit status;
ambiguous and not-found results never expose a midpoint CEP or a borrowed
Maxwell chemical potential.
"""
function CEPResult(;
        found::Bool=false,
        T_cep_MeV::Real=NaN,
        mu_cep_MeV::Real=NaN,
        uncertainty_T_MeV::Real=NaN,
        T_bracket_low_MeV::Real=NaN,
        T_bracket_high_MeV::Real=NaN,
        bracket_width_T_MeV::Real=NaN,
        result_status::Union{Nothing, Symbol}=nothing,
        T_last_first_order_MeV::Real=NaN,
        mu_last_first_order_MeV::Real=NaN,
        T_first_monotone_MeV::Real=NaN,
        ambiguity_width_T_MeV::Real=NaN,
        temperature_resolution_target_MeV::Real=NaN,
        eval_count::Int=0,
        unknown_count::Int=0,
        reason::Union{Nothing, String}=nothing,
        method::Symbol=:none)
    T_cep_MeV = Float64(T_cep_MeV)
    mu_cep_MeV = Float64(mu_cep_MeV)
    uncertainty_T_MeV = Float64(uncertainty_T_MeV)
    T_bracket_low_MeV = Float64(T_bracket_low_MeV)
    T_bracket_high_MeV = Float64(T_bracket_high_MeV)
    bracket_width_T_MeV = Float64(bracket_width_T_MeV)
    T_last_first_order_MeV = Float64(T_last_first_order_MeV)
    mu_last_first_order_MeV = Float64(mu_last_first_order_MeV)
    T_first_monotone_MeV = Float64(T_first_monotone_MeV)
    ambiguity_width_T_MeV = Float64(ambiguity_width_T_MeV)
    temperature_resolution_target_MeV = Float64(temperature_resolution_target_MeV)

    status = isnothing(result_status) ? (found ? :resolved : :not_found) : result_status
    status in CEP_RESULT_STATUSES || throw(ArgumentError(
        "result_status must be one of $(CEP_RESULT_STATUSES), got $(status)",
    ))
    eval_count >= 0 || throw(ArgumentError("eval_count must be nonnegative, got $(eval_count)"))
    unknown_count >= 0 || throw(ArgumentError("unknown_count must be nonnegative, got $(unknown_count)"))
    (isnan(temperature_resolution_target_MeV) ||
        (isfinite(temperature_resolution_target_MeV) && temperature_resolution_target_MeV > 0)) ||
        throw(ArgumentError(
            "temperature_resolution_target_MeV must be NaN or finite and positive, got $(temperature_resolution_target_MeV)",
        ))

    # Preserve old `CEPResult(found=true, ...)` construction verbatim. Any
    # explicit modern non-resolved status (and a legacy `found=false` result)
    # is canonicalized so downstream artifacts cannot publish a midpoint as
    # an ambiguous/not-found CEP.
    if status != :resolved && (!isnothing(result_status) || !found)
        found = false
        T_cep_MeV = NaN
        mu_cep_MeV = NaN
        uncertainty_T_MeV = NaN
    else
        if !isnothing(result_status) && status == :resolved
            (isfinite(T_cep_MeV) && isfinite(mu_cep_MeV)) || throw(ArgumentError(
                "resolved CEP results require finite T_cep_MeV and mu_cep_MeV",
            ))
            found = true
        else
            found = status == :resolved ? true : found
        end
    end

    if !isfinite(T_bracket_low_MeV) && isfinite(T_last_first_order_MeV)
        T_bracket_low_MeV = T_last_first_order_MeV
    end
    if !isfinite(T_bracket_high_MeV) && isfinite(T_first_monotone_MeV)
        T_bracket_high_MeV = T_first_monotone_MeV
    end
    if !isfinite(ambiguity_width_T_MeV) &&
       isfinite(T_bracket_low_MeV) && isfinite(T_bracket_high_MeV)
        ambiguity_width_T_MeV = T_bracket_high_MeV - T_bracket_low_MeV
    end
    if status == :ambiguous && isfinite(T_last_first_order_MeV) && isfinite(T_first_monotone_MeV) &&
       T_first_monotone_MeV < T_last_first_order_MeV
        throw(ArgumentError(
            "ambiguous CEP evidence interval must be ordered: " *
            "T_first_monotone_MeV >= T_last_first_order_MeV",
        ))
    end

    return CEPResult(
        found,
        T_cep_MeV,
        mu_cep_MeV,
        uncertainty_T_MeV,
        T_bracket_low_MeV,
        T_bracket_high_MeV,
        bracket_width_T_MeV,
        status,
        T_last_first_order_MeV,
        mu_last_first_order_MeV,
        T_first_monotone_MeV,
        ambiguity_width_T_MeV,
        temperature_resolution_target_MeV,
        eval_count,
        unknown_count,
        reason,
        method,
    )
end

Base.@kwdef struct FirstOrderSweepResult
    temperatures_MeV::Vector{Float64} = Float64[]
    statuses::Vector{Symbol} = Symbol[]
    mu_transitions_MeV::Vector{Float64} = Float64[]
    area_residuals::Vector{Float64} = Float64[]
    reasons::Vector{String} = String[]
    first_point_fallback::Bool = false
    fallback_start_T_MeV::Float64 = NaN
    unknown_count::Int = 0
    forced_invalid_count::Int = 0
end

"""Request-scoped verification contract for the opt-in hybrid rho policy.

The fields are deliberately versioned and diagnostic: they make the Stage-C
guard and the endpoint routes reproducible without exposing a new numerical
tolerance knob.  `:bounded_zero_density_v1` remains the legacy explicit route;
`:three_crossing_endpoint_local_v2` uses the actual right Maxwell crossing
bracket and active left-bracket midpoints.  The default production path does
not consume this contract.
"""
Base.@kwdef struct RhoHybridVerificationConfig
    local_step::Float64 = 0.003125
    targeted_cap::Int = 12
    guard_rule::Symbol = :extrema_outer_samples_v1
    comparison_epsilon::Float64 = 32eps(Float64)
    point_ranking_version::Symbol = :stage_b_features_v1
    candidate_policy::Symbol = :unique_three_crossing_topology_v1
    endpoint_policy::Symbol = :bounded_zero_density_v1
end

Base.@kwdef struct ProductionPipelineConfig
    T_start::Float64 = NaN
    T_end::Float64 = NaN
    dT_initial::Float64 = 5.0
    temperature_resolution_target_MeV::Float64 = 0.1
    # Deprecated compatibility alias; use temperature_resolution_target_MeV.
    cep_tol_MeV::Float64 = 0.1
    cep_max_bisect_iter::Int = 20
    area_tol_good::Float64 = 1e-4
    area_tol_bad::Float64 = 5e-4
    unknown_budget::Int = 5
    max_refine_level_rho::Int = 2
    adaptive_rho::Bool = true
    adaptive_slope_tol::Float64 = 5.0
    adaptive_min_gap::Float64 = 0.002
    adaptive_max_points::Int = 32
    adaptive_digits::Int = 6
    rho_geometry_convergence::Bool = true
    rho_position_tol_MeV::Float64 = 0.05
    rho_density_tol::Float64 = 0.005
    rho_maxwell_area_tol::Float64 = 1e-4
    rho_refinement_policy::Symbol = :uniform_nested
    rho_support_fine_step::Float64 = 0.025
    rho_support_targeted_cap::Int = 12
    rho_support_config::RhoSupportRefinement.RhoSupportConfig = RhoSupportRefinement.RhoSupportConfig()
    rho_hybrid_verification::RhoHybridVerificationConfig = RhoHybridVerificationConfig()
    adaptive_temperature::Bool = false
    temperature_max_refine_level::Int = 2
    temperature_position_tol_MeV::Float64 = 0.10
    temperature_density_tol::Float64 = 0.01
    temperature_maxwell_area_tol::Float64 = 1e-4
    crossover_T_max_MeV::Float64 = NaN
end

Base.@kwdef struct PromotionResult
    passed::Bool = false
    baseline_id::Union{Nothing, String} = nothing
    failed_checks::Vector{String} = String[]
    reference_dir::Union{Nothing, String} = nothing
end

Base.@kwdef struct PhasePipelineResult
    model_kind::Symbol = :PNJL
    model_variant::String = "default"
    xi::Float64 = 0.0
    run_id::String = ""
    cep::CEPResult = CEPResult()
    first_order_boundary::Vector{NamedTuple} = NamedTuple[]
    spinodal::Vector{NamedTuple} = NamedTuple[]
    crossover_line::Vector{NamedTuple} = NamedTuple[]
    diagnostics::Dict{String, Any} = Dict{String, Any}()
    config_snapshot::Dict{String, Any} = Dict{String, Any}()
    artifact_paths::Dict{String, String} = Dict{String, String}()
    promotion_status::PromotionResult = PromotionResult()
end
