Base.@kwdef struct CEPResult
    found::Bool = false
    T_cep_MeV::Float64 = NaN
    # Phase solvers work in the quark chemical potential mu_q. Keep the
    # historical field name for API compatibility; public artifacts should also
    # expose explicit muq_cep_MeV and muB_cep_MeV names.
    mu_cep_MeV::Float64 = NaN
    uncertainty_T_MeV::Float64 = NaN
    eval_count::Int = 0
    unknown_count::Int = 0
    reason::Union{Nothing, String} = nothing
    method::Symbol = :none
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

Base.@kwdef struct ProductionPipelineConfig
    T_start::Float64 = NaN
    T_end::Float64 = NaN
    dT_initial::Float64 = 5.0
    cep_tol_MeV::Float64 = 0.5
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
