Base.@kwdef struct CEPResult
    found::Bool = false
    T_cep_MeV::Float64 = NaN
    mu_cep_MeV::Float64 = NaN
    uncertainty_T_MeV::Float64 = NaN
    eval_count::Int = 0
    unknown_count::Int = 0
    reason::Union{Nothing, String} = nothing
    method::Symbol = :none
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
