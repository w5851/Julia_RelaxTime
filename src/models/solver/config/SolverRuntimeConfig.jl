"""SolverRuntimeConfig

运行时配置类型：
- 先覆盖 FixedRho 与 FixedEntropy 关键路径；
- 入口解析时完成类型校验与默认值折叠，减少 Dict 管道分叉。
"""

struct FixedRhoRuntimeConfig
    seed_guess::Vector{Float64}
    seed_candidates::Vector{Vector{Float64}}
    primary_method::Symbol
    trust_region_fallback::Bool
    fallback_method::Symbol
    continuity_seed::Bool
    evaluate_all_attempts::Bool
    fixedrho_joint_solve::Bool
    residual_norm_max::Float64
    xi::Float64
    p_num::Int
    t_num::Int
end

struct FixedEntropyRuntimeConfig
    seed_guess::Vector{Float64}
    seed_candidates::Vector{Vector{Float64}}
    primary_method::Symbol
    trust_region_fallback::Bool
    fallback_method::Symbol
    continuity_seed::Bool
    evaluate_all_attempts::Bool
    residual_norm_max::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    rho0::Float64
    mu0::Float64
end

struct FixedSigmaRuntimeConfig
    seed_guess::Vector{Float64}
    seed_candidates::Vector{Vector{Float64}}
    primary_method::Symbol
    trust_region_fallback::Bool
    fallback_method::Symbol
    continuity_seed::Bool
    evaluate_all_attempts::Bool
    residual_norm_max::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    rho0::Float64
    mu0::Float64
end

struct FixedAsymmetricRhoRuntimeConfig
    seed_guess::Vector{Float64}
    seed_candidates::Vector{Vector{Float64}}
    primary_method::Symbol
    trust_region_fallback::Bool
    fallback_method::Symbol
    continuity_seed::Bool
    evaluate_all_attempts::Bool
    residual_norm_max::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    rho0::Float64
    mu0::Vector{Float64}
end

@inline function _runtime_bool_option(kwargs::Dict{Symbol,Any}, key::Symbol, default::Bool)::Bool
    value = get(kwargs, key, default)
    value isa Bool || throw(ArgumentError("$(key) must be Bool, got $(typeof(value))"))
    return value
end

@inline function _runtime_symbol_option(kwargs::Dict{Symbol,Any}, key::Symbol, default::Symbol)::Symbol
    value = get(kwargs, key, default)
    value isa Symbol || throw(ArgumentError("$(key) must be Symbol, got $(typeof(value))"))
    return value
end

@inline function _runtime_positive_real_option(kwargs::Dict{Symbol,Any}, key::Symbol, default::Real)::Float64
    value = get(kwargs, key, default)
    value isa Real || throw(ArgumentError("$(key) must be Real, got $(typeof(value))"))
    value > 0 || throw(ArgumentError("$(key) must be positive, got $(value)"))
    return Float64(value)
end

@inline function _runtime_positive_int_option(kwargs::Dict{Symbol,Any}, key::Symbol, default::Int)::Int
    value = get(kwargs, key, default)
    value isa Integer || throw(ArgumentError("$(key) must be Integer, got $(typeof(value))"))
    value > 0 || throw(ArgumentError("$(key) must be positive, got $(value)"))
    return Int(value)
end

@inline function _runtime_real_option(kwargs::Dict{Symbol,Any}, key::Symbol, default::Real)::Float64
    value = get(kwargs, key, default)
    value isa Real || throw(ArgumentError("$(key) must be Real, got $(typeof(value))"))
    return Float64(value)
end

function _fixedrho_runtime_config_from_kwargs(mode::FixedRho, kwargs::Dict{Symbol,Any})
    _ = mode
    seed_guess_raw = get(kwargs, :seed_guess, nothing)
    seed_guess_raw === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedRho forward_solve"))
    seed_guess = Float64.(seed_guess_raw)
    seed_candidates = _resolve_optional_seed_candidates(kwargs)

    primary_method = if haskey(kwargs, :nlsolve_method)
        _runtime_symbol_option(kwargs, :nlsolve_method, :trust_region)
    else
        _runtime_bool_option(kwargs, :continuity_seed, false) ? :newton : :trust_region
    end

    return FixedRhoRuntimeConfig(
        seed_guess,
        seed_candidates,
        primary_method,
        _runtime_bool_option(kwargs, :trust_region_fallback, true),
        _runtime_symbol_option(kwargs, :fallback_method, :trust_region),
        _runtime_bool_option(kwargs, :continuity_seed, false),
        _runtime_bool_option(kwargs, :evaluate_all_attempts, false),
        _runtime_bool_option(kwargs, :fixedrho_joint_solve, true),
        _runtime_positive_real_option(kwargs, :residual_norm_max, 1e-6),
        _runtime_real_option(kwargs, :xi, 0.0),
        _runtime_positive_int_option(kwargs, :p_num, 24),
        _runtime_positive_int_option(kwargs, :t_num, 8),
    )
end

function _fixedentropy_runtime_config_from_kwargs(mode::FixedEntropy, kwargs::Dict{Symbol,Any})
    _ = mode
    seed_guess_raw = get(kwargs, :seed_guess, nothing)
    seed_guess_raw === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedEntropy forward_solve"))
    seed_guess = Float64.(seed_guess_raw)
    seed_candidates = _resolve_optional_seed_candidates(kwargs)

    primary_method = if haskey(kwargs, :nlsolve_method)
        _runtime_symbol_option(kwargs, :nlsolve_method, :trust_region)
    else
        _runtime_bool_option(kwargs, :continuity_seed, false) ? :newton : :trust_region
    end

    mu0 = if haskey(kwargs, :mu0)
        Float64(kwargs[:mu0])
    elseif length(seed_guess) >= 8
        mean(seed_guess[6:8])
    else
        0.2
    end

    return FixedEntropyRuntimeConfig(
        seed_guess,
        seed_candidates,
        primary_method,
        _runtime_bool_option(kwargs, :trust_region_fallback, true),
        _runtime_symbol_option(kwargs, :fallback_method, :trust_region),
        _runtime_bool_option(kwargs, :continuity_seed, false),
        _runtime_bool_option(kwargs, :evaluate_all_attempts, false),
        _runtime_positive_real_option(kwargs, :residual_norm_max, 1e-6),
        _runtime_real_option(kwargs, :xi, 0.0),
        _runtime_positive_int_option(kwargs, :p_num, 24),
        _runtime_positive_int_option(kwargs, :t_num, 8),
        _runtime_real_option(kwargs, :rho0, rho0),
        mu0,
    )
end

function _fixedsigma_runtime_config_from_kwargs(mode::FixedSigma, kwargs::Dict{Symbol,Any})
    _ = mode
    seed_guess_raw = get(kwargs, :seed_guess, nothing)
    seed_guess_raw === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedSigma forward_solve"))
    seed_guess = Float64.(seed_guess_raw)
    seed_candidates = _resolve_optional_seed_candidates(kwargs)

    primary_method = if haskey(kwargs, :nlsolve_method)
        _runtime_symbol_option(kwargs, :nlsolve_method, :trust_region)
    else
        _runtime_bool_option(kwargs, :continuity_seed, false) ? :newton : :trust_region
    end

    mu0 = if haskey(kwargs, :mu0)
        Float64(kwargs[:mu0])
    elseif length(seed_guess) >= 8
        mean(seed_guess[6:8])
    else
        0.2
    end

    return FixedSigmaRuntimeConfig(
        seed_guess,
        seed_candidates,
        primary_method,
        _runtime_bool_option(kwargs, :trust_region_fallback, true),
        _runtime_symbol_option(kwargs, :fallback_method, :trust_region),
        _runtime_bool_option(kwargs, :continuity_seed, false),
        _runtime_bool_option(kwargs, :evaluate_all_attempts, false),
        _runtime_positive_real_option(kwargs, :residual_norm_max, 1e-6),
        _runtime_real_option(kwargs, :xi, 0.0),
        _runtime_positive_int_option(kwargs, :p_num, 24),
        _runtime_positive_int_option(kwargs, :t_num, 8),
        _runtime_real_option(kwargs, :rho0, rho0),
        mu0,
    )
end

function _fixedasymrho_runtime_config_from_kwargs(mode::FixedAsymmetricRho, kwargs::Dict{Symbol,Any})
    _ = mode
    seed_guess_raw = get(kwargs, :seed_guess, nothing)
    seed_guess_raw === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedAsymmetricRho forward_solve"))
    seed_guess = Float64.(seed_guess_raw)
    seed_candidates = _resolve_optional_seed_candidates(kwargs)

    primary_method = if haskey(kwargs, :nlsolve_method)
        _runtime_symbol_option(kwargs, :nlsolve_method, :trust_region)
    else
        _runtime_bool_option(kwargs, :continuity_seed, false) ? :newton : :trust_region
    end

    mu0 = if haskey(kwargs, :mu0)
        Float64.(kwargs[:mu0])
    elseif length(seed_guess) >= 8
        Float64.(seed_guess[6:8])
    else
        Float64[0.2, 0.2, 0.2]
    end

    return FixedAsymmetricRhoRuntimeConfig(
        seed_guess,
        seed_candidates,
        primary_method,
        _runtime_bool_option(kwargs, :trust_region_fallback, true),
        _runtime_symbol_option(kwargs, :fallback_method, :trust_region),
        _runtime_bool_option(kwargs, :continuity_seed, false),
        _runtime_bool_option(kwargs, :evaluate_all_attempts, false),
        _runtime_positive_real_option(kwargs, :residual_norm_max, 1e-6),
        _runtime_real_option(kwargs, :xi, 0.0),
        _runtime_positive_int_option(kwargs, :p_num, 24),
        _runtime_positive_int_option(kwargs, :t_num, 8),
        _runtime_real_option(kwargs, :rho0, rho0),
        mu0,
    )
end
