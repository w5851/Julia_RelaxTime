"""
    PrimaryStrategy

统一主策略配置：将主方法、多 seed 与 fallback 聚合到单一入口。
"""
Base.@kwdef struct PrimaryStrategy
    primary_method::Symbol = :newton
    use_fallback::Bool = true
    fallback_method::Symbol = :trust_region
    use_multiseed::Bool = false
    seed_strategy::Any = nothing
end

@inline function _resolve_primary_strategy_kwargs(kwargs)
    haskey(kwargs, :primary_strategy) || return kwargs

    strategy = kwargs[:primary_strategy]
    strategy isa PrimaryStrategy || throw(ArgumentError("primary_strategy must be PrimaryStrategy, got $(typeof(strategy))"))

    merged = Dict{Symbol,Any}(pairs(kwargs))
    delete!(merged, :primary_strategy)

    haskey(merged, :nlsolve_method) || (merged[:nlsolve_method] = strategy.primary_method)
    haskey(merged, :trust_region_fallback) || (merged[:trust_region_fallback] = strategy.use_fallback)
    haskey(merged, :fallback_method) || (merged[:fallback_method] = strategy.fallback_method)
    haskey(merged, :auto_multiseed_fallback) || (merged[:auto_multiseed_fallback] = strategy.use_multiseed)

    if strategy.seed_strategy !== nothing && !haskey(merged, :seed_strategy)
        merged[:seed_strategy] = strategy.seed_strategy
    end

    return (; pairs(merged)...)
end
