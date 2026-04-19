struct DiffTarget{F}
    name::Symbol
    evaluator::F
end

@inline DiffTarget(name::Symbol) = DiffTarget(name, _ctx -> throw(ErrorException("DiffTarget $(name) evaluator is not implemented")))

const _DIFF_TARGET_REGISTRY = Dict{Symbol, DiffTarget}(
    :pressure => DiffTarget(:pressure),
    :entropy => DiffTarget(:entropy),
    :rho_norm => DiffTarget(:rho_norm),
    :energy => DiffTarget(:energy),
)

@inline function diff_target(name::Symbol)
    haskey(_DIFF_TARGET_REGISTRY, name) || throw(ArgumentError("unknown diff target: $(name)"))
    return _DIFF_TARGET_REGISTRY[name]
end
