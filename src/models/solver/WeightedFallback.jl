"""
    solve_weighted_block_fallback(mode::FixedAsymmetricRho, T_fm; initial_seed, kwargs...) -> Union{SolverResult, Nothing}

Hybrid 连续性链路失败时的 weighted-block 兜底。
"""

@inline function _weighted_stage_schedule()
    return (
        (w6=0.10, w7=0.10, w8=1.0, method=:trust_region, iterations=600),
        (w6=0.50, w7=0.50, w8=1.0, method=:trust_region, iterations=1000),
        (w6=1.00, w7=1.00, w8=1.0, method=:trust_region, iterations=1400),
        (w6=1.00, w7=1.00, w8=1.0, method=:newton,       iterations=1000),
    )
end

function _run_weighted_stage(raw_residual!::Function, x0::Vector{Float64};
    w6::Float64,
    w7::Float64,
    w8::Float64,
    method::Symbol,
    iterations::Int,
)
    weights = Float64[1.0, 1.0, 1.0, 1.0, 1.0, w6, w7, w8]

    weighted_residual! = (F, x) -> begin
        rawF = similar(F)
        raw_residual!(rawF, x)
        @inbounds for i in 1:8
            F[i] = weights[i] * rawF[i]
        end
        return nothing
    end

    res = nlsolve(weighted_residual!, x0;
        autodiff=:forward,
        method=method,
        xtol=1e-9,
        ftol=1e-9,
        iterations=iterations,
    )

    x = Vector{Float64}(res.zero)
    rawF = zeros(8)
    raw_residual!(rawF, x)
    raw_norm = sqrt(sum(abs2, rawF))

    return (x=x, raw_norm=raw_norm)
end

function solve_weighted_block_fallback(mode::FixedAsymmetricRho, T_fm::Real;
    initial_seed::AbstractVector{<:Real},
    max_seed_candidates::Int=3,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    model_kind::Symbol=:PNJL,
    residual_norm_max::Real=1e-6,
    nlsolve_kwargs...)

    base_seed = if length(initial_seed) >= 8
        Float64.(initial_seed[1:8])
    else
        extend_seed(Float64.(initial_seed), mode)
    end

    seed_candidates = Vector{Vector{Float64}}()
    push!(seed_candidates, copy(base_seed))
    extra = max(max_seed_candidates - 1, 0)
    for seed in Iterators.take(seed_catalog(mode, [T_fm]), extra)
        s8 = length(seed) >= 8 ? Float64.(seed[1:8]) : extend_seed(Float64.(seed), mode)
        push!(seed_candidates, s8)
    end

    uniq = Dict{String, Vector{Float64}}()
    for s in seed_candidates
        key = join(round.(s; digits=6), ",")
        if !haskey(uniq, key)
            uniq[key] = s
        end
    end
    seed_candidates = collect(values(uniq))

    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num, t_num=t_num, model_kind=model_kind)
    raw_residual! = build_residual!(mode, params)

    best_x = copy(base_seed)
    best_raw = Inf

    for seed in seed_candidates
        x = copy(seed)
        for cfg in _weighted_stage_schedule()
            stage = try
                _run_weighted_stage(raw_residual!, x;
                    w6=cfg.w6,
                    w7=cfg.w7,
                    w8=cfg.w8,
                    method=cfg.method,
                    iterations=cfg.iterations,
                )
            catch
                nothing
            end
            stage === nothing && continue

            x = stage.x
            if isfinite(stage.raw_norm) && stage.raw_norm < best_raw
                best_raw = stage.raw_norm
                best_x = copy(stage.x)
            end

            if isfinite(best_raw) && best_raw <= 1e-8
                early_seed = DefaultSeed(best_x, best_x, :hadron)
                early_result = try
                    solve(mode, T_fm;
                        xi=xi,
                        seed_strategy=early_seed,
                        p_num=p_num,
                        t_num=t_num,
                        model_kind=model_kind,
                        allow_legacy_fallback=true,
                        nlsolve_method=:trust_region,
                        trust_region_fallback=true,
                        residual_norm_max=residual_norm_max,
                        nlsolve_kwargs...)
                catch
                    nothing
                end
                if early_result !== nothing && early_result.converged
                    return early_result
                end
            end
            if best_raw <= 1e-8
                break
            end
        end
        if best_raw <= 1e-8
            break
        end
    end

    final_seed = DefaultSeed(best_x, best_x, :hadron)
    result = try
        solve(mode, T_fm;
            xi=xi,
            seed_strategy=final_seed,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind,
            allow_legacy_fallback=true,
            nlsolve_method=:trust_region,
            trust_region_fallback=true,
            residual_norm_max=residual_norm_max,
            nlsolve_kwargs...)
    catch
        nothing
    end

    return result
end
