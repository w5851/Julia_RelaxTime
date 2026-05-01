"""
Phase D2: 验证稳定粒子极限 meson density 基线的量级与趋势。

当前固定验证窗口：
- xi = 0
- muB = 0
- T = 150:10:190 MeV

验证目标：
1. n_pi > 0, n_K > 0
2. n_pi > n_K
3. n_pi(T), n_K(T), K/pi(T) 随 T 单调上升

输出：
- CSV: data/outputs/results/relaxtime/analysis/meson_density_stable_baseline.csv
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models: solve_gap_and_meson_density_point

const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "analysis", "meson_density_stable_baseline.csv")

@inline function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    return string(x)
end

function main()
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
    mkpath(dirname(output))

    Ts_MeV = collect(150.0:10.0:190.0)
    rows = NamedTuple[]
    continuation_state = nothing

    for T_MeV in Ts_MeV
        T_fm = T_MeV / ħc_MeV_fm
        res = solve_gap_and_meson_density_point(
            T_fm,
            0.0;
            xi=0.0,
            mesons=(:pi, :K),
            continuation_state=continuation_state,
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(; iterations=20),
            mass_kwargs=(; iterations=20),
            density_kwargs=(; num_q_nodes=96),
        )
        md = res.meson_density
        push!(rows, (
            T_MeV=T_MeV,
            m_pi=md.m_pi,
            m_K=md.m_K,
            n_pi=md.n_pi,
            n_K=md.n_K,
            kpi_ratio=md.kpi_ratio,
            positive=(md.n_pi > 0.0 && md.n_K > 0.0),
            pi_gt_k=(md.n_pi > md.n_K),
        ))
        continuation_state = res.continuation_state
    end

    n_pi_mono = all(diff([r.n_pi for r in rows]) .> 0.0)
    n_K_mono = all(diff([r.n_K for r in rows]) .> 0.0)
    ratio_mono = all(diff([r.kpi_ratio for r in rows]) .> 0.0)
    all_positive = all(r.positive for r in rows)
    all_pi_gt_k = all(r.pi_gt_k for r in rows)

    open(output, "w") do io
        println(io, "T_MeV,m_pi,m_K,n_pi,n_K,kpi_ratio,positive,pi_gt_k")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (:T_MeV, :m_pi, :m_K, :n_pi, :n_K, :kpi_ratio, :positive, :pi_gt_k)), ','))
        end
        println(io, "# checks,n_pi_monotonic=$(n_pi_mono),n_K_monotonic=$(n_K_mono),kpi_monotonic=$(ratio_mono),all_positive=$(all_positive),all_pi_gt_k=$(all_pi_gt_k)")
    end

    println("Wrote validation CSV: ", output)
    println("Checks:")
    println("  n_pi monotonic: ", n_pi_mono)
    println("  n_K monotonic: ", n_K_mono)
    println("  K/pi monotonic: ", ratio_mono)
    println("  all positive: ", all_positive)
    println("  all n_pi > n_K: ", all_pi_gt_k)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
