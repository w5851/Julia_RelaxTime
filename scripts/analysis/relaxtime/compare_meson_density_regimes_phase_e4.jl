"""
Phase E4: 在同一 meson workflow / continuation 口径下，并排比较
stable / BW / Phase-E3 三条介子数密度链。

目标：
- 固定同一组 `meson_point`
- 对比稳定粒子极限、BW 过渡代理与 Phase-E3 最小相移双积分
- 输出正式比较 CSV，供 Phase E4 收口说明使用
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(@__DIR__, "validate_meson_density_bw_minimal.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models: solve_gap_and_meson_point,
               solve_meson_density_from_meson_point,
               solve_phase_shift_meson_density_from_meson_point

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "meson_density_phase_e4_compare.csv",
)

const DEFAULT_TS_MEV = collect(208.0:2.0:220.0)

@inline function _fmt(x)
    x isa Bool && return (x ? "true" : "false")
    return string(x)
end

function _parse_float_list(raw::AbstractString)
    vals = Float64[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("temperature list cannot be empty"))
    return vals
end

function _selected_temperatures()
    raw = strip(get(ENV, "MESON_E4_T_VALUES", ""))
    isempty(raw) && return DEFAULT_TS_MEV
    return _parse_float_list(raw)
end

function main()
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
    mkpath(dirname(output))

    Ts_MeV = _selected_temperatures()
    continuation_state = nothing
    rows = NamedTuple[]

    for T_MeV in Ts_MeV
        T_fm = T_MeV / ħc_MeV_fm
        meson_point = solve_gap_and_meson_point(
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
        )

        stable = solve_meson_density_from_meson_point(meson_point; num_q_nodes=256)
        phase = solve_phase_shift_meson_density_from_meson_point(
            meson_point;
            qmax=12.0,
            q_nodes=48,
            omega_min=0.05,
            omega_max=10.0,
            omega_nodes=48,
            eta=1e-6,
        )

        mpi = meson_point.meson_results[:pi]
        mK = meson_point.meson_results[:K]
        n_pi_bw = bw_meson_number_density(Float64(mpi.mass), Float64(mpi.gamma), T_fm; degeneracy=stable.d_pi)
        n_K_bw = bw_meson_number_density(Float64(mK.mass), Float64(mK.gamma), T_fm; degeneracy=stable.d_K)
        kpi_bw = iszero(n_pi_bw) ? NaN : n_K_bw / n_pi_bw

        push!(rows, (
            T_MeV=T_MeV,
            m_pi=stable.m_pi,
            gamma_pi=Float64(mpi.gamma),
            m_K=stable.m_K,
            gamma_K=Float64(mK.gamma),
            n_pi_stable=stable.n_pi,
            n_pi_bw=n_pi_bw,
            n_pi_phase=phase.n_pi,
            bw_over_stable_pi=iszero(stable.n_pi) ? NaN : n_pi_bw / stable.n_pi,
            phase_over_stable_pi=iszero(stable.n_pi) ? NaN : phase.n_pi / stable.n_pi,
            n_K_stable=stable.n_K,
            n_K_bw=n_K_bw,
            n_K_phase=phase.n_K,
            bw_over_stable_K=iszero(stable.n_K) ? NaN : n_K_bw / stable.n_K,
            phase_over_stable_K=iszero(stable.n_K) ? NaN : phase.n_K / stable.n_K,
            kpi_stable=stable.kpi_ratio,
            kpi_bw=kpi_bw,
            kpi_phase=phase.kpi_ratio,
            delta_kpi_bw=kpi_bw - stable.kpi_ratio,
            delta_kpi_phase=phase.kpi_ratio - stable.kpi_ratio,
            pi_shell_phase=phase.pi_density.omega_shell_at_qmax,
            K_shell_phase=phase.k_density.omega_shell_at_qmax,
        ))

        continuation_state = meson_point.continuation_state
    end

    open(output, "w") do io
        println(io, "T_MeV,m_pi,gamma_pi,m_K,gamma_K,n_pi_stable,n_pi_bw,n_pi_phase,bw_over_stable_pi,phase_over_stable_pi,n_K_stable,n_K_bw,n_K_phase,bw_over_stable_K,phase_over_stable_K,kpi_stable,kpi_bw,kpi_phase,delta_kpi_bw,delta_kpi_phase,pi_shell_phase,K_shell_phase")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :m_pi, :gamma_pi, :m_K, :gamma_K,
                :n_pi_stable, :n_pi_bw, :n_pi_phase,
                :bw_over_stable_pi, :phase_over_stable_pi,
                :n_K_stable, :n_K_bw, :n_K_phase,
                :bw_over_stable_K, :phase_over_stable_K,
                :kpi_stable, :kpi_bw, :kpi_phase,
                :delta_kpi_bw, :delta_kpi_phase,
                :pi_shell_phase, :K_shell_phase
            )), ','))
        end
    end

    println("Wrote Phase E4 comparison CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
