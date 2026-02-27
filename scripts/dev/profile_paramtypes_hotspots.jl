"""
Parameter Types 热点 profiling（Phase D）

用途：
- 评估参数归一化/转换开销在代表性调用路径中的占比
- 给“是否值得做结构体直通特化”提供量化依据

运行：
  julia --project=. scripts/dev/profile_paramtypes_hotspots.jl
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "MesonMassWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TransportCoefficients.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxationTime.jl"))

using .TransportWorkflow
using .MesonMassWorkflow
using .TransportCoefficients
using .RelaxationTime
using Main.ParameterTypes: as_namedtuple

@inline function avg_elapsed_us(f, n::Int)
    t = @elapsed for _ in 1:n
        f()
    end
    return t / n * 1e6
end

function main()
    T = 0.15
    mu = 0.0
    xi = 0.0

    println("[Phase F] preparing equilibrium baseline...")
    base = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
        T,
        mu;
        xi=xi,
        thermo_backend=:models,
        solver_backend=:models,
        p_num=8,
        t_num=4,
        seed_state=TransportWorkflow.PNJL.HADRON_SEED_5,
        models_residual_norm_max=1e-4,
    )

    inputs = TransportWorkflow._transport_inputs_from_equilibrium(base, T, mu;
        xi=xi,
        thermo_backend=:legacy,
        p_num=8,
        t_num=4,
    )

    q_struct = inputs.quark_params
    t_struct = inputs.thermo_params
    q_nt = as_namedtuple(q_struct)
    t_nt = as_namedtuple(t_struct)

    println("[Phase F] warm-up...")
    TransportWorkflow.normalize_quark_params(q_struct)
    TransportWorkflow.normalize_thermo_params(t_struct)
    TransportWorkflow.normalize_quark_params(q_nt)
    TransportWorkflow.normalize_thermo_params(t_nt)
    TransportWorkflow._A_from_equilibrium(T, q_struct, t_struct;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )
    TransportWorkflow._A_from_equilibrium(T, q_nt, t_nt;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    n_fast = 400_000
    n_mid = 300
    n_solve = 80
    n_tc = 200
    n_rt = 80_000

    tau_one = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    tc_cfg = TransportCoefficients.TransportIntegrationConfig(p_nodes=8, p_max=8.0, cos_nodes=4)
    tc_req = TransportCoefficients.TransportRequest(q_struct, t_struct; tau=tau_one, integration=tc_cfg)

    TransportCoefficients.transport_coefficients(tc_req)
    TransportCoefficients.transport_coefficients(
        (m=tc_req.quark.m, μ=tc_req.quark.μ),
        (T=tc_req.thermo.T, Φ=tc_req.thermo.Φ, Φbar=tc_req.thermo.Φbar, ξ=tc_req.thermo.ξ);
        tau=tc_req.tau,
        config=tc_req.integration,
    )

    rt_rates_sample = (
        uu_to_uu=1.0,
        ud_to_ud=2.0,
        us_to_us=3.0,
        udbar_to_udbar=4.0,
        dubar_to_dubar=5.0,
        uubar_to_uubar=6.0,
        uubar_to_ddbar=7.0,
        usbar_to_usbar=8.0,
        subar_to_subar=9.0,
        uubar_to_ssbar=10.0,
        ss_to_ss=11.0,
        ssbar_to_uubar=12.0,
        ssbar_to_ssbar=13.0,
        ubardbar_to_ubardbar=14.0,
        ubarubar_to_ubarubar=15.0,
        ubarsbar_to_ubarsbar=16.0,
        sbarsbar_to_sbarsbar=17.0,
    )
    rt_densities = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    rt_k = (K_σπ=1.0, K_σK=1.0, K_σ=1.0, K_δπ=1.0, K_δK=1.0)

    RelaxationTime.relaxation_times(q_struct, t_struct, rt_k; densities=rt_densities, existing_rates=rt_rates_sample)
    RelaxationTime.relaxation_times(q_nt, t_nt, rt_k; densities=rt_densities, existing_rates=rt_rates_sample)


    println("[Phase F] measuring...")
    norm_struct_acc = Ref(0.0)
    norm_struct_us = avg_elapsed_us(() -> begin
        qn = TransportWorkflow.normalize_quark_params(q_struct)
        tn = TransportWorkflow.normalize_thermo_params(t_struct)
        norm_struct_acc[] += qn.m.u + tn.T
        nothing
    end, n_fast)

    norm_nt_acc = Ref(0.0)
    norm_nt_us = avg_elapsed_us(() -> begin
        qn = TransportWorkflow.normalize_quark_params(q_nt)
        tn = TransportWorkflow.normalize_thermo_params(t_nt)
        norm_nt_acc[] += qn.m.u + tn.T
        nothing
    end, n_fast)

    A_struct_acc = Ref(0.0)
    A_struct_us = avg_elapsed_us(() -> begin
        A = TransportWorkflow._A_from_equilibrium(T, q_struct, t_struct;
            a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
        )
        A_struct_acc[] += A.u
        nothing
    end, n_mid)

    A_nt_acc = Ref(0.0)
    A_nt_us = avg_elapsed_us(() -> begin
        A = TransportWorkflow._A_from_equilibrium(T, q_nt, t_nt;
            a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
        )
        A_nt_acc[] += A.u
        nothing
    end, n_mid)

    solve_acc = Ref(0.0)
    solve_us = avg_elapsed_us(() -> begin
        s = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            solver_backend=:models,
            p_num=8,
            t_num=4,
            seed_state=TransportWorkflow.PNJL.HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )
        solve_acc[] += s.x_state[1]
        nothing
    end, n_solve)

    tc_req_acc = Ref(0.0)
    tc_req_us = avg_elapsed_us(() -> begin
        out = TransportCoefficients.transport_coefficients(tc_req)
        tc_req_acc[] += out.eta + out.sigma
        nothing
    end, n_tc)

    tc_manual_nt_acc = Ref(0.0)
    tc_manual_nt_us = avg_elapsed_us(() -> begin
        out = TransportCoefficients.transport_coefficients(
            (m=tc_req.quark.m, μ=tc_req.quark.μ),
            (T=tc_req.thermo.T, Φ=tc_req.thermo.Φ, Φbar=tc_req.thermo.Φbar, ξ=tc_req.thermo.ξ);
            tau=tc_req.tau,
            config=tc_req.integration,
        )
        tc_manual_nt_acc[] += out.eta + out.sigma
        nothing
    end, n_tc)

    rt_struct_acc = Ref(0.0)
    rt_struct_us = avg_elapsed_us(() -> begin
        out = RelaxationTime.relaxation_times(
            q_struct,
            t_struct,
            rt_k;
            densities=rt_densities,
            existing_rates=rt_rates_sample,
        )
        rt_struct_acc[] += out.tau.u
        nothing
    end, n_rt)

    rt_nt_acc = Ref(0.0)
    rt_nt_us = avg_elapsed_us(() -> begin
        out = RelaxationTime.relaxation_times(
            q_nt,
            t_nt,
            rt_k;
            densities=rt_densities,
            existing_rates=rt_rates_sample,
        )
        rt_nt_acc[] += out.tau.u
        nothing
    end, n_rt)

    share_norm_struct_vs_A = norm_struct_us / A_struct_us
    share_norm_nt_vs_A = norm_nt_us / A_nt_us
    share_norm_struct_vs_solve = norm_struct_us / solve_us
    share_norm_nt_vs_solve = norm_nt_us / solve_us

    println("\n=== ParameterTypes Hotspot Baseline (us/call) ===")
    println("normalize(struct inputs): ", round(norm_struct_us; digits=4))
    println("normalize(nt inputs):     ", round(norm_nt_us; digits=4))
    println("A_from_equilibrium(struct): ", round(A_struct_us; digits=4))
    println("A_from_equilibrium(nt):     ", round(A_nt_us; digits=4))
    println("equilibrium solve baseline: ", round(solve_us; digits=4))
    println("transport_coeff(req struct): ", round(tc_req_us; digits=4))
    println("transport_coeff(explicit nt): ", round(tc_manual_nt_us; digits=4))
    println("relaxation_times(struct, existing): ", round(rt_struct_us; digits=4))
    println("relaxation_times(nt, existing):     ", round(rt_nt_us; digits=4))

    println("\n=== Relative Share ===")
    println("normalize(struct) / A(struct): ", round(share_norm_struct_vs_A; digits=4))
    println("normalize(nt) / A(nt):         ", round(share_norm_nt_vs_A; digits=4))
    println("normalize(struct) / solve:     ", round(share_norm_struct_vs_solve; digits=4))
    println("normalize(nt) / solve:         ", round(share_norm_nt_vs_solve; digits=4))
    println("transport(req)/transport(explicit nt): ", round(tc_req_us / tc_manual_nt_us; digits=4))
    println("relaxation(struct)/relaxation(nt): ", round(rt_struct_us / rt_nt_us; digits=4))

    println("\n=== Suggested Next Targets ===")
    println("1) AverageScatteringRate: continue reducing nested normalization in hot loops")
    println("2) AverageScatteringRate/TotalCrossSection: evaluate coalescing opportunities with model-ready fixtures")
end

main()
