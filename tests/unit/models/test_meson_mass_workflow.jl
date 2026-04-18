# MesonMassWorkflow.jl 单元测试
#
# 测试内容：
# 1. 模块加载
# 2. DEFAULT_MESONS 常量
# 3. solve_gap_and_meson_point / build_equilibrium_params 接口存在

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _MW = Models.MesonMassWorkflow

# ============================================================================

@testset "MesonMassWorkflow" begin

    @testset "DEFAULT_MESONS 常量" begin
        @test isdefined(_MW, :DEFAULT_MESONS)
        mesons = _MW.DEFAULT_MESONS
        @test mesons isa Tuple
        @test :pi in mesons
        @test :K in mesons
        @test :sigma_pi in mesons
        @test length(mesons) >= 6
    end

    @testset "solve_gap_and_meson_point 接口存在" begin
        @test isdefined(_MW, :solve_gap_and_meson_point)
        @test _MW.solve_gap_and_meson_point isa Function
    end

    @testset "build_equilibrium_params 接口存在" begin
        @test isdefined(_MW, :build_equilibrium_params)
        @test _MW.build_equilibrium_params isa Function
    end

    @testset "meson root diagnostics expose governance fields" begin
        T_fm = 240.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0
        res = _MW.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi,),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )
        diag = res.meson_results[:pi].root_diagnostics
        @test hasproperty(diag, :governance_candidate_count)
        @test hasproperty(diag, :governance_selection_reason)
        @test diag.governance_candidate_count >= 0
    end

    @testset "mixed diagnostics expose continuity and second-pass fields" begin
        T_fm = 250.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0
        res = _MW.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=-0.3,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
        )
        diag = res.meson_results[:eta_prime].root_diagnostics
        @test hasproperty(diag, :continuity_penalty_weight)
        @test hasproperty(diag, :continuity_penalty_distance)
        @test hasproperty(diag, :second_pass_triggered)
        @test hasproperty(diag, :second_pass_candidate_count)
        @test hasproperty(diag, :second_pass_eta_prime_fine_gamma)
        @test hasproperty(diag, :continuity_guard_applied)
        @test hasproperty(diag, :nudged_restart_applied)
        @test hasproperty(diag, :roi_rescue_attempted)
        @test hasproperty(diag, :roi_rescue_applied)
        @test hasproperty(diag, :roi_rescue_candidate_count)
        @test hasproperty(diag, :pre_joint_mass)
        @test hasproperty(diag, :pre_joint_gamma)
        @test hasproperty(diag, :pre_joint_residual)
        @test hasproperty(diag, :post_joint_mass)
        @test hasproperty(diag, :post_joint_gamma)
        @test hasproperty(diag, :post_joint_residual)
    end

    @testset "eta_prime ROI rescue diagnostics fields exist in hard edge window" begin
        T_fm = 255.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0
        res = _MW.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=-0.3,
            mesons=(:eta_prime,),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
        )
        diag = res.meson_results[:eta_prime].root_diagnostics
        @test hasproperty(diag, :roi_rescue_attempted)
        @test hasproperty(diag, :roi_rescue_candidate_count)
    end

    @testset "joint pair objective can keep T255 eta_prime in ROI under continuation" begin
        muq_fm = 0.0
        meson_seed_state = nothing
        mixed_seed_tracking_state = nothing
        etap = nothing

        for T_MeV in (245.0, 250.0, 255.0)
            T_fm = T_MeV / Main.Constants_PNJL.ħc_MeV_fm
            res = _MW.solve_gap_and_meson_point(
                T_fm,
                muq_fm;
                xi=-0.3,
                mesons=(:eta, :eta_prime),
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=25,),
                mass_kwargs=(iterations=25,),
                meson_seed_state=meson_seed_state,
                mixed_seed_tracking_state=mixed_seed_tracking_state,
            )
            etap = res.meson_results[:eta_prime]
            meson_seed_state = res.meson_seed_state
            mixed_seed_tracking_state = res.mixed_seed_tracking
        end

        @test etap !== nothing
        @test 2.80 <= etap.mass <= 2.95
        @test 2.20 <= etap.gamma <= 2.60
        @test etap.residual <= 0.03
    end

    @testset "global-fallback diagnostics fields exist" begin
        T_fm = 255.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0
        res = _MW.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=-0.3,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
        )
        diag = res.meson_results[:eta_prime].root_diagnostics
        @test hasproperty(diag, :global_fallback_attempted)
        @test hasproperty(diag, :global_fallback_applied)
        @test hasproperty(diag, :global_fallback_candidate_count)
        @test hasproperty(diag, :b2_guard_applied)
        if diag.global_fallback_attempted
            @test diag.global_fallback_candidate_count <= 4
        end
    end

    @testset "force_global_fallback diagnostics fields exist" begin
        T_fm = 200.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0
        res = _MW.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=-0.3,
            mesons=(:eta_prime,),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
            force_global_fallback=true,
        )
        diag = res.meson_results[:eta_prime].root_diagnostics
        @test hasproperty(diag, :global_fallback_attempted)
        @test hasproperty(diag, :global_fallback_applied)
        @test hasproperty(diag, :global_fallback_candidate_count)
    end

    @testset "B2 guard reduces T260 gamma back-jump under forced fallback" begin
        muq_fm = 0.0
        meson_seed_state = nothing
        mixed_seed_tracking_state = nothing
        g255 = NaN
        g260 = NaN

        for T_MeV in (245.0, 250.0, 255.0, 260.0)
            T_fm = T_MeV / Main.Constants_PNJL.ħc_MeV_fm
            res = _MW.solve_gap_and_meson_point(
                T_fm,
                muq_fm;
                xi=-0.3,
                mesons=(:eta, :eta_prime),
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=25,),
                mass_kwargs=(iterations=25,),
                meson_seed_state=meson_seed_state,
                mixed_seed_tracking_state=mixed_seed_tracking_state,
                force_global_fallback=true,
            )
            etap = res.meson_results[:eta_prime]
            if T_MeV == 255.0
                g255 = etap.gamma
            elseif T_MeV == 260.0
                g260 = etap.gamma
            end
            meson_seed_state = res.meson_seed_state
            mixed_seed_tracking_state = res.mixed_seed_tracking
        end

        @test isfinite(g255)
        @test isfinite(g260)
        @test g255 >= 2.3
        @test g260 >= 2.3
    end

    @testset "blackbox fallback helper exists" begin
        @test isdefined(_MW, :_blackbox_global_fallback_seed)
    end
end
