using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _MDW = Models.MesonDensityWorkflow

@testset "MesonDensityWorkflow" begin
    @testset "接口存在" begin
        @test isdefined(_MDW, :solve_meson_density_from_meson_point)
        @test isdefined(_MDW, :solve_gap_and_meson_density_point)
        @test isdefined(_MDW, :solve_strict_bw_meson_density_from_meson_point)
        @test isdefined(_MDW, :solve_gap_and_strict_bw_meson_density_point)
        @test isdefined(_MDW, :solve_phase_shift_meson_density_from_meson_point)
        @test isdefined(_MDW, :solve_gap_and_phase_shift_meson_density_point)
        @test _MDW.solve_meson_density_from_meson_point isa Function
        @test _MDW.solve_gap_and_meson_density_point isa Function
        @test _MDW.solve_strict_bw_meson_density_from_meson_point isa Function
        @test _MDW.solve_gap_and_strict_bw_meson_density_point isa Function
        @test _MDW.solve_phase_shift_meson_density_from_meson_point isa Function
        @test _MDW.solve_gap_and_phase_shift_meson_density_point isa Function
    end

    @testset "后处理消费 meson workflow 返回值" begin
        T_fm = 170.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0
        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )

        density = Models.solve_meson_density_from_meson_point(
            meson_point;
            num_q_nodes=128,
        )

        @test density.T_fm ≈ meson_point.thermo_params.T
        @test density.m_pi ≈ meson_point.meson_results[:pi].mass
        @test density.m_K ≈ meson_point.meson_results[:K].mass
        @test density.n_pi > 0.0
        @test density.n_K > 0.0
        @test 0.0 < density.kpi_ratio < 1.0
    end

    @testset "strict BW 入口复用 meson workflow 主链" begin
        T_fm = 210.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0

        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )

        density = Models.solve_strict_bw_meson_density_from_meson_point(
            meson_point;
            qmax=12.0,
            q_nodes=16,
            omega_max=10.0,
            omega_nodes=16,
        )

        full = Models.solve_gap_and_strict_bw_meson_density_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
            density_kwargs=(; qmax=12.0, q_nodes=16, omega_max=10.0, omega_nodes=16),
        )

        @test density.m_pi ≈ meson_point.meson_results[:pi].mass
        @test density.m_K ≈ meson_point.meson_results[:K].mass
        @test density.gamma_pi ≈ meson_point.meson_results[:pi].gamma
        @test density.gamma_K ≈ meson_point.meson_results[:K].gamma
        @test isfinite(density.n_pi)
        @test isfinite(density.n_K)
        @test density.n_pi > 0.0
        @test density.n_K > 0.0
        @test 0.0 < density.kpi_ratio < 2.0

        @test hasproperty(full, :strict_bw_meson_density)
        @test full.strict_bw_meson_density.n_pi ≈ density.n_pi rtol=1e-10
        @test full.strict_bw_meson_density.n_K ≈ density.n_K rtol=1e-10
        @test full.strict_bw_meson_density.gamma_pi ≈ density.gamma_pi rtol=1e-10
        @test full.strict_bw_meson_density.gamma_K ≈ density.gamma_K rtol=1e-10
    end

    @testset "strict BW Stage2 q-pole 入口可运行" begin
        T_fm = 210.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0

        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )

        density = Models.solve_strict_bw_meson_density_from_meson_point(
            meson_point;
            stage=:stage2_qpole,
            qmax=12.0,
            q_nodes=4,
            omega_max=10.0,
            omega_nodes=8,
            solver_iterations=12,
            pole_residual_norm_max=1e-4,
            pole_require_converged=false,
        )

        @test density.stage == :strict_bw_stage2_qpole
        @test isfinite(density.n_pi)
        @test isfinite(density.n_K)
        @test density.n_pi >= 0.0
        @test density.n_K >= 0.0
        @test length(density.pi_density.q_values) == 4
        @test length(density.k_density.q_values) == 4
        @test all(isfinite, density.pi_density.residual_norms)
        @test all(isfinite, density.k_density.residual_norms)
    end

    @testset "完整入口复用 meson workflow 主链" begin
        T_fm = 170.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0

        full = Models.solve_gap_and_meson_density_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
            density_kwargs=(; num_q_nodes=128),
        )

        post = Models.solve_meson_density_from_meson_point(full; num_q_nodes=128)

        @test hasproperty(full, :meson_density)
        @test full.meson_density.m_pi ≈ full.meson_results[:pi].mass
        @test full.meson_density.m_K ≈ full.meson_results[:K].mass
        @test full.meson_density.n_pi ≈ post.n_pi rtol=1e-10
        @test full.meson_density.n_K ≈ post.n_K rtol=1e-10
        @test full.meson_density.kpi_ratio ≈ post.kpi_ratio rtol=1e-10
    end

    @testset "扫描 continuation 下与 meson workflow 质量结果一致" begin
        continuation_meson = nothing
        continuation_density = nothing

        for T_MeV in (150.0, 170.0, 190.0)
            T_fm = T_MeV / Main.Constants_PNJL.ħc_MeV_fm

            meson_only = Models.solve_gap_and_meson_point(
                T_fm,
                0.0;
                xi=0.0,
                mesons=(:pi, :K),
                continuation_state=continuation_meson,
                mixed_branch_align=:strict_sign_binding,
                p_num=8,
                t_num=4,
                solver_kwargs=(; iterations=20),
                mass_kwargs=(; iterations=20),
            )

            with_density = Models.solve_gap_and_meson_density_point(
                T_fm,
                0.0;
                xi=0.0,
                mesons=(:pi, :K),
                continuation_state=continuation_density,
                mixed_branch_align=:strict_sign_binding,
                p_num=8,
                t_num=4,
                solver_kwargs=(; iterations=20),
                mass_kwargs=(; iterations=20),
                density_kwargs=(; num_q_nodes=64),
            )

            @test with_density.meson_results[:pi].mass ≈ meson_only.meson_results[:pi].mass rtol=1e-10
            @test with_density.meson_results[:K].mass ≈ meson_only.meson_results[:K].mass rtol=1e-10
            @test with_density.meson_density.m_pi ≈ meson_only.meson_results[:pi].mass rtol=1e-10
            @test with_density.meson_density.m_K ≈ meson_only.meson_results[:K].mass rtol=1e-10

            continuation_meson = meson_only.continuation_state
            continuation_density = with_density.continuation_state
        end
    end

    @testset "Phase-E3 相移数密度入口复用 meson workflow 主链" begin
        T_fm = 210.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0

        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )

        density = Models.solve_phase_shift_meson_density_from_meson_point(
            meson_point;
            qmax=12.0,
            q_nodes=16,
            omega_max=10.0,
            omega_nodes=16,
        )

        full = Models.solve_gap_and_phase_shift_meson_density_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
            density_kwargs=(; qmax=12.0, q_nodes=16, omega_max=10.0, omega_nodes=16),
        )

        @test density.m_pi ≈ meson_point.meson_results[:pi].mass
        @test density.m_K ≈ meson_point.meson_results[:K].mass
        @test isfinite(density.n_pi)
        @test isfinite(density.n_K)
        @test density.n_pi > 0.0
        @test density.n_K > 0.0
        @test 0.0 < density.kpi_ratio < 2.0
        @test density.qmax == 12.0
        @test density.q_nodes == 16
        @test density.omega_max == 10.0
        @test density.omega_nodes == 16

        @test hasproperty(full, :phase_shift_meson_density)
        @test full.phase_shift_meson_density.n_pi ≈ density.n_pi rtol=1e-10
        @test full.phase_shift_meson_density.n_K ≈ density.n_K rtol=1e-10
    end

    @testset "Phase-E3 generalized BU reference 复用同一 workflow 入口" begin
        T_fm = 210.0 / Main.Constants_PNJL.ħc_MeV_fm
        muq_fm = 0.0

        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )

        current = Models.solve_phase_shift_meson_density_from_meson_point(
            meson_point;
            scheme=:current,
            qmax=12.0,
            q_nodes=12,
            omega_max=10.0,
            omega_nodes=12,
        )
        gbu = Models.solve_phase_shift_meson_density_from_meson_point(
            meson_point;
            scheme=:gbu_reference,
            qmax=12.0,
            q_nodes=12,
            omega_max=10.0,
            omega_nodes=12,
        )

        full = Models.solve_gap_and_phase_shift_meson_density_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
            density_kwargs=(; scheme=:gbu_reference, qmax=12.0, q_nodes=12, omega_max=10.0, omega_nodes=12),
        )

        @test current.scheme == :phase_shift_current
        @test gbu.scheme == :phase_shift_gbu_reference
        @test hasproperty(full, :phase_shift_meson_density)
        @test full.phase_shift_meson_density.scheme == :phase_shift_gbu_reference
        @test full.phase_shift_meson_density.pi_density.scheme == :phase_shift_gbu_reference
        @test full.phase_shift_meson_density.k_density.scheme == :phase_shift_gbu_reference
        @test isfinite(gbu.n_pi)
        @test isfinite(gbu.n_K)
        @test gbu.n_pi > 0.0
        @test gbu.n_K > 0.0
        @test !isapprox(current.n_pi, gbu.n_pi; rtol=1e-6, atol=1e-10)
        @test !isapprox(current.n_K, gbu.n_K; rtol=1e-6, atol=1e-10)
        @test full.phase_shift_meson_density.n_pi ≈ gbu.n_pi rtol=1e-10
        @test full.phase_shift_meson_density.n_K ≈ gbu.n_K rtol=1e-10
    end
end
