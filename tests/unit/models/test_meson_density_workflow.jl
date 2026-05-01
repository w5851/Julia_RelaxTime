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
        @test _MDW.solve_meson_density_from_meson_point isa Function
        @test _MDW.solve_gap_and_meson_density_point isa Function
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
end
