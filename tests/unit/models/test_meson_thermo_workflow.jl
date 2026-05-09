using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _MTW = Models.MesonThermoWorkflow

@testset "MesonThermoWorkflow" begin
    @testset "接口存在" begin
        @test isdefined(_MTW, :solve_meson_thermo_from_meson_point)
        @test isdefined(_MTW, :solve_gap_and_meson_thermo_point)
        @test isdefined(_MTW, :solve_strict_bw_meson_thermo_from_meson_point)
        @test isdefined(_MTW, :solve_gap_and_strict_bw_meson_thermo_point)
        @test isdefined(_MTW, :solve_phase_shift_meson_thermo_from_meson_point)
        @test isdefined(_MTW, :solve_gap_and_phase_shift_meson_thermo_point)
        @test isdefined(_MTW, :build_meson_thermo_contract_row)
    end

    @testset "stable meson thermo workflow 生成总热力学字段" begin
        T_fm = 170.0 / Main.Constants_PNJL.ħc_MeV_fm
        point = Models.solve_gap_and_meson_thermo_point(
            T_fm,
            0.0;
            xi=0.0,
            mesons=(:pi, :K),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
            thermo_kwargs=(; num_q_nodes=64, p_num=8, t_num=4),
            temperature_step_fm=0.5 / Main.Constants_PNJL.ħc_MeV_fm,
        )

        thermo = point.meson_thermo
        row = Models.build_meson_thermo_contract_row(point)

        @test thermo.workflow == :stable_meson_pressure
        @test thermo.P_meson > 0.0
        @test thermo.P_total > thermo.P_quark_meanfield
        @test isfinite(thermo.entropy)
        @test isfinite(thermo.epsilon)
        @test isfinite(thermo.trace_anomaly)
        @test row.workflow == "stable_meson_pressure"
        @test row.channel_set == "pi,K"
        @test row.P_meson ≈ thermo.P_meson rtol=1e-12
    end

    @testset "phase-shift meson thermo workflow 支持 current / gbu 两口径" begin
        T_fm = 210.0 / Main.Constants_PNJL.ħc_MeV_fm
        meson_point = Models.solve_gap_and_meson_point(
            T_fm,
            0.0;
            xi=0.0,
            mesons=(:pi, :K),
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=20,),
            mass_kwargs=(iterations=20,),
        )

        current = Models.solve_phase_shift_meson_thermo_from_meson_point(
            meson_point;
            scheme=:current,
            qmax=4.0,
            q_nodes=6,
            omega_max=3.0,
            omega_nodes=6,
            p_num=8,
            t_num=4,
        )
        gbu = Models.solve_phase_shift_meson_thermo_from_meson_point(
            meson_point;
            scheme=:gbu_reference,
            qmax=4.0,
            q_nodes=6,
            omega_max=3.0,
            omega_nodes=6,
            p_num=8,
            t_num=4,
        )

        @test current.workflow == :phase_shift_current
        @test gbu.workflow == :phase_shift_gbu_reference
        @test current.phase_shift_variant == :phase_shift_current
        @test gbu.phase_shift_variant == :phase_shift_gbu_reference
        @test isfinite(current.P_meson)
        @test isfinite(gbu.P_meson)
        @test current.P_meson > 0.0
        @test gbu.P_meson > 0.0
        @test current.P_meson ≈ current.P_meson_qp + current.P_meson_ld rtol=1e-12
        @test gbu.P_meson ≈ gbu.P_meson_qp + gbu.P_meson_ld rtol=1e-12
        @test current.ld_cutoff ≈ 4.0 rtol=1e-12
        @test current.ld_cutoff_mode == :match_qmax
        @test current.ld_threshold_mode == :omega_lt_q
        @test !isapprox(current.P_meson, gbu.P_meson; rtol=1e-6, atol=1e-10)

        row = Models.build_meson_thermo_contract_row((; meson_point..., phase_shift_meson_thermo=current))
        @test row.P_meson_qp ≈ current.P_meson_qp rtol=1e-12
        @test row.P_meson_ld ≈ current.P_meson_ld rtol=1e-12
        @test row.P_meson ≈ row.P_meson_qp + row.P_meson_ld rtol=1e-12
        @test row.ld_cutoff ≈ current.ld_cutoff rtol=1e-12
        @test row.ld_cutoff_mode == "match_qmax"
        @test row.ld_threshold_mode == "omega_lt_q"
    end
end
