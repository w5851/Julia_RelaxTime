# DualBranchScan.jl 单元测试
#
# 测试内容：
# 1. BranchPoint / DualBranchResult / PhaseTransitionInfo 结构
# 2. run_dual_branch_scan 接口存在
# 3. find_phase_transition / merge_branches 接口存在

using Test
using StaticArrays: SA

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _DBS = Models.DualBranchScan

# ============================================================================

@testset "DualBranchScan" begin

    @testset "BranchPoint 结构" begin
        @test isdefined(_DBS, :BranchPoint)
        fnames = fieldnames(_DBS.BranchPoint)
        @test :mu_mev in fnames
        @test :converged in fnames
        @test :omega in fnames
        @test :pressure in fnames
    end

    @testset "DualBranchResult 结构" begin
        @test isdefined(_DBS, :DualBranchResult)
        fnames = fieldnames(_DBS.DualBranchResult)
        @test :T_mev in fnames
        @test :hadron_branch in fnames
        @test :quark_branch in fnames
    end

    @testset "PhaseTransitionInfo 结构" begin
        @test isdefined(_DBS, :PhaseTransitionInfo)
        fnames = fieldnames(_DBS.PhaseTransitionInfo)
        @test :found in fnames
        @test :mu_c in fnames
    end

    @testset "run_dual_branch_scan 接口存在" begin
        @test isdefined(_DBS, :run_dual_branch_scan)
        @test _DBS.run_dual_branch_scan isa Function
    end

    @testset "find_phase_transition 接口存在" begin
        @test isdefined(_DBS, :find_phase_transition)
        @test _DBS.find_phase_transition isa Function
    end

    @testset "merge_branches 接口存在" begin
        @test isdefined(_DBS, :merge_branches)
        @test _DBS.merge_branches isa Function
    end

    @testset "auto backend 路由规则" begin
        @test _DBS._effective_solver_backend(:models, :auto) == :models
        @test _DBS._effective_solver_backend(:legacy, :auto; auto_pnjl_backend=:models) == :models
        @test _DBS._effective_solver_backend(:legacy, :auto; auto_pnjl_backend=:legacy) == :legacy
        @test _DBS._effective_solver_backend(:legacy, :models; auto_pnjl_backend=:legacy) == :models
    end

    @testset "semantic_mode selector 解析" begin
        h = _DBS.BranchPoint(100.0, true, -10.0, 10.0, 0.2, 0.3, 0.4, SA[0.1, 0.1, 0.1, 0.2, 0.2], SA[0.3, 0.3, 0.5], 10, 1e-9)
        q = _DBS.BranchPoint(100.0, true, -9.0, 9.0, 0.25, 0.35, 0.45, SA[0.12, 0.1, 0.1, 0.25, 0.25], SA[0.28, 0.28, 0.48], 12, 2e-9)

        default_ground = _DBS._resolve_branch_selector(:ground_state, nothing)
        @test default_ground(h, q) == :hadron

        default_manifold = _DBS._resolve_branch_selector(:constrained_manifold, nothing)
        @test default_manifold(h, q) == :hadron

        custom_selector = (hh, qq) -> hh.pressure > qq.pressure ? :hadron : :quark
        chosen = _DBS._resolve_branch_selector(:ground_state, custom_selector)
        @test chosen(h, q) == :hadron

        @test_throws ArgumentError _DBS._resolve_branch_selector(:unknown_mode, nothing)
    end

    @testset "physical branch supports custom selector" begin
        h = _DBS.BranchPoint(120.0, true, -8.0, 8.0, 0.2, 0.3, 0.4, SA[0.1, 0.1, 0.1, 0.2, 0.2], SA[0.3, 0.3, 0.5], 8, 1e-9)
        q = _DBS.BranchPoint(120.0, true, -9.0, 9.0, 0.25, 0.35, 0.45, SA[0.5, 0.5, 0.5, 0.8, 0.8], SA[0.2, 0.2, 0.4], 9, 2e-9)

        hadron = Union{Nothing, _DBS.BranchPoint}[h]
        quark = Union{Nothing, _DBS.BranchPoint}[q]

        by_ground = _DBS._select_physical_branch(hadron, quark, _DBS._resolve_branch_selector(:ground_state, nothing))
        @test by_ground[1] === q

        by_manifold = _DBS._select_physical_branch(hadron, quark, _DBS._resolve_branch_selector(:constrained_manifold, nothing))
        @test by_manifold[1] === h

        invalid_selector = (_h, _q) -> :invalid
        @test_throws ArgumentError _DBS._select_physical_branch(hadron, quark, invalid_selector)
    end

    @testset "models 路径禁止 legacy solver 开关" begin
        @test _DBS._reject_legacy_solver_kwargs((; solver=:newton)) === nothing
        @test_throws ArgumentError _DBS._reject_legacy_solver_kwargs((; use_problem_spec=false))
        @test_throws ArgumentError _DBS._reject_legacy_solver_kwargs((; allow_legacy_path=true))
        @test_throws ArgumentError _DBS._reject_legacy_solver_kwargs((; warn_on_legacy_path=false))
        @test_throws ArgumentError _DBS._reject_legacy_solver_kwargs((; problem_spec=:dummy))
    end
end
