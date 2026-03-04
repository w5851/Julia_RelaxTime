# DualBranchScan.jl 单元测试
#
# 测试内容：
# 1. BranchPoint / DualBranchResult / PhaseTransitionInfo 结构
# 2. run_dual_branch_scan 接口存在
# 3. find_phase_transition / merge_branches 接口存在

using Test

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
end
