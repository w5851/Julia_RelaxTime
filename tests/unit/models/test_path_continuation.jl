using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

if !isdefined(Main, :PathContinuationDummyModel)
    @eval Main struct PathContinuationDummyModel <: Models.AbstractQCDModel end
end

function _toy_point(param_value::Real, pressure::Real;
    residual_norm::Real=1e-8,
    converged::Bool=true,
    step_index::Int=1,
)
    return Models.BranchPoint(
        Float64(param_value),
        Float64[Float64(param_value), Float64(pressure)],
        Float64[Float64(param_value)],
        Float64[0.0, 0.0, 0.0],
        Float64(pressure),
        Float64(param_value),
        Float64(residual_norm),
        converged,
        :toy,
        step_index,
    )
end

function _toy_branch(branch_id::Symbol, pressures::AbstractVector{<:Real};
    residual_norm::Real=1e-8,
    converged::Bool=true,
)
    points = Models.BranchPoint[
        _toy_point(idx, p; residual_norm=residual_norm, converged=converged, step_index=idx)
        for (idx, p) in enumerate(pressures)
    ]
    return Models.ContinuationBranch(branch_id, Symbol(branch_id, "_anchor"), :toy, points, :complete, NamedTuple())
end

@testset "PathContinuation" begin
    @test isdefined(Models, :FixedAsymmetricRhoPath)
    @test isdefined(Models, :SeedContinuation)
    @test isdefined(Models, :PALCContinuation)
    @test isdefined(Models, :GroundStateBranchPolicy)
    @test isdefined(Models, :solve_path)
    @test isdefined(Models, :apply_branch_policy)

    @testset "FixedAsymmetricRhoPath validates parameters" begin
        path = Models.FixedAsymmetricRhoPath(1.0, [0.1, 0.2]; ud_ratio_target=0.876, s_target=0.0, xi=0.0, p_num=8, t_num=4)
        @test path.T_fm == 1.0
        @test path.rho_values == [0.1, 0.2]
        @test path.ud_ratio_target == 0.876
        @test path.p_num == 8
        @test path.t_num == 4

        @test_throws ArgumentError Models.FixedAsymmetricRhoPath(0.0, [0.1])
        @test_throws ArgumentError Models.FixedAsymmetricRhoPath(1.0, Float64[])
        @test_throws ArgumentError Models.FixedAsymmetricRhoPath(1.0, [-0.1])
        @test_throws ArgumentError Models.FixedAsymmetricRhoPath(1.0, [0.1]; ud_ratio_target=0.0)
        @test_throws ArgumentError Models.FixedAsymmetricRhoPath(1.0, [0.1]; p_num=0)
        @test_throws ArgumentError Models.FixedAsymmetricRhoPath(1.0, [0.1]; t_num=0)
    end

    @testset "GroundStateBranchPolicy selects pressure max under residual constraints" begin
        low = _toy_branch(:low, [8.0, 7.0])
        high = _toy_branch(:high, [10.0, 11.0])
        invalid = _toy_branch(:invalid_high, [12.0, 13.0]; residual_norm=1e-3)

        result = Models.apply_branch_policy(
            Models.ContinuationBranch[low, high, invalid],
            Models.GroundStateBranchPolicy(1e-3, 1e-6),
        )

        @test length(result.branches) == 3
        @test length(result.selections) == 2
        @test all(selection -> selection.selected_branch_id === :high, result.selections)
        @test all(selection -> selection.selection_reason === :pressure_max_under_constraints, result.selections)
        @test all(selection -> selection.candidate_branch_count == 3, result.selections)
        @test result.selections[1].runner_up_branch_id === :low
        @test result.selections[1].pressure_gap == 2.0
    end

    @testset "ContinuationBranchPolicy preserves requested branch" begin
        low = _toy_branch(:low, [8.0, 7.0])
        high = _toy_branch(:high, [10.0, 11.0])

        result = Models.apply_branch_policy(
            Models.ContinuationBranch[low, high],
            Models.ContinuationBranchPolicy(:low),
        )

        @test length(result.selections) == 2
        @test all(selection -> selection.selected_branch_id === :low, result.selections)
        @test all(selection -> selection.selection_reason === :continuation_branch, result.selections)
        @test result.selections[1].selected_pressure == 8.0
        @test result.selections[1].runner_up_branch_id === :high
        @test result.selections[1].pressure_gap == -2.0
        @test_throws ArgumentError Models.apply_branch_policy(
            Models.ContinuationBranch[low, high],
            Models.ContinuationBranchPolicy(:missing),
        )
    end

    @testset "AllBranchesPolicy returns all branches without single selection" begin
        low = _toy_branch(:low, [8.0])
        high = _toy_branch(:high, [10.0])

        result = Models.apply_branch_policy(
            Models.ContinuationBranch[low, high],
            Models.AllBranchesPolicy(),
        )

        @test [branch.branch_id for branch in result.branches] == [:low, :high]
        @test isempty(result.selections)
    end

    @testset "PALCContinuation is config-only in root environment" begin
        model = Main.PathContinuationDummyModel()
        path = Models.FixedAsymmetricRhoPath(1.0, [0.1])
        @test_throws Models.UnsupportedCapabilityError Models.solve_path(
            model,
            path;
            continuation_strategy=Models.PALCContinuation(),
        )
    end

    @testset "PathSolveResult converts to stable named tuples" begin
        branch = _toy_branch(:high, [10.0])
        policy_result = Models.apply_branch_policy(
            Models.ContinuationBranch[branch],
            Models.GroundStateBranchPolicy(),
        )
        diagnostics = Models.PathDiagnostics(:seed_continuation, :composite_anchor, :ground_state, 0.0, 1, 1, 0, 0, String[])
        path = Models.FixedAsymmetricRhoPath(1.0, [0.1])
        result = Models.PathSolveResult(path, policy_result.branches, policy_result.selections, Models.AnchorRoot[], diagnostics)

        nt = Models.to_namedtuple(result)
        @test haskey(nt, :branches)
        @test haskey(nt, :selections)
        @test haskey(nt, :anchors)
        @test haskey(nt, :diagnostics)
        @test nt.branches[1].branch_id === :high
        @test nt.selections[1].selected_branch_id === :high
        @test nt.diagnostics.continuation_backend === :seed_continuation
    end
end
