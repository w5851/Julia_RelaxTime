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
    @test isdefined(Models, :SeedPoolAnchor)
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

    @testset "GroundStateBranchPolicy marks pressure-degenerate selections" begin
        high = _toy_branch(:high, [10.0])
        near = _toy_branch(:near, [9.9995])

        result = Models.apply_branch_policy(
            Models.ContinuationBranch[high, near],
            Models.GroundStateBranchPolicy(1e-3, 1e-6),
        )

        @test length(result.selections) == 1
        @test result.selections[1].selected_branch_id === :high
        @test result.selections[1].runner_up_branch_id === :near
        @test result.selections[1].selection_reason === :pressure_degenerate_under_constraints
        @test result.selections[1].pressure_gap ≈ 5e-4 atol=1e-12
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

    @testset "PALC rejects models without supported equilibrium residuals" begin
        model = Main.PathContinuationDummyModel()
        path = Models.FixedAsymmetricRhoPath(1.0, [0.1])
        @test_throws Models.UnsupportedCapabilityError Models.solve_path(
            model,
            path;
            continuation_strategy=Models.PALCContinuation(),
        )
    end

    @testset "New path units, input validation and scan dispatch" begin
        mu = Models.FixedMuPath(1.0, [1.8, 1.5, 1.5]; p_num=8, t_num=4)
        rho = Models.FixedRhoPath(1.0, [0.0, 0.5]; p_num=8, t_num=4)
        @test Models.to_namedtuple(mu).mu_values == [1.8, 1.5, 1.5]
        @test Models.to_namedtuple(rho).path_kind === :fixed_rho
        @test_throws ArgumentError Models.FixedMuPath(1.0, [NaN])
        @test_throws ArgumentError Models.FixedRhoPath(1.0, [-0.1])
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state, nothing, 21) === :palc
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state,
            nothing, 23; constraint_mode=:fixed_asymmetric_rho) === :palc
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state,
            nothing, 31; constraint_mode=:fixed_mu) === :palc
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state,
            nothing, 30; constraint_mode=:fixed_mu) === :legacy
        @test Models.resolve_scan_continuation(:legacy, :PNJL, :tensor_gauss, :ground_state,
            nothing, 31; constraint_mode=:fixed_mu) === :legacy
        @test Models.resolve_scan_continuation(:auto, :RPNJL, :tensor_gauss, :ground_state,
            nothing, 31; constraint_mode=:fixed_mu) === :legacy
        @test Models.resolve_scan_continuation(:auto, :PNJL, :rs_reduced_adaptive, :ground_state,
            nothing, 31; constraint_mode=:fixed_mu) === :legacy
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state,
            nothing, 31; constraint_mode=:fixed_mu, controls_supported=false) === :legacy
        @test Models.resolve_scan_continuation(:auto, :RPNJL, :tensor_gauss, :ground_state, nothing, 23) === :legacy
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :constrained_manifold, nothing, 23) === :legacy
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state, identity, 23) === :legacy
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state, nothing, 20) === :legacy
        @test Models.resolve_scan_continuation(:palc, :PNJL, :tensor_gauss, :ground_state, nothing, 21) === :palc
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state, nothing, 2) === :legacy
        @test Models.resolve_scan_continuation(:auto, :NJL, :tensor_gauss, :ground_state, nothing, 21) === :legacy
        @test Models.resolve_scan_continuation(:legacy, :PNJL, :tensor_gauss, :ground_state, nothing, 21) === :legacy
        @test_throws ArgumentError Models.resolve_scan_continuation(:palc, :PNJL, :rs_reduced_adaptive, :ground_state, nothing, 21)
        @test_throws ArgumentError Models.resolve_scan_continuation(:unknown, :PNJL, :tensor_gauss, :ground_state, nothing, 21)
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state,
            nothing, 31; controls_supported=false) === :legacy
        @test_throws ArgumentError Models.resolve_scan_continuation(:palc, :PNJL, :tensor_gauss,
            :ground_state, nothing, 31; controls_supported=false)
        @test Models._scan_path_controls_supported((; iterations=300, residual_norm_max=1e-7))
        @test !Models._scan_path_controls_supported((; solver=:custom))
        @test !Models._scan_path_controls_supported((; evaluate_all_attempts=false))
    end

    @testset "Path controls reject unsupported options before solving" begin
        model = Models.create_model(:PNJL)
        path = Models.FixedMuPath(1.0, [0.0, 0.1])
        @test_throws ArgumentError Models._path_options(model, path; invalid_option=1)
        @test_throws ArgumentError Models._path_options(model, path; thermo_quadrature_rtol=-1.0)
        @test Models._path_options(model, path; iterations=30).iteration_limit == 30
        @test all(seed -> length(seed) == 5, Models._path_catalog(path, 0.1))
        @test all(seed -> length(seed) == 8, Models._path_catalog(Models.FixedRhoPath(1.0, [0.5]), 0.5))
        @test all(seed -> length(seed) == 8, Models._path_catalog(Models.FixedAsymmetricRhoPath(1.0, [0.5]), 0.5))
        @test_throws ArgumentError Models._path_options(model, path; evaluate_all_attempts=false)
    end

    @testset "Temperature paths retain constraint targets and input order" begin
        modes = (Models.FixedEntropy(0.5), Models.FixedSigma(12.0),
            Models.FixedMuBConservedCharges(4.5))
        for mode in modes
            path = Models.ConstraintTPath(mode, [0.7, 0.5, 0.6]; p_num=24, t_num=8)
            @test Models._path_values(path) == [0.7, 0.5, 0.6]
            @test Models._path_temperature(path, 0.6) == 0.6
            @test Models._path_mode(path, 0.6) === mode
            @test Models.to_namedtuple(path).T_values == [0.7, 0.5, 0.6]
            @test all(seed -> length(seed) == 8, Models._path_catalog(path, 0.6))
            @test_throws ArgumentError Models.ConstraintTPath(mode, [0.0])
            @test_throws ArgumentError Models.ConstraintTPath(mode, [NaN])
            @test_throws ArgumentError Models.ConstraintTPath(mode, Float64[])
            @test_throws ArgumentError Models.ConstraintTPath(mode, [0.5]; xi=-1.0)
        end
        @test_throws ArgumentError Models.SingleSeedAnchor(NaN, ones(8))
        @test_throws ArgumentError Models.SingleSeedAnchor(0.5, [Inf])
        model = Models.create_model(:PNJL)
        path = Models.ConstraintTPath(Models.FixedSigma(12.0), [0.5, 0.6])
        @test_throws ArgumentError Models.solve_path(model, path; continuation=:palc,
            anchor_strategy=Models.SingleSeedAnchor(0.5, ones(5)))
        @test_throws ArgumentError Models.solve_path(model, path; continuation=:palc,
            anchor_strategy=Models.SingleSeedAnchor(0.7, ones(8)))
        @test_throws ArgumentError Models.solve_path(model, path; continuation=:palc,
            anchor_strategy=Models.SingleSeedAnchor(0.5, ones(8)), anchor_values=[0.5])
    end

    @testset "FixedMu scan preserves its historical density unit" begin
        value = Models.SolverResult(Models.FixedMu(), true, zeros(5), zeros(5), zeros(3),
            -1.0, 1.0, 2.0, 3.0, 4.0, ones(3), 0, 1e-9, 0.0)
        adapted = Models.TmuScan._fixedmu_scan_units(value)
        @test adapted.rho_norm == 2.0 * Main.Constants_PNJL.ρ0_inv_fm3
        @test adapted.pressure == value.pressure
        @test adapted.solution == value.solution
        @test Models.TmuScan._fixedmu_scan_units(nothing) === nothing
    end

    @testset "Fold crossings remain separate target seeds" begin
        curve = [(param=0.0, solution=[-2.0]), (param=2.0, solution=[0.0]),
            (param=0.0, solution=[2.0])]
        seeds = Models._palc_target_seeds(curve, 1.0)
        @test length(seeds) == 2
        @test seeds == [[-1.0], [1.0]]
        @test isempty(Models._palc_target_seeds(curve, 3.0))
        @test length(Models._palc_target_seeds(curve, 3.0; nearest=true)) == 1
        first_point = _toy_point(1.0, 2.0)
        second_point = _toy_point(1.0, 1.5)
        folded = Models.ContinuationBranch(:fold, :anchor, :bifurcationkit,
            [first_point, second_point], :complete, (;))
        selection = only(Models.apply_branch_policy([folded], Models.GroundStateBranchPolicy()).selections)
        @test selection.candidate_branch_count == 2
        @test selection.pressure_gap == 0.5
    end

    @testset "Pseudo-arclength crosses a fold on an analytic equation" begin
        residual = (x, value) -> [x[1]^2 - value]
        work = Models.PALCWorkTelemetry()
        strategy = Models.PALCContinuation(ds=0.05, dsmax=0.1, max_steps=100, newton_tol=1e-10)
        curve = Models._palc_curve(residual, (solution=[-1.0],), 1.0, (0.0, 1.1), strategy, work)
        @test any(point -> point.solution[1] > 0.5, curve)
        @test maximum(abs(point.solution[1]^2 - point.param) for point in curve) < 1e-7
        @test length(Models._palc_target_seeds(curve, 0.5)) >= 2
        @test work.continuation_residual_calls > 0
        @test work.continuation_jacobians > 0
    end

    @testset "SeedContinuation diagnostics do not claim unused anchor discovery" begin
        model = Main.PathContinuationDummyModel()
        path = Models.FixedAsymmetricRhoPath(1.0, [0.1])
        @test_throws ArgumentError Models.solve_path(
            model,
            path;
            anchor_strategy=Models.MultiSeedAnchor(),
        )
    end

    @testset "Path branch jump metrics count state discontinuities" begin
        points = Models.BranchPoint[
            Models.BranchPoint(0.1, [0.0, 0.0], [0.0], [0.0, 0.0, 0.0], 1.0, 0.1, 1e-8, true, :toy, 1),
            Models.BranchPoint(0.2, [0.1, 0.1], [0.1], [0.0, 0.0, 0.0], 1.1, 0.2, 1e-8, true, :toy, 2),
            Models.BranchPoint(0.3, [2.0, 2.0], [2.0], [0.0, 0.0, 0.0], 1.2, 0.3, 1e-8, true, :toy, 3),
        ]
        branch = Models.ContinuationBranch(:toy, :toy_anchor, :toy, points, :complete, NamedTuple())

        jumps = Models._path_branch_jump_metrics(Models.ContinuationBranch[branch]; tol=0.5)

        @test jumps.count == 1
        @test jumps.max_jump ≈ sqrt(1.9^2 + 1.9^2)
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
        @test nt.path.path_kind === :fixed_asymmetric_rho
        @test nt.path.rho_values == [0.1]
        @test nt.branches[1].branch_id === :high
        @test nt.selections[1].selected_branch_id === :high
        @test nt.diagnostics.continuation_backend === :seed_continuation
    end
end
