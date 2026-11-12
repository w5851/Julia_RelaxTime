using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

if !isdefined(Main, :PathContinuationDummyModel)
    @eval Main struct PathContinuationDummyModel <: Models.AbstractQCDModel end
end

if !isdefined(Main, :ConstraintPathRecordingModel)
    @eval Main struct ConstraintPathRecordingModel <: Models.AbstractQCDModel
        calls::Vector{NamedTuple}
    end
end

if !isdefined(Main, :EndpointPathTestModel)
    @eval Main struct EndpointPathTestModel <: Models.AbstractQCDModel
        reject_attempts::Int
        calls::Vector{Symbol}
    end
end

if !isdefined(Main, :CorrectorPathTestModel)
    @eval Main struct CorrectorPathTestModel <: Models.AbstractQCDModel end
end

function Models._path_candidate(::Main.CorrectorPathTestModel, path, value, solution, residual, iterations, options)
    residual_norm = maximum(abs, residual(solution,value))
    return (; solution, iterations, residual_norm, converged=residual_norm < options.residual_tol)
end

Models._path_options(::Main.EndpointPathTestModel,path;kwargs...) = nothing
Models._path_residual(::Main.EndpointPathTestModel,path,options) = (x,value) -> x .- value
function Models._path_root(model::Main.EndpointPathTestModel,path,value,seed,residual,options,work;
    source::Symbol=:correction,work_telemetry=nothing,workspace=nothing)
    push!(model.calls,source)
    work.root_solves += 1
    source === :anchor && (work.anchor_solves += 1)
    source === :correction && (work.correction_solves += 1)
    source === :fallback && (work.fallback_solves += 1)
    (length(model.calls) <= model.reject_attempts || length(model.calls) == -model.reject_attempts) && return nothing
    return (;solution=fill(value,length(seed)),x_state=fill(value,5),mu_vec=fill(value,3),
        pressure=value,rho_norm=value,residual_norm=0.0,converged=true)
end

function Models.solve(model::Main.ConstraintPathRecordingModel, mode::Models.FixedSigma,
    T_fm::Real; kwargs...)
    push!(model.calls, (;T_fm, kwargs...))
    solution = copy(kwargs[:seed_guess])
    solution[6:8] .+= 0.1
    return Models.SolverResult(mode, true, solution, solution[1:5], solution[6:8],
        -1.0, 1.0, 1.0, 2.0, 3.0, ones(3), 1, 1e-9, 0.0)
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
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state, nothing, 21) === :natural
        @test Models.resolve_scan_continuation(:auto, :PNJL, :tensor_gauss, :ground_state,
            nothing, 23; constraint_mode=:fixed_asymmetric_rho) === :natural
        @test Models.resolve_scan_continuation(:natural,:PNJL,:tensor_gauss,:ground_state,nothing,3) === :natural
        @test_throws ArgumentError Models.resolve_scan_continuation(:natural,:PNJL,:tensor_gauss,
            :ground_state,nothing,31;constraint_mode=:fixed_mu)
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
        @test_throws ArgumentError Models._path_options(model, path; corrector_method=:unknown)
    end

    @testset "Optimized defaults follow the qualified path and preserve explicit controls" begin
        model = Models.create_model(:PNJL)
        for path in (Models.FixedMuPath(1.0,[0.0,0.1]), Models.FixedRhoPath(1.0,[0.1,0.2]),
            Models.FixedAsymmetricRhoPath(1.0,[0.1,0.2]),
            Models.ConstraintTPath(Models.FixedSigma(12.0),[0.5,0.6]),
            Models.ConstraintTPath(Models.FixedMuBConservedCharges(4.5),[0.5,0.6]))
            options = Models._path_options(model,path)
            @test options.corrector_method === :newton
            strategy = Models._resolved_palc_strategy(model,path,Models.PALCContinuation())
            @test strategy.target_predictor === :cubic
            @test strategy.jacobian === (path isa Models.ConstraintTPath{Models.FixedSigma} ?
                :forwarddiff : :finite_difference)
            explicit = Models._resolved_palc_strategy(model,path,
                Models.PALCContinuation(jacobian=:finite_difference,target_predictor=:linear))
            @test explicit.jacobian === :finite_difference && explicit.target_predictor === :linear
            @test Models._path_options(model,path;corrector_method=:trust_region).corrector_method === :trust_region
        end
        entropy = Models.ConstraintTPath(Models.FixedEntropy(0.5),[0.5,0.6])
        rpnjl = Models.create_model(:RPNJL)
        for (m,p) in ((model,entropy),(rpnjl,Models.FixedMuPath(1.0,[0.0,0.1])))
            @test Models._path_options(m,p).corrector_method === :trust_region
            strategy = Models._resolved_palc_strategy(m,p,Models.PALCContinuation())
            @test strategy.jacobian === :finite_difference && strategy.target_predictor === :linear
        end
    end

    @testset "Endpoint starts preserve path constraints and explicit catalogs" begin
        model = Models.create_model(:PNJL)
        path = Models.FixedMuPath(120/197.327, [1.4,1.6]; p_num=24,t_num=8)
        @test Models._default_path_anchor(model,path) isa Models.EndpointAnchor
        @test Models._default_path_anchor(Models.create_model(:RPNJL),path) isa Models.SeedPoolAnchor
        @test Models._default_path_anchor(model,Models.FixedMuPath(1.0,[1.5])) isa Models.SeedPoolAnchor
        attempts = Models._endpoint_anchor_attempts(path)
        @test first(attempts).value == maximum(path.mu_values)
        @test all(a -> length(a.seed)==5, attempts)
        @test first(attempts).seed == Models.get_seed(Models.DefaultSeed(phase_hint=:quark),
            [path.T_fm,maximum(path.mu_values)],Models.FixedMu())
        thermal = Models.ConstraintTPath(Models.FixedSigma(12.0),[0.8,0.5,0.6])
        @test first(Models._endpoint_anchor_attempts(thermal)).value == 0.8
        @test all(a -> length(a.seed)==8, Models._endpoint_anchor_attempts(thermal))
        entropy = Models.ConstraintTPath(Models.FixedEntropy(0.5),[0.5,0.6])
        @test Models._default_path_anchor(model,entropy) isa Models.SeedPoolAnchor
        @test_throws ArgumentError Models.solve_path(model,path;
            continuation_strategy=Models.PALCContinuation(),anchor_strategy=Models.EndpointAnchor(),
            anchor_values=[1.5])
    end

    @testset "Endpoint recovery runs only after a failed start or incomplete path" begin
        path = Models.FixedMuPath(1.0,[0.1,0.2])
        model = Main.EndpointPathTestModel(1,Symbol[])
        work = Models.PALCWorkTelemetry()
        result = Models.solve_path(model,path;continuation_strategy=Models.PALCContinuation(),
            anchor_strategy=Models.EndpointAnchor(),path_work=work)
        @test result.diagnostics.anchor_strategy === :endpoint_anchor
        @test result.diagnostics.failure_count == 0
        @test length(result.anchors) == 1
        @test work.anchor_solves == 2
        @test work.fallback_solves == 0
        @test [s.param_value for s in result.selections] == [0.1,0.2]

        failed = Main.EndpointPathTestModel(typemax(Int),Symbol[])
        result = Models.solve_path(failed,path;continuation_strategy=Models.PALCContinuation(),
            anchor_strategy=Models.EndpointAnchor())
        @test result.diagnostics.anchor_strategy === :endpoint_catalog_recovery
        @test result.diagnostics.failure_count == 2
        @test isempty(result.anchors)
        @test :fallback in failed.calls

        short = Main.EndpointPathTestModel(0,Symbol[])
        result = Models.solve_path(short,Models.FixedMuPath(1.0,[0.1,1.0,2.0]);
            continuation_strategy=Models.PALCContinuation(max_steps=1),
            anchor_strategy=Models.EndpointAnchor())
        @test result.diagnostics.anchor_strategy === :endpoint_catalog_recovery
        @test result.diagnostics.failure_count == 0
        @test length(result.anchors) > 1
        @test all(s -> s.selected_branch_id !== nothing,result.selections)
    end

    @testset "Natural density continuation preserves target order and full state dimension" begin
        model = Main.EndpointPathTestModel(0,Symbol[])
        path = Models.FixedRhoPath(1.0,[0.1,0.3,0.2,0.1])
        result = Models.solve_path(model,path;continuation_strategy=Models.NaturalContinuation())
        @test result.diagnostics.continuation_backend === :natural_parameter
        @test result.diagnostics.failure_count == 0
        @test length(result.anchors) == 1
        @test length(model.calls) == 3
        @test all(p -> length(p.solution)==8,only(result.branches).points)
        @test [p.param_value for p in only(result.branches).points] == [0.3,0.2,0.1]
        @test Models._path_values(result.path) == [0.1,0.3,0.2,0.1]
        @test_throws ArgumentError Models.solve_path(model,Models.FixedMuPath(1.0,[0.1,0.2]);
            continuation_strategy=Models.NaturalContinuation())
        @test_throws ArgumentError Models.solve_path(model,path;
            continuation_strategy=Models.NaturalContinuation(),anchor_strategy=Models.SeedPoolAnchor())
        interrupted = Main.EndpointPathTestModel(-2,Symbol[])
        recovered = Models.solve_path(interrupted,path;continuation_strategy=Models.NaturalContinuation())
        @test recovered.diagnostics.continuation_backend === :natural_palc_recovery
        @test recovered.diagnostics.failure_count == 0
        @test all(s -> s.selected_branch_id !== nothing,recovered.selections)
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
            model = Models.create_model(:PNJL)
            @test Models._constraint_path_route(:palc, model, path; iterations=300,
                path_work=Models.PALCWorkTelemetry()) === :palc
            @test Models._constraint_path_route(:auto, model, path) === :legacy
            long_path = Models.ConstraintTPath(mode, collect(range(0.5, 0.7; length=21)))
            expected = mode isa Models.FixedEntropy ? :legacy : :palc
            @test Models._constraint_path_route(:auto, model, long_path) === expected
            @test Models._constraint_path_route(:auto, model, long_path; solver=:custom) === :legacy
            @test Models._constraint_path_route(:auto, Models.create_model(:RPNJL), long_path) === :legacy
            @test_throws ArgumentError Models.ConstraintTPath(mode, [0.0])
            @test_throws ArgumentError Models.ConstraintTPath(mode, [NaN])
            @test_throws ArgumentError Models.ConstraintTPath(mode, Float64[])
            @test_throws ArgumentError Models.ConstraintTPath(mode, [0.5]; xi=-1.0)
        end
        @test_throws ArgumentError Models.SingleSeedAnchor(NaN, ones(8))
        @test_throws ArgumentError Models.ConstraintTPath(Models.FixedEntropy(-0.5), [0.5])
        @test_throws ArgumentError Models.ConstraintTPath(Models.FixedSigma(0.0), [0.5])
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

    @testset "Legacy temperature path forwards custom controls and the previous root" begin
        model = Main.ConstraintPathRecordingModel(NamedTuple[])
        path = Models.ConstraintTPath(Models.FixedSigma(12.0), [0.5, 0.6])
        seed = [-1.0,-1.0,-2.0,0.2,0.2,1.0,1.0,1.0]
        result = Models.solve_path(model, path; continuation=:legacy,
            seed_guess=seed, evaluate_all_attempts=true, iterations=21, mu0=0.4)
        @test model.calls[1].seed_guess == seed
        @test model.calls[2].seed_guess == result.branches[1].points[1].solution
        @test all(call -> call.evaluate_all_attempts && call.iterations == 21 && call.mu0 == 0.4, model.calls)
        @test [call.T_fm for call in model.calls] == [0.5,0.6]
        @test result.diagnostics.failure_count == 0
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

    @testset "Cubic target prediction stays on one monotone segment" begin
        points = [(param=p, solution=[p^3]) for p in 0.0:0.25:1.0]
        cubic = only(Models._palc_target_seeds(points, 0.4; predictor=:cubic))
        linear = only(Models._palc_target_seeds(points, 0.4; predictor=:linear))
        @test cubic[1] ≈ 0.4^3 atol=1e-14
        @test abs(cubic[1]-0.4^3) < abs(linear[1]-0.4^3)
        folded = [(param=0.0, solution=[-2.0]), (param=1.0, solution=[-1.0]),
            (param=2.0, solution=[0.0]), (param=1.0, solution=[1.0]),
            (param=0.0, solution=[2.0])]
        @test Models._palc_target_seeds(folded, 1.5; predictor=:cubic) == [[-0.5], [0.5]]
        @test_throws ArgumentError Models.PALCContinuation(target_predictor=:unknown)
        @test_throws ArgumentError Models.PALCContinuation(jacobian=:unknown)
    end

    @testset "PALC Jacobians reuse storage and update the parameter" begin
        residual = (x,p) -> [x[1]^2+p*x[2], sin(x[2])+p*x[1]]
        x = [0.4, -0.2]
        for backend in (:finite_difference, :forwarddiff)
            work = Models.PALCWorkTelemetry()
            jac = Models._palc_jacobian(residual, x, 0.3, work, backend)
            storage = jac(x, (;value=0.3))
            @test storage ≈ [0.8 0.3; 0.3 cos(-0.2)] atol=1e-9
            @test jac(x, (;value=0.7)) === storage
            @test storage ≈ [0.8 0.7; 0.7 cos(-0.2)] atol=1e-9
            @test work.continuation_jacobians == 2
        end
    end

    @testset "Root workspace invalidates a changed parameter at the same seed" begin
        residual = (x,p) -> [x[1]^2-p, x[2]-p]
        work = Models.PALCWorkTelemetry()
        workspace = Models.PathRootWorkspace(residual, 1.0, [1.0,1.0], work)
        first_root = Models._workspace_root!(workspace, 1.0, [1.0,1.0], 50)
        second_root = Models._workspace_root!(workspace, 4.0, [1.0,1.0], 50)
        @test first_root.zero ≈ [1.0,1.0] atol=1e-9
        @test second_root.zero ≈ [2.0,4.0] atol=1e-9
        third_root = Models._workspace_root!(workspace, 1.0, [1.0,1.0], 50)
        @test third_root.f_calls == first_root.f_calls
        @test third_root.g_calls == first_root.g_calls
        @test work.workspace_builds == 1
        @test first_root.zero == [1.0,1.0]
        local_root = Models._newton_corrector!(workspace, 4.0, [1.8,3.5], 50)
        @test local_root.accepted
        @test local_root.zero ≈ [2.0,4.0] atol=1e-9
        @test local_root.g_calls == local_root.iterations
        already_solved = Models._newton_corrector!(workspace, 4.0, [2.0,4.0], 50)
        @test already_solved.accepted
        @test already_solved.g_calls == 0
        singular = Models._newton_corrector!(workspace, 4.0, [0.0,4.0], 50)
        @test !singular.accepted
    end

    @testset "Warm endpoint is a guess with ordinary recovery" begin
        path = Models.FixedRhoPath(0.7, [0.1,0.2,0.3])
        model = Main.EndpointPathTestModel(1, Symbol[])
        result = Models.solve_path(model, path; continuation_strategy=Models.NaturalContinuation(),
            anchor_strategy=Models.EndpointAnchor(fill(0.25, 8)))
        @test result.diagnostics.failure_count == 0
        @test result.diagnostics.anchor_solve_count == 2
        @test Models._reusable_endpoint_anchor(result, Models.FixedRhoPath(0.8,[0.1,0.3])).seed == fill(0.3,8)
        @test Models._reusable_endpoint_anchor(result, Models.FixedRhoPath(0.8,[0.1,0.3];xi=0.4)).seed === nothing
        @test Models._reusable_endpoint_anchor(result, Models.FixedRhoPath(0.8,[0.1,0.4])).seed === nothing
        @test Models._reusable_endpoint_anchor(result, Models.FixedMuPath(0.8,[0.1,0.3])).seed === nothing
        @test_throws ArgumentError Models._endpoint_anchor_attempts(path, Models.EndpointAnchor(ones(5)))
    end

    @testset "Local corrector failure restores trust-region from the original seed" begin
        model = Main.CorrectorPathTestModel()
        path = Models.FixedMuPath(1.0,[1.0])
        residual = (x,p) -> [x[1]^3-p]
        work = Models.PALCWorkTelemetry()
        options = (;iteration_limit=300,corrector_method=:newton,residual_tol=1e-6)
        seed = [0.01]
        workspace = Models.PathRootWorkspace(residual,1.0,seed,work)
        candidate = Models._path_root(model,path,1.0,seed,residual,options,work;workspace)
        @test candidate !== nothing && candidate.converged
        @test candidate.solution ≈ [1.0] atol=1e-8
        @test work.corrector_fallbacks == 1
        @test work.newton_corrections == 1
        @test seed == [0.01]
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

    @testset "Connected pitchfork branches require branch switching" begin
        residual = (x, value) -> [x[1] * (x[1]^2 - value)]
        work = Models.PALCWorkTelemetry()
        curve = Models._palc_curve(residual, (solution=[0.0],), -0.5, (-0.6, 1.0),
            Models.PALCContinuation(ds=0.03, dsmax=0.1, max_steps=100), work)
        @test maximum(point.param for point in curve) > 0.5
        @test all(point -> abs(point.solution[1]) < 1e-10, curve)
        @test residual([1.0], 1.0) == [0.0]
        @test residual([-1.0], 1.0) == [0.0]
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
