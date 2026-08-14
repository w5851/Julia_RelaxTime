using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "GenericRootEngine" begin
    @test isdefined(Models, :RootProblemSpec)
    @test isdefined(Models, :RootPolicy)
    @test isdefined(Models, :ContinuationState)
    @test isdefined(Models, :solve_root_with_policy)
    @test isdefined(Models, :solve_root_continuation)

    @testset "solve_root_with_policy solves scalar root" begin
        spec = Models.RootProblemSpec(
            (F, x, _) -> (F[1] = x[1]^2 - 2.0),
            nothing,
            1,
            :single,
        )
        result = Models.solve_root_with_policy(spec, [1.0])
        @test result.converged
        @test result.quality_tag in (:good, :fallback)
        @test isapprox(result.x[1], sqrt(2.0); atol=1e-8)
        @test length(result.diagnostics.attempts) >= 1
    end

    @testset "solve_root_continuation tracks across points" begin
        spec_factory = p -> Models.RootProblemSpec(
            (F, x, ctx) -> (F[1] = x[1] - ctx.target),
            (target=Float64(p),),
            1,
            :single,
        )
        tracker = Models.ContinuationState()
        results = Models.solve_root_continuation([1.0, 2.0, 3.0], spec_factory; tracker=tracker, x0=[0.0])

        @test length(results) == 3
        @test isapprox(results[end].x[1], 3.0; atol=1e-10)
        @test haskey(tracker.seed_by_branch, :default)
        @test isapprox(tracker.seed_by_branch[:default][1], 3.0; atol=1e-10)
    end

    @testset "domain_quality can veto quality gate" begin
        spec = Models.RootProblemSpec(
            (F, x, _) -> (F[1] = x[1] + 1.0),
            nothing,
            1,
            :single,
        )
        dq = (x, _) -> (ok=x[1] > 0.0, score=abs(x[1]), reason=:positive_required)
        result = Models.solve_root_with_policy(spec, [0.0]; domain_quality=dq)
        @test result.converged
        @test result.quality_tag in (:degraded, :bad)
    end

    @testset "callback solver path supports fallback" begin
        solve_once = function (method::Symbol, seed::Vector{Float64})
            if method === :newton
                return (mass=seed[1], gamma=seed[2], converged=false, residual_norm=1.0)
            end
            return (mass=2.0, gamma=0.0, converged=true, residual_norm=1e-12)
        end

        result = Models.solve_root_with_policy(solve_once, [1.0, 0.0])
        @test result.converged
        @test result.quality_tag in (:good, :fallback)
        @test isapprox(result.x[1], 2.0; atol=1e-12)
    end

    @testset "callback solver does not repeat identical fallback method" begin
        calls = Ref(0)
        solve_once = function (method::Symbol, seed::Vector{Float64})
            calls[] += 1
            return (mass=seed[1], gamma=seed[2], converged=false, residual_norm=1.0)
        end

        result = Models.solve_root_with_policy(
            solve_once,
            [1.0, 0.0];
            policy=Models.RootPolicy(
                primary_method=:trust_region,
                fallback_method=:trust_region,
                use_fallback=true,
            ),
        )
        @test calls[] == 1
        @test length(result.diagnostics.attempts) == 1
        @test !result.converged
    end

    @testset "callback solver prefers lower score over lower residual" begin
        solve_once = function (method::Symbol, seed::Vector{Float64})
            if method === :newton
                return (mass=10.0, gamma=0.0, converged=true, residual_norm=1e-12, score=10.0)
            end
            return (mass=11.0, gamma=0.0, converged=true, residual_norm=1e-8, score=0.1)
        end

        result = Models.solve_root_with_policy(solve_once, [1.0, 0.0])
        @test result.quality_tag == :fallback
        @test isapprox(result.x[1], 11.0; atol=1e-12)
        selected = result.diagnostics.attempts[result.diagnostics.selected_attempt]
        @test selected.score ≈ 0.1
    end
end
