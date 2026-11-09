using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "FixedAsymmetricRho single-pass candidate search" begin
    M = Models
    mode = M.FixedAsymmetricRho(0.375, 0.876, 0.0)
    T = 120.0 / Main.Constants_PNJL.ħc_MeV_fm
    seeds = M.get_all_seeds(M.MultiSeed(), [T], mode)
    mu0 = [0.2, 0.2, 0.2]
    plan = [(; seed) for seed in vcat(seeds, seeds)]
    flattened = M._fixedasymrho_single_pass_plan(mode, T, plan, mu0)
    expected = Set(Tuple(s) for a in plan for s in M._fixedasymrho_inner_seeds(mode, T, a.seed, mu0))
    @test Set(Tuple(a.seed) for a in flattened) == expected
    @test length(flattened) == length(expected)
    @test M._fixedasymrho_single_pass_plan(mode, T, reverse(plan), mu0) == flattened
    @test all(a.method === :trust_region for a in flattened)
    for a in flattened
        @test M._fixedasymrho_inner_seeds(mode, T, a.seed, mu0; single_seed=true) == [a.seed]
    end

    model = M.create_model(:PNJL)
    @test_throws ArgumentError M._fixedasymrho_problem_spec_forward_solve(model, mode, T;
        seed_guess=seeds[1], fixedasymrho_seed_policy=:unknown)
    @test_throws ArgumentError M._fixedasymrho_problem_spec_forward_solve(model, mode, T;
        seed_guess=seeds[1], fixedasymrho_seed_policy=:single_pass, evaluate_all_attempts=false)

    # Exercise the real orchestrator with a cheap multi-root solver stand-in.
    # A smaller residual on the low-pressure root must not hide the other root.
    observed = NamedTuple[]
    calls = Vector{Vector{Float64}}()
    selector = candidates -> begin
        append!(observed, candidates)
        M.select_pressure_max_candidate(candidates)
    end
    kwargs = Dict{Symbol,Any}(:seed_guess => seeds[1], :seed_candidates => seeds,
        :mu0 => mu0, :rho0 => 0.16, :residual_norm_max => 1e-6,
        :evaluate_all_attempts => true, :fixedasymrho_seed_policy => :single_pass,
        :selector => selector)
    fake_solve = local_kwargs -> begin
        @test local_kwargs[:single_seed]
        seed = local_kwargs[:seed_guess]
        push!(calls, copy(seed))
        high = seed[1] > -1.0
        M._build_mode_result_from_outer_state(
            [-1.0, -1.0, -2.0, 0.2, 0.2], mu0,
            high ? 22.0 : 18.0, 0.375, 0.2, 0.4, [1.5, 1.5, 2.5],
            high ? 1e-8 : 1e-12; converged=true, iterations=2)
    end
    result = M._governed_nonrho_problem_spec_forward_solve(model, mode, T, kwargs,
        "FixedAsymmetricRho", :asym_attempt_origin, fake_solve)
    @test result.converged
    @test result.pressure == 22.0
    @test any(c.pressure == 18.0 for c in observed)
    @test any(c.pressure == 22.0 for c in observed)
    @test length(calls) == length(Set(Tuple(s) for s in calls))
    @test result.candidate_count == length(calls)
end
