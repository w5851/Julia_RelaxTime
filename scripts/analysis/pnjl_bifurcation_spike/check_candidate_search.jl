# Bounded remote checks for the two candidate-visibility counterexamples.
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "run_matched_pilot.jl"))

const P = MatchedFixedAsymRhoPilot

function check_candidate_search(output)
    original = P.JSON.parsefile(joinpath(Main.REPO_ROOT, "data", "outputs", "results",
        "analysis", "palc_matched_20261108", "pilot.json"))
    report = Any[]
    model = Models.create_model(:PNJL)
    @testset "single-pass numerical candidate visibility" begin
        for (case, target) in (("t120_rho035", 0.375), ("t130_rho080", 0.75))
            cfg = P.config(case, output; profile="local")
            T = cfg.T_MeV / P.U.HBARC_MEV_FM
            mode = P.mode_at(cfg, target)
            old = original["cases"][case]["methods"]["palc"]["repetitions"][end]["result"]
            old_roots = [r for r in old["candidates"] if r["accepted"] && abs(r["rho"]-target) <= 1e-12]
            low = old_roots[argmin(r["pressure"] for r in old_roots)]
            high = old_roots[argmax(r["pressure"] for r in old_roots)]
            seeds = vcat([Float64.(low["solution"]), Float64.(high["solution"])],
                Models.get_all_seeds(Models.MultiSeed(), [T], mode))
            measurements = Any[]
            selected_solutions = Vector{Float64}[]
            for permutation in (seeds, reverse(seeds), circshift(seeds, 2))
                candidates = NamedTuple[]
                selector = roots -> begin
                    append!(candidates, roots)
                    Models.select_pressure_max_candidate(roots)
                end
                telemetry = Models.SolverWorkTelemetry()
                selected = Models.solve_constraint(model, mode, T;
                    seed_guess=seeds[1], seed_candidates=permutation, mu0=[0.2, 0.2, 0.2],
                    selector, fixedasymrho_seed_policy=:single_pass, evaluate_all_attempts=true,
                    xi=0.0, p_num=8, t_num=4, rho0=Models.ρ0,
                    residual_norm_max=P.RESIDUAL_TOL, iterations=300, work_telemetry=telemetry)
                accepted = filter(c -> c.converged && c.hard_constraint_ok, candidates)
                @test selected.converged
                @test selected.residual_norm <= P.RESIDUAL_TOL
                @test any(norm(c.solution-Float64.(low["solution"])) <= P.ROOT_DISTANCE for c in accepted)
                @test any(norm(c.solution-Float64.(high["solution"])) <= P.ROOT_DISTANCE for c in accepted)
                @test selected.pressure >= high["pressure"]-1e-6
                @test telemetry.fixedasymrho_inner_solves == selected.candidate_count
                @test telemetry.fixedasymrho_inner_solves == telemetry.governed_attempts
                points = NamedTuple[]
                work = P.PilotWork()
                for root in accepted
                    P.add_distinct!(points, P.checked_point(model, cfg, target, root.solution, true, work))
                end
                push!(selected_solutions, copy(selected.solution))
                push!(measurements, (; selected_pressure=selected.pressure,
                    selected_solution=selected.solution, roots=points,
                    work=Models.solver_work_snapshot(telemetry)))
            end
            @test all(norm(x-first(selected_solutions)) <= P.ROOT_DISTANCE for x in selected_solutions)

            # Default dispatch must preserve the previous nested implementation.
            common = (; seeds, mu0=[0.2,0.2,0.2], xi=0.0, p_num=8, t_num=4,
                iterations=300, residual_norm_max=P.RESIDUAL_TOL, evaluate_all_attempts=true)
            default = Models.solve_multi(model, mode, T; common...)
            nested = Models.solve_multi(model, mode, T; common..., fixedasymrho_seed_policy=:nested)
            @test default.converged == nested.converged
            @test norm(default.solution-nested.solution) <= P.ROOT_DISTANCE
            @test abs(default.pressure-nested.pressure) <= 1e-6
            push!(report, (; case, target, permutations=measurements))
        end
    end
    mkpath(dirname(output))
    P.U._write_json(output, (; schema="fixedasymrho_candidate_checks_v1",
        source_sha=readchomp(`git -C $(Main.REPO_ROOT) rev-parse HEAD`), cases=report))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || throw(ArgumentError("usage: check_candidate_search.jl OUTPUT_JSON"))
    check_candidate_search(ARGS[1])
end
