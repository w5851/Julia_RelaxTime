#!/usr/bin/env julia
"""Compare exact zero, near-zero and ordinary-temperature PNJL costs.

Run on quote-remote after zero_temperature_gap.jl. Compilation is excluded from
steady-state timing; common-state kernels and root-search policies are separate.
"""
module ZeroTemperatureGapBenchmark

using BenchmarkTools
using Dates
using ForwardDiff
using JSON3
using LinearAlgebra
using SHA
using StaticArrays
using Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(ROOT, "src", "models", "Models.jl"))
end
const M = Main.Models
const HC = 197.327
const QUADRATURE = (; thermo_quadrature_policy=:rs_reduced_adaptive,
    thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-11,
    thermo_quadrature_maxevals=200000)
const SOLVER = M.NLsolveGapSolver(ftol=1e-11, xtol=1e-12)

function timing(f; samples=30, seconds=1.0)
    f()
    # Record at least the requested sample count for expensive catalog searches.
    # The previous fixed one-second budget could leave only one observation.
    elapsed = @elapsed f()
    budget = max(seconds, 2samples * elapsed)
    trial = @benchmark $f() samples=samples evals=1 seconds=budget
    estimate = median(trial)
    return (; median_us=estimate.time / 1000, minimum_us=minimum(trial.times) / 1000,
        p90_us=quantile(trial.times, 0.9) / 1000, samples=length(trial.times),
        memory_bytes=estimate.memory, allocations=estimate.allocs)
end

function pressure_chemical4(masses, mu, T, xi)
    delta = M.PNJLChiBTaylorDiff.taylor_variable(0.0, 4)
    value = -M.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(masses,
        0.1, 0.1, mu .+ delta, T, xi; rtol=1e-8, atol=1e-11, maxevals=200000)
    return ntuple(i -> M.PNJLChiBTaylorDiff.nth_derivative_from_series(value, i - 1), 5)
end

function solver_work(model, T, mu, xi, seed; zero_temperature=false)
    if zero_temperature
        potential = M._zero_temperature_gap_potential(model, mu; xi)
        residual! = M.ZeroTemperatureGapResidual(potential)
        initial = collect(seed[1:3])
    else
        kwargs = (; xi, QUADRATURE...)
        residual! = M.GapResidualInplace{typeof(model), Float64, typeof(mu), typeof(kwargs), 5}(
            model, T, mu, kwargs)
        initial = collect(seed)
    end
    result = M.nlsolve(residual!, initial; autodiff=:forward, method=SOLVER.method,
        ftol=SOLVER.ftol, xtol=SOLVER.xtol, iterations=1000)
    result.f_converged || error("timed single-seed solve did not converge")
    return (; iterations=result.iterations, residual_calls=result.f_calls,
        jacobian_calls=result.g_calls, residual_norm=result.residual_norm)
end

function run(output; samples=30)
    BLAS.set_num_threads(1)
    model = M.create_model(:PNJL)
    temperatures = (0.0, 0.01, 0.1, 1.0, 5.0, 60.0, 150.0, 250.0)
    masses = SVector(1.5, 1.5, 2.5)
    mu = SVector(1.8, 1.8, 1.8)
    kernels = []
    for T_MeV in temperatures
        T = T_MeV / HC
        pressure() = -M.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
            masses, 0.1, 0.1, mu, T, 0.0; rtol=1e-8, atol=1e-11, maxevals=200000)
        response() = pressure_chemical4(masses, mu, T, 0.0)
        push!(kernels, (; T_MeV, pressure=pressure(), pressure_mu_0_to_4=response(),
            pressure_timing=timing(pressure; samples),
            chemical4_timing=timing(response; samples)))
        println("kernel T=", T_MeV)
        flush(stdout)
    end

    # Dense matter retains the same branch from 60 MeV down to zero.
    # Ordinary-temperature measurements use independently converged states.
    mu = SVector(2.2, 2.2, 2.2)
    cases = []
    previous = nothing
    for T_MeV in (250.0, 150.0, 60.0, 5.0, 1.0, 0.1, 0.01, 0.0)
        T = T_MeV / HC
        state = M.solve_gap(model, T, mu; initial_guess=previous,
            solver=SOLVER, residual_norm_max=1e-9, QUADRATURE...)
        previous = state
        exact = M.state_vector(state)
        # A reproducible nearby seed, not the converged root itself. Identical
        # relative displacement for each temperature and for both dimensions.
        seed = collect(exact)
        seed[1:3] .*= 1.02
        T > 0 && (seed[4:5] .*= 1.02)
        single() = M.solve_gap(model, T, mu; initial_guess=seed,
            solver=SOLVER, residual_norm_max=1e-9, QUADRATURE...)
        selected() = M.solve_gap(model, T, mu;
            solver=SOLVER, residual_norm_max=1e-9, QUADRATURE...)
        pressure() = M.model_pressure(model, state, mu, T; QUADRATURE...)
        density() = M.number_densities(model, state, T, mu; QUADRATURE...)
        if T == 0
            potential = M._zero_temperature_gap_potential(model, mu)
            residual = () -> ForwardDiff.gradient(potential, state.phi)
            jacobian = () -> ForwardDiff.hessian(potential, state.phi)
        else
            residual = () -> M.gap_residual(model, exact, T, mu; QUADRATURE...)
            jacobian = () -> ForwardDiff.jacobian(
                x -> M.gap_residual(model, x, T, mu; QUADRATURE...), exact)
        end
        chosen = selected()
        single_state = single()
        numerical_drift = maximum(abs, single_state.phi - state.phi)
        chosen_drift = maximum(abs, chosen.phi - state.phi)
        push!(cases, (; T_MeV, mu_q_inv_fm=mu[1], state=collect(exact),
            seed, masses=collect(M.calculate_mass_vec(model, state.phi)),
            pressure=pressure(), residual_norm=maximum(abs, residual()),
            jacobian_condition=cond(jacobian()), numerical_drift, chosen_drift,
            work=solver_work(model, T, mu, 0.0, seed; zero_temperature=T == 0),
            pressure_timing=timing(pressure; samples), density_timing=timing(density; samples),
            residual_timing=timing(residual; samples), jacobian_timing=timing(jacobian; samples),
            single_seed_timing=timing(single; samples), default_selection_timing=timing(selected; samples)))
        println("solver T=", T_MeV)
        flush(stdout)
    end
    sources = ("src/models/pnjl_physics/PNJLModel.jl",
        "src/models/solver/runtime/ZeroTemperatureGapSolver.jl",
        "src/models/pnjl_physics/PNJLIntegrals.jl",
        "benchmark/pnjl/bench_zero_temperature_gap.jl")
    report = (; schema="pnjl_zero_temperature_performance_v1", created_utc=string(now(UTC)),
        execution=(; julia=string(VERSION), threads=Threads.nthreads(), blas_threads=BLAS.get_num_threads(),
            cpu=Sys.cpu_info()[1].model, blas=string(BLAS.get_config()), word_size=Sys.WORD_SIZE),
        source_hashes=Dict(path => bytes2hex(sha256(read(joinpath(ROOT, path)))) for path in sources),
        controls=(; samples, evals=1, quadrature=QUADRATURE, solver_ftol=SOLVER.ftol,
            seed_relative_displacement=0.02, timed_scope="warmed calls, excludes Julia startup and JIT",
            cpu_affinity_note=get(ENV, "PNJL_BENCHMARK_CPU_AFFINITY", "unrestricted")),
        kernels, solvers=cases)
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        JSON3.pretty(io, report)
    end
    println("wrote ", output)
    return report
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: bench_zero_temperature_gap.jl OUTPUT.json")
    ZeroTemperatureGapBenchmark.run(ARGS[1])
end
