#!/usr/bin/env julia
"""Independent zero-temperature gap and same-branch T -> 0 checks.

Run on quote-remote; only representative points are recomputed.
"""
module ZeroTemperatureGapStudy

using Dates
using FastGaussQuadrature: gausslegendre
using ForwardDiff
using JSON3
using LinearAlgebra
using SHA
using StaticArrays

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(ROOT, "src", "models", "Models.jl"))
end
const M = Main.Models
const HC = 197.327
const SOLVER = M.NLsolveGapSolver(ftol=1e-12, xtol=1e-12)
const QUADRATURE = (; thermo_quadrature_policy=:rs_reduced_adaptive,
    thermo_quadrature_rtol=1e-10, thermo_quadrature_atol=1e-13,
    thermo_quadrature_maxevals=200000)

function file_hash(path)
    return bytes2hex(sha256(read(joinpath(ROOT, path))))
end

function snapshot(model, state, T, mu; xi=0.0)
    densities = M.number_densities(model, state, T, mu; xi, QUADRATURE...)
    residual = M.gap_residual(model, state, T, mu; xi, QUADRATURE...)
    return (; state=collect(M.state_vector(state)),
        masses=collect(M.calculate_mass_vec(model, state.phi)),
        pressure=M.model_pressure(model, state, mu, T; xi, QUADRATURE...),
        net_densities=collect(densities.quark - densities.antiquark),
        residual_norm=maximum(abs, residual))
end

function differences(value, reference)
    return (; chiral_max_abs=maximum(abs, value.state[1:3] - reference.state[1:3]),
        mass_max_abs=maximum(abs, value.masses - reference.masses),
        pressure_abs=abs(value.pressure - reference.pressure),
        pressure_relative=abs(value.pressure - reference.pressure) / max(abs(reference.pressure), 1e-12),
        density_max_abs=maximum(abs, value.net_densities - reference.net_densities),
        density_relative=maximum(abs, value.net_densities - reference.net_densities) /
            max(maximum(abs, reference.net_densities), 1e-12))
end

function zero_chiral_response(model, state, mu_vec; xi=0.0)
    phi = state.phi
    pressure_phi_mu = function (y)
        st = M.MeanFieldState(SVector(y[1], y[2], y[3]); Phi=zero(y[1]), PhiBar=zero(y[1]))
        return -M.omega(model, st, 0.0, mu_vec .+ y[4]; xi)
    end
    hessian = ForwardDiff.hessian(pressure_phi_mu, SVector(phi[1], phi[2], phi[3], 0.0))
    chiral_jacobian = -hessian[1:3, 1:3]
    mixed = hessian[1:3, 4]
    dphi_dmu = chiral_jacobian \ mixed
    pressure_mu2 = hessian[4, 4] + dot(mixed, dphi_dmu)
    mass_jacobian = ForwardDiff.jacobian(p -> M.calculate_mass_vec(model, p), phi)
    return (; chiral_eigenvalues=collect(eigvals(Symmetric(chiral_jacobian))),
        jacobian_condition=cond(chiral_jacobian),
        dphi_dmu=collect(dphi_dmu), dmass_dmu=collect(mass_jacobian * dphi_dmu),
        pressure_mu2)
end

# Independent Gauss rules integrate E(q) and M/E(q), avoiding derivatives of
# the production closed-form sea and vacuum expressions.
function integrated_zero_reference(model, state, mu_vec; xi=0.0, nodes=128)
    xs, weights = gausslegendre(nodes)
    function integral(f, upper)
        upper == 0 && return 0.0
        half = upper / 2
        return half * sum(weights[i] * f(half * (xs[i] + 1)) for i in eachindex(xs))
    end
    masses = M.calculate_mass_vec(model, state.phi)
    cutoff = model.params.Λ_inv_fm
    nc = model.params.N_color
    angular = inv(sqrt(1 + xi))
    mass_gradient = zeros(3)
    vacuum, medium = 0.0, 0.0
    for i in 1:3
        mass, mu = masses[i], abs(mu_vec[i])
        qf = mu > mass ? sqrt(mu^2 - mass^2) : 0.0
        vacuum -= nc / pi^2 * integral(q -> q^2 * sqrt(q^2 + mass^2), cutoff)
        medium -= nc / pi^2 * angular * integral(q -> q^2 * (mu - sqrt(q^2 + mass^2)), qf)
        mass_gradient[i] = -nc / pi^2 * integral(q -> q^2 * mass / sqrt(q^2 + mass^2), cutoff) +
            nc / pi^2 * angular * integral(q -> q^2 * mass / sqrt(q^2 + mass^2), qf)
    end
    chiral_gradient = ForwardDiff.gradient(phi -> M.calculate_chiral(model, phi), state.phi)
    mass_jacobian = ForwardDiff.jacobian(phi -> M.calculate_mass_vec(model, phi), state.phi)
    residual = chiral_gradient + mass_jacobian' * mass_gradient
    return (; pressure=-(M.calculate_chiral(model, state.phi) + vacuum + medium),
        residual=collect(residual), residual_norm=maximum(abs, residual))
end

function zero_candidates(model, mu; xi=0.0)
    states = typeof(M.MeanFieldState(SVector(0.0, 0.0, 0.0)))[]
    failed_seeds = []
    for seed in M._zero_temperature_seeds()
        state = try
            M.solve_gap(model, 0.0, mu; xi, initial_guess=seed,
                solver=SOLVER, residual_norm_max=1e-10)
        catch err
            err isa InterruptException && rethrow()
            push!(failed_seeds, (; seed, error=sprint(showerror, err)))
            continue
        end
        any(other -> maximum(abs, other.phi - state.phi) < 1e-7, states) && continue
        push!(states, state)
    end
    return (; roots=[merge(snapshot(model, state, 0.0, mu; xi),
        zero_chiral_response(model, state, mu; xi)) for state in states], failed_seeds)
end

function run(output)
    BLAS.set_num_threads(1)
    model = M.create_model(:PNJL)
    points = ((mu=1.25, xi=0.0), (mu=1.8, xi=0.0), (mu=2.2, xi=0.0),
        (mu=2.2, xi=-0.3), (mu=2.2, xi=0.3))
    rows = []
    for point in points
        mu = SVector(point.mu, point.mu, point.mu)
        xi = point.xi
        default_zero = M.solve_gap(model, 0.0, mu; xi, solver=SOLVER, residual_norm_max=1e-10)
        zero_summary = snapshot(model, default_zero, 0.0, mu; xi)
        oracle = integrated_zero_reference(model, default_zero, mu; xi)
        oracle_fine = integrated_zero_reference(model, default_zero, mu; xi, nodes=192)
        # Start independently at 60 MeV, then follow the same stationary branch.
        finite_states = []
        seed = nothing
        for T_MeV in (60.0, 20.0, 5.0, 1.0, 0.1, 0.01)
            state = M.solve_gap(model, T_MeV / HC, mu; xi, initial_guess=seed,
                solver=SOLVER, residual_norm_max=1e-10, QUADRATURE...)
            context = M.PNJLChiBTaylorDiff.build_taylor_base_context(T_MeV / HC, mu;
                model, xi, base_state=M.state_vector(state), base_residual_tol=1e-13,
                polish_max_iterations=8, thermo_quadrature_kwargs=QUADRATURE)
            state = M.MeanFieldState(context.x0)
            seed = state
            push!(finite_states, (; T_MeV, snapshot(model, state, T_MeV / HC, mu; xi)...))
        end
        same_branch_zero = M.solve_gap(model, 0.0, mu; xi, initial_guess=seed,
            solver=SOLVER, residual_norm_max=1e-10)
        same_zero_summary = snapshot(model, same_branch_zero, 0.0, mu; xi)
        response = zero_chiral_response(model, same_branch_zero, mu; xi)
        comparison = [merge(row, (; versus_same_branch_zero=differences(row, same_zero_summary),
            versus_default_zero=differences(row, zero_summary))) for row in finite_states]
        # Directly re-solve nearby zero-temperature points for the implicit response.
        derivative_checks = map((1e-3, 3e-4, 1e-4)) do h
            lower = M.solve_gap(model, 0.0, mu .- h; xi, initial_guess=same_branch_zero,
                solver=SOLVER, residual_norm_max=1e-10)
            upper = M.solve_gap(model, 0.0, mu .+ h; xi, initial_guess=same_branch_zero,
                solver=SOLVER, residual_norm_max=1e-10)
            lo, hi = snapshot(model, lower, 0.0, mu .- h; xi), snapshot(model, upper, 0.0, mu .+ h; xi)
            return (; h, dmass_dmu=(hi.masses - lo.masses) / (2h),
                pressure_mu2=sum(hi.net_densities - lo.net_densities) / (2h))
        end
        push!(rows, (; mu_q_inv_fm=point.mu, xi, default_zero=zero_summary,
            same_branch_zero=same_zero_summary, response, oracle, oracle_fine,
            finite_temperature=comparison, derivative_checks,
            candidates=zero_candidates(model, mu; xi)))
        println("zero/finite comparison mu=", point.mu, " xi=", xi)
        flush(stdout)
    end
    files = ("src/models/pnjl_physics/PNJLModel.jl",
        "src/models/solver/runtime/ZeroTemperatureGapSolver.jl",
        "src/models/pnjl_physics/PNJLIntegrals.jl",
        "scripts/analysis/pnjl/zero_temperature_gap.jl")
    report = (; schema="pnjl_zero_temperature_gap_v1", created_utc=string(now(UTC)),
        execution=(; julia=string(VERSION), threads=Threads.nthreads(), blas_threads=BLAS.get_num_threads()),
        source_hashes=Dict(path => file_hash(path) for path in files),
        controls=(; temperatures_MeV=[0.0, 0.01, 0.1, 1.0, 5.0, 20.0, 60.0],
            quadrature=QUADRATURE, solver_ftol=SOLVER.ftol, polish_residual_tol=1e-13),
        units=(; chemical_potential="fm^-1", masses="fm^-1", condensates="fm^-3",
            pressure="fm^-4", number_density="fm^-3", pressure_mu2="fm^-2"), rows)
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        JSON3.pretty(io, report)
    end
    println("wrote ", output)
    return report
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: zero_temperature_gap.jl OUTPUT.json")
    ZeroTemperatureGapStudy.run(ARGS[1])
end
