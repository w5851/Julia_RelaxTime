#=
Isolated full-RS A angular-reduction feasibility and cost probe.

Run with Julia 1.12.5 on quote-remote:
    julia --project=. benchmark/relaxtime/bench_full_rs_a_reduction.jl OUTPUT_DIR

The prototype below does not replace a production entrypoint. Both quadratures
retain the ordinary on-shell denominator, and the vacuum contribution is read
from the existing production implementation.
=#
using BenchmarkTools
using JSON
using Printf
using SHA
using Statistics

const PROJECT_ROOT_A_PROBE = normpath(joinpath(@__DIR__, "..", ".."))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, joinpath(PROJECT_ROOT_A_PROBE, "src", "relaxtime", "RelaxTime.jl"))
end
using Main.RelaxTime.GaussLegendre: gauleg
using Main.RelaxTime.OneLoopIntegrals: const_integral_term_A
using Main.RelaxTime.OneLoopIntegralsCorrection: A_aniso
using Main.PNJLQuarkDistributions_Aniso: distribution_aniso

function semi_infinite_grid_a_probe(n::Int, scale::Float64)
    x, w = gauleg(0.0, 1.0, n)
    return scale .* x ./ (1 .- x), scale .* w ./ (1 .- x).^2
end

"""q^2 times the angular weight on [-u,u], after q=p*sqrt(1+xi*cos(theta)^2)."""
function radial_weight_a_probe(q::Float64, m::Float64, xi::Float64, u::Float64=1.0)
    isfinite(xi) && xi > -1.0 || throw(ArgumentError("xi must be finite and greater than -1"))
    q >= 0.0 && m >= 0.0 || throw(ArgumentError("q and m must be nonnegative"))
    0.0 <= u <= 1.0 || throw(ArgumentError("u must lie in [0,1]"))
    (q == 0.0 || u == 0.0) && return 0.0
    e = hypot(q, m)
    xi == 0.0 && return 2.0 * q^2 * u / e
    if xi > 0.0
        d = hypot(q, m * sqrt(1.0 + xi * u^2))
        z = q * sqrt(xi) * u / d
        return 2.0 * q^2 * u / d * (z == 0.0 ? 1.0 : atan(z) / z)
    end
    # Equivalent to atanh(q*sqrt(-xi)*u/sqrt(q^2+m^2*(1+xi*u^2))).
    # asinh avoids rounding that atanh argument to one as xi approaches -1.
    d = e * sqrt(1.0 + xi * u^2)
    z = q * sqrt(-xi) * u / d
    return 2.0 * q^2 * u / d * (z == 0.0 ? 1.0 : asinh(z) / z)
end

function thermal_sum_a_probe(q::Float64, case)
    return distribution_aniso(:quark, q, case.m, case.mu, case.T, case.phi, case.phibar, 0.0, 0.0) +
        distribution_aniso(:antiquark, q, case.m, case.mu, case.T, case.phi, case.phibar, 0.0, 0.0)
end

function thermal_1d_a_probe(case, q, w)
    integral = 0.0
    for i in eachindex(q, w)
        integral += w[i] * radial_weight_a_probe(q[i], case.m, case.xi) * thermal_sum_a_probe(q[i], case)
    end
    return 2.0 * integral
end

function thermal_2d_a_probe(case, p, wp, c, wc)
    integral = 0.0
    for i in eachindex(p, wp)
        prefactor = p[i]^2 / hypot(p[i], case.m)
        angular = 0.0
        for j in eachindex(c, wc)
            f = distribution_aniso(:quark, p[i], case.m, case.mu, case.T, case.phi, case.phibar, case.xi, c[j]) +
                distribution_aniso(:antiquark, p[i], case.m, case.mu, case.T, case.phi, case.phibar, case.xi, c[j])
            angular += wc[j] * prefactor * f
        end
        integral += wp[i] * angular
    end
    return 2.0 * integral
end

function thermal_1d_finite_a_probe(case, pmax::Float64, n::Int)
    # The original spherical p cutoff becomes a q-dependent angular domain.
    edge = sqrt(1.0 + case.xi)
    qbreak, qmax = pmax * min(1.0, edge), pmax * max(1.0, edge)
    q, w = gauleg(0.0, qbreak, n)
    total = thermal_1d_a_probe(case, q, w)
    qmax == qbreak && return total
    t, wt = gauleg(0.0, 1.0, n)
    width = qmax - qbreak
    for i in eachindex(t, wt)
        # Square the distance from the square-root endpoint of the angular domain.
        qi = case.xi > 0.0 ? qbreak + width * t[i]^2 : qmax - width * t[i]^2
        u2 = ((qi / pmax)^2 - 1.0) / case.xi
        -32eps() <= u2 <= 1.0 + 32eps() || error("finite-cutoff angular domain out of bounds: $u2")
        u = sqrt(clamp(u2, 0.0, 1.0))
        weight = case.xi > 0.0 ?
            radial_weight_a_probe(qi, case.m, case.xi) - radial_weight_a_probe(qi, case.m, case.xi, u) :
            radial_weight_a_probe(qi, case.m, case.xi, u)
        total += 2.0 * wt[i] * (2.0 * width * t[i]) * weight * thermal_sum_a_probe(qi, case)
    end
    return total
end

function angular_bigfloat_reference_a_probe(q::Float64, m::Float64, xi::Float64)
    return setprecision(256) do
        qb, mb, xb = BigFloat(q), BigFloat(m), BigFloat(xi)
        qb == 0 && return 0.0
        if xb == 0
            return Float64(2qb^2 / sqrt(qb^2 + mb^2))
        end
        d = sqrt(qb^2 + mb^2 * (1 + xb))
        value = xb > 0 ? 2qb / sqrt(xb) * atan(qb * sqrt(xb) / d) :
            2qb / sqrt(-xb) * atanh(qb * sqrt(-xb) / d)
        return Float64(value)
    end
end

function write_probe_json(path, value)
    open(path, "w") do io
        JSON.print(io, value, 2)
        println(io)
    end
end

function main_a_probe(outdir::String)
    ispath(outdir) && error("Output directory exists; preserve the previous probe: $outdir")
    mkpath(outdir)
    failures = String[]
    bounds = NamedTuple[]
    for xi in (-1.0 + eps(), -0.999999, -0.9, -1e-12, 0.0, 1e-12, 0.5, 100.0),
        m in (0.0, 0.3, 3.0), q in (0.0, 1e-12, 0.1, 3.0, 1e5)
        value = radial_weight_a_probe(q, m, xi)
        reference = angular_bigfloat_reference_a_probe(q, m, xi)
        passed = isfinite(value) && isapprox(value, reference; rtol=2e-13, atol=1e-25)
        push!(bounds, (;xi, m, q, value, reference, passed))
        passed || push!(failures, "angular weight xi=$xi m=$m q=$q")
    end
    for xi in (-2.0, -1.0, NaN, Inf)
        rejected = try
            radial_weight_a_probe(1.0, 0.3, xi)
            false
        catch error
            error isa ArgumentError
        end
        rejected || push!(failures, "invalid xi accepted: $xi")
    end
    write_probe_json(joinpath(outdir, "angular_limits.json"), bounds)

    backgrounds = [
        (name="massless", m=0.0, mu=0.1, T=0.25, phi=0.2, phibar=0.3),
        (name="light", m=0.3, mu=0.5, T=1.0, phi=0.5, phibar=0.4),
        (name="heavy", m=2.0, mu=1.4, T=0.65, phi=0.3, phibar=0.2),
    ]
    q256, w256 = semi_infinite_grid_a_probe(256, 5.0)
    q512, w512 = semi_infinite_grid_a_probe(512, 5.0)
    c128, wc128 = gauleg(-1.0, 1.0, 128)
    c256, wc256 = gauleg(-1.0, 1.0, 256)
    equivalence = NamedTuple[]
    for background in backgrounds, xi in (-0.9, -0.5, -0.22, -1e-10, 0.0, 1e-10, 0.37, 0.5, 5.0)
        case = merge(background, (;xi))
        one = thermal_1d_a_probe(case, q512, w512)
        one_coarse = thermal_1d_a_probe(case, q256, w256)
        two = thermal_2d_a_probe(case, q512, w512, c256, wc256)
        two_coarse = thermal_2d_a_probe(case, q256, w256, c128, wc128)
        passed = all(isfinite, (one, two, one_coarse, two_coarse)) &&
            isapprox(one, two; rtol=1e-9, atol=1e-11) &&
            isapprox(one, one_coarse; rtol=1e-9, atol=1e-11) &&
            isapprox(two, two_coarse; rtol=1e-9, atol=1e-11)
        push!(equivalence, (;domain="semi_infinite", name=case.name, xi, pmax=nothing,
            thermal_1d=one, thermal_2d=two, one_refinement=one-one_coarse,
            two_refinement=two-two_coarse, relative_error=abs(one/two-1), passed))
        passed || push!(failures, "semi-infinite $(case.name) xi=$xi")
        println("EQUIVALENCE ", case.name, " xi=", xi, " passed=", passed); flush(stdout)
    end
    for background in backgrounds, xi in (-0.9, -0.5, 0.0, 0.5, 3.0), pmax in (0.5, 3.0, 25.0)
        case = merge(background, (;xi))
        p, wp = gauleg(0.0, pmax, 256)
        pc, wpc = gauleg(0.0, pmax, 128)
        one = thermal_1d_finite_a_probe(case, pmax, 256)
        one_coarse = thermal_1d_finite_a_probe(case, pmax, 128)
        two = thermal_2d_a_probe(case, p, wp, c256, wc256)
        two_coarse = thermal_2d_a_probe(case, pc, wpc, c128, wc128)
        passed = all(isfinite, (one, two, one_coarse, two_coarse)) &&
            isapprox(one, two; rtol=1e-9, atol=1e-11) &&
            isapprox(one, one_coarse; rtol=1e-9, atol=1e-11) &&
            isapprox(two, two_coarse; rtol=1e-9, atol=1e-11)
        push!(equivalence, (;domain="finite_spherical_p", name=case.name, xi, pmax,
            thermal_1d=one, thermal_2d=two, one_refinement=one-one_coarse,
            two_refinement=two-two_coarse, relative_error=abs(one/two-1), passed))
        passed || push!(failures, "finite $(case.name) xi=$xi pmax=$pmax")
    end
    write_probe_json(joinpath(outdir, "equivalence.json"), equivalence)

    q32, w32 = semi_infinite_grid_a_probe(32, 5.0)
    c16, wc16 = gauleg(-1.0, 1.0, 16)
    timings = NamedTuple[]
    for background in backgrounds[2:3], xi in (-0.5, -0.22, 0.37, 0.5)
        case = merge(background, (;xi))
        production_call = () -> A_aniso(case.m, case.mu, case.T, case.phi, case.phibar,
            case.xi, q32, w32, c16, wc16)
        reduced_call = () -> thermal_1d_a_probe(case, q32, w32) - 4.0 * const_integral_term_A(case.m)
        production_value, reduced_value = production_call(), reduced_call()
        reference = thermal_1d_a_probe(case, q512, w512) - 4.0 * const_integral_term_A(case.m)
        old_trial = @benchmark $production_call() samples=200 evals=10 seconds=0.25
        new_trial = @benchmark $reduced_call() samples=200 evals=10 seconds=0.25
        old_time, new_time = median(old_trial).time * 1e-9, median(new_trial).time * 1e-9
        push!(timings, (;name=case.name, xi, production_value, reduced_value, reference,
            production_absolute_error=abs(production_value-reference),
            reduced_absolute_error=abs(reduced_value-reference),
            production_seconds=old_time, reduced_seconds=new_time, speedup=old_time/new_time,
            production_allocated_bytes=old_trial.memory, reduced_allocated_bytes=new_trial.memory))
    end
    write_probe_json(joinpath(outdir, "timings.json"), timings)
    summary = (;complete=true, prototype_only=true, production_entry_changed=false,
        numerical_host=get(ENV,"COMPUTERNAME","unknown"), julia_version=string(VERSION),
        julia_threads=Threads.nthreads(), angular_checks=length(bounds),
        integral_checks=length(equivalence), timing_cases=length(timings),
        passed=isempty(failures), failures,
        maximum_thermal_relative_error=maximum(r.relative_error for r in equivalence),
        median_kernel_speedup=median(r.speedup for r in timings),
        maximum_existing_A_seconds=maximum(r.production_seconds for r in timings),
        whole_transport_speedup_measured=false,
        probe_sha256=bytes2hex(sha256(read(@__FILE__))))
    write_probe_json(joinpath(outdir, "summary.json"), summary)
    println(JSON.json(summary)); flush(stdout)
    return isempty(failures)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: bench_full_rs_a_reduction.jl OUTPUT_DIR")
    main_a_probe(abspath(ARGS[1])) || exit(1)
end
