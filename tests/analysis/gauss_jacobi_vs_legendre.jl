#!/usr/bin/env julia

using FastGaussQuadrature
using Printf

"""
Compare convergence of Gauss–Legendre vs Gauss–Jacobi on an integrable
endpoint singularity of the form (s - a)^(-1/2).

Integrand: f(s) = (s - a)^(-1/2) * (1 + c * (s - a))
Exact: ∫ f(s) ds from s=a to b = 2√(b-a) + (2/3) * c * (b-a)^(3/2)
"""
function run_experiment(; a=25.0, b=26.0, c=0.3, orders=(2,4,8,16,32,64))
    Δ = b - a
    exact = 2 * sqrt(Δ) + (2/3) * c * Δ^(3/2)
    @printf("Interval [%.2f, %.2f], Δ=%.3f, c=%.3f, exact=%.10f\n", a, b, Δ, c, exact)
    println("n\tmethod\t\tapprox\t\trel_error")

    # Gauss–Legendre (no special weighting)
    for n in orders
        x, w = gausslegendre(n)
        s = @. (Δ / 2) * (x + 1) + a
        jac = Δ / 2
        approx = sum(@. w * ((s - a)^(-0.5) * (1 + c * (s - a)))) * jac
        rel_err = abs(approx - exact) / exact
        @printf("%d\tLegendre\t%.10f\t%.3e\n", n, approx, rel_err)
    end

    # Gauss–Jacobi with weight (1 + x)^(-1/2) mapped to (s - a)^(-1/2)
    for n in orders
        x, w = gaussjacobi(n, 0.0, -0.5)  # weight = (1 + x)^(-1/2)
        s = @. (Δ / 2) * (x + 1) + a
        # Integral = sqrt(Δ/2) * Σ w_i * g(s_i)
        prefactor = sqrt(Δ / 2)
        approx = prefactor * sum(@. w * (1 + c * (s - a)))
        rel_err = abs(approx - exact) / exact
        @printf("%d\tJacobi\t\t%.10f\t%.3e\n", n, approx, rel_err)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_experiment()
end
