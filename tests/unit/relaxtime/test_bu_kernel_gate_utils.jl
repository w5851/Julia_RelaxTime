using Test

const _BU_GATE_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "bu_kernel_gate_utils.jl",
))
if !isdefined(Main, :BUKernelGateUtils)
    Base.include(Main, _BU_GATE_PATH)
end
using Main.BUKernelGateUtils: bose_occupation, finite_window_bu_identity,
                              nonuniform_three_point_derivative

function _midpoint_rule(a::Float64, b::Float64, n::Int)
    h = (b - a) / n
    return (collect(range(a + h / 2.0; step=h, length=n)), fill(h, n))
end

@testset "finite-window BU identity gate" begin
    T = 0.31
    mu = 0.07
    a = 0.12
    b = 2.4

    @test bose_occupation(0.8, mu, T) > 0.0
    @test_throws ArgumentError bose_occupation(0.8, mu, 0.0)
    @test_throws ArgumentError bose_occupation(mu, mu, T)

    x_irregular = [0.1, 0.27, 0.8, 1.1, 2.3]
    @test nonuniform_three_point_derivative(x_irregular, x_irregular .^ 2) ≈ 2.0 .* x_irregular atol=1e-13
    @test nonuniform_three_point_derivative([0.2, 0.9], [1.0, 2.4]) == [2.0, 2.0]
    @test_throws ArgumentError nonuniform_three_point_derivative([0.1, 0.1], [1.0, 2.0])

    omega, weights = _midpoint_rule(a, b, 1024)

    constant_F = fill(0.8, length(omega))
    constant_result = finite_window_bu_identity(
        omega, weights, constant_F, zeros(length(omega));
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=0.8, F_max=0.8,
    )
    @test constant_result.boundary_lower < 0.0
    @test constant_result.boundary_upper > 0.0
    @test constant_result.byparts_total == constant_result.byparts_bulk + constant_result.boundary
    @test constant_result.derivative == 0.0
    @test constant_result.closure_abs < 2e-4

    linear_F(x) = 0.2 + 0.6 * x
    linear_result = finite_window_bu_identity(
        omega, weights, linear_F.(omega), fill(0.6, length(omega));
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=linear_F(a), F_max=linear_F(b),
    )
    @test linear_result.closure_abs < 2e-4

    delta(x) = 0.3 + 0.4 * x + 0.05 * sin(1.7 * x)
    ddelta(x) = 0.4 + 0.085 * cos(1.7 * x)
    current = finite_window_bu_identity(
        omega, weights, delta.(omega), ddelta.(omega);
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=delta(a), F_max=delta(b),
    )
    gbu_F(x) = delta(x) - 0.5 * sin(2.0 * delta(x))
    gbu_dF(x) = 2.0 * sin(delta(x))^2 * ddelta(x)
    gbu = finite_window_bu_identity(
        omega, weights, gbu_F.(omega), gbu_dF.(omega);
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=gbu_F(a), F_max=gbu_F(b),
    )
    @test current.closure_abs < 2e-4
    @test gbu.closure_abs < 2e-4

    omega_coarse, weights_coarse = _midpoint_rule(a, b, 32)
    coarse = finite_window_bu_identity(
        omega_coarse, weights_coarse, gbu_F.(omega_coarse), gbu_dF.(omega_coarse);
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=gbu_F(a), F_max=gbu_F(b),
    )
    @test gbu.closure_abs < coarse.closure_abs

    @test_throws ArgumentError finite_window_bu_identity(
        omega, weights, constant_F, zeros(length(omega));
        T=0.0, mu=mu, omega_min=a, omega_max=b, F_min=0.8, F_max=0.8,
    )
    @test_throws ArgumentError finite_window_bu_identity(
        omega, weights, constant_F, zeros(length(omega));
        T=T, mu=a, omega_min=a, omega_max=b, F_min=0.8, F_max=0.8,
    )
    @test_throws ArgumentError finite_window_bu_identity(
        omega, weights, constant_F, zeros(length(omega));
        T=T, mu=mu, omega_min=b, omega_max=a, F_min=0.8, F_max=0.8,
    )
    @test_throws ArgumentError finite_window_bu_identity(
        omega, weights[1:end-1], constant_F, zeros(length(omega));
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=0.8, F_max=0.8,
    )
    bad_F = copy(constant_F)
    bad_F[1] = NaN
    @test_throws ArgumentError finite_window_bu_identity(
        omega, weights, bad_F, zeros(length(omega));
        T=T, mu=mu, omega_min=a, omega_max=b, F_min=0.8, F_max=0.8,
    )
end
