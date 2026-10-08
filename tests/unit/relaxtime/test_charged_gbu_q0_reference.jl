using Test

isdefined(Main, :ChargedGBUContourScan) || include(joinpath(@__DIR__, "..", "..", "..",
    "scripts", "analysis", "relaxtime", "run_charged_gbu_contour_scan.jl"))
const Q0_REFERENCE_TEST = Main.ChargedGBUContourScan.Reference
# Other contract tests can reload Main.ChargedGBUResearchWorkflow. Extend the
# profile module captured by the reference under test, even after that reload.
const Q0_REFERENCE_P = Q0_REFERENCE_TEST.P
const Q0_REFERENCE_Y = Q0_REFERENCE_TEST.Y

struct Q0SyntheticProfile{K,T}
    kernel::K
    threshold::T
end

struct Q0SyntheticLandauProfile{K,T}
    kernel::K
    threshold::T
end

function Q0_REFERENCE_P.polarization(p::Q0SyntheticLandauProfile, z)
    x = real(z); k = p.kernel
    # A synthetic Landau cut with the exact external-static occupation zero.
    # No real zero or negative static inverse is involved in this example.
    rho = abs(x) < k.landau ? x * (x-k.shift) * (k.landau^2-x^2) : 0.
    return complex(0., rho) / (4k.bg.coupling[k.ch])
end

function Q0_REFERENCE_P.polarization(p::Q0SyntheticProfile, z)
    x = real(z)
    # One root on each analytic gap, then a pi -> 0 continuum on [2,4].
    f = abs(x) <= 2 ? complex((1 - x^2) / 3) :
        cis(-pi * clamp((4 - x) / 2, 0., 1.))
    return (1 - f) / (4p.kernel.bg.coupling[p.kernel.ch])
end

function q0_synthetic_profile(; shift=.2, q=0.)
    bg = (T=.8, coupling=Dict(:pi_plus=>1.))
    kernel = (q=q, bg=bg, ch=:pi_plus, shift=shift, threshold=2., landau=.2,
        split=40., edges=[-40., -4., -2., -.2, 0., .2, 2., 4., 40.])
    return Q0SyntheticProfile(kernel, (U=2., H=.25, positive=1., negative=-1.))
end

@testset "q0 reference internal-lambda coordinates and units" begin
    for shift in (-.3, 0., .4), lambda in (0., 1., 3.)
        @test Q0_REFERENCE_TEST.reference_coordinate(lambda-shift, 0., shift) ≈ lambda atol=1e-14
    end
    @test Q0_REFERENCE_TEST.reference_coordinate(2.6, 2., .4) ≈ sqrt(5.)
    @test Q0_REFERENCE_TEST.reference_coordinate(.1, 1., .2) === nothing
    @test Q0_REFERENCE_TEST.reference_coordinate(-2., 1., 0.) === nothing
    @test Q0_REFERENCE_TEST.boosted_energy(1., 2., .3) ≈ sqrt(5.)-.3
    @test_throws ArgumentError Q0_REFERENCE_TEST.reference_coordinate(1., -1., 0.)
    @test_throws ArgumentError Q0_REFERENCE_TEST.reference_coordinate(Inf, 1., 0.)
end

@testset "Independent q0 gap count and unaltered positive-frequency q=0 limit" begin
    p = q0_synthetic_profile(); ref = Q0_REFERENCE_TEST.reference(p)
    @test ref.rootsets[1] ≈ [-1.] atol=1e-12
    @test ref.rootsets[2] ≈ [1.] atol=1e-12
    actual = Q0_REFERENCE_TEST.shell(ref, 0.)
    ordinary = Q0_REFERENCE_Y.shell(p)
    @test actual.roots ≈ ordinary.roots
    @test actual.threshold_phase == ordinary.threshold_phase == Float64(pi)
    @test actual.density == ordinary.density == 0.
    @test Q0_REFERENCE_TEST.cut_integral(ref, 1.8, 24., 0., .8) ≈
        Q0_REFERENCE_Y.cut_integral(p, 1.8, 24., .8) atol=1e-12
    @test_throws ArgumentError Q0_REFERENCE_TEST.reference(q0_synthetic_profile(q=1.))
    unsafe = Q0_REFERENCE_TEST.reference(q0_synthetic_profile(shift=1.2))
    @test_throws ArgumentError Q0_REFERENCE_TEST.shell(unsafe, .1)
end

@testset "Boosted discrete root, signed GBU continuum and boundary bookkeeping" begin
    p = q0_synthetic_profile(); ref = Q0_REFERENCE_TEST.reference(p); q = 1.3
    s = Q0_REFERENCE_TEST.shell(ref, q; nodes=96)
    T = p.kernel.bg.T; shift = p.kernel.shift
    lower = hypot(2., q)-shift; upper = hypot(4., q)-shift
    g(w) = inv(expm1(w/T))
    derivative = real(Q0_REFERENCE_Y.I.C.mapped_integral([lower, upper]; nodes=96) do w
        lambda0 = sqrt((w+shift)^2-q^2)
        delta = pi*(4-lambda0)/2
        # d[delta-sin(2delta)/2]/dw carries this coordinate derivative.
        -g(w)*sin(delta)^2*(w+shift)/lambda0
    end)
    measure = q^2/(2pi^2)
    @test s.roots ≈ [hypot(1.,q)-shift]
    @test s.root_count == 1 && s.negative_root_count == 1
    @test s.landau == 0.
    @test s.pair ≈ measure*derivative atol=1e-12
    @test s.pair < 0 < s.density
    @test s.bound ≈ measure*g(only(s.roots)) atol=1e-12
    @test s.density ≈ s.bound+s.landau+s.pair atol=1e-12
    @test Q0_REFERENCE_TEST.cut_phase(ref, q-shift-1e-7, q) == 0.
    @test s.route == "q0_lambda_reference" && !s.production_authorized
    @test_throws ArgumentError Q0_REFERENCE_TEST.shell(ref, -1.)
    @test_throws ArgumentError Q0_REFERENCE_TEST.shell(ref, q; upper=lower)
end

@testset "Lambda extrapolation does not automatically preserve the Bose endpoint" begin
    base = q0_synthetic_profile(shift=.4)
    p = Q0SyntheticLandauProfile(merge(base.kernel, (landau=.8,)), base.threshold)
    ref = Q0_REFERENCE_TEST.reference(p)
    @test iszero(imag(Q0_REFERENCE_P.inverse(p, p.kernel.shift)))
    @test all(isempty, ref.rootsets)
    @test Q0_REFERENCE_Y.shell(p).static_inverse == 1.
    @test !iszero(Q0_REFERENCE_TEST.cut_phase(ref, 0., .1))
    @test_throws ArgumentError Q0_REFERENCE_TEST.shell(ref, .1; endpoint_policy="strict_zero_limit")
    finite = Q0_REFERENCE_TEST.shell(ref, .1)
    @test isfinite(finite.density) && finite.endpoint_warning
    @test finite.static_inverse == 1. && abs(finite.static_inverse_imag) > 1e-8
    @test finite.warning_code == "q0_extrapolated_static_imaginary"
    @test finite.omega_lower_inv_fm == 1e-5 && finite.endpoint_policy == "finite_window"
    @test Q0_REFERENCE_TEST.cut_phase(ref, 0., .5) == 0.
end

@testset "Stable finite-window derivative and explicit boundary identity" begin
    r = Q0_REFERENCE_TEST
    for phase in (-.3, 1e-8, .3, pi), lo in (1e-3, 1e-5, 1e-7)
        s = r.finite_window_integral(w -> phase, lo, .2, .7)
        @test s.derivative == 0.
        @test s.lower_boundary ≈ inv(expm1(lo/.7))*Q0_REFERENCE_Y.O.gbu_weight(phase)/pi
    end
    for lo in (1e-3, 1e-5, 1e-7)
        phase(w) = .3 + .2w
        s = r.finite_window_integral(phase, lo, .2, .7; nodes=128)
        # Independent analytic dW/dw; do not differentiate the new implementation.
        edges = exp.(range(log(lo), log(.2); length=33)); edges[1]=lo; edges[end]=.2
        exact = real(Q0_REFERENCE_Y.I.C.mapped_integral(edges; nodes=128) do w
            inv(expm1(w/.7))*.4sin(phase(w))^2/pi
        end)
        @test s.derivative ≈ exact rtol=1e-10 atol=1e-12
        bulk = real(Q0_REFERENCE_Y.I.C.mapped_integral(edges; nodes=128) do w
            g=inv(expm1(w/.7)); g*(1+g)/.7*Q0_REFERENCE_Y.O.gbu_weight(phase(w))/pi
        end)
        @test s.derivative ≈ bulk+s.upper_boundary-s.lower_boundary atol=2e-9
    end
    for anchor in (-pi, -.3, -1e-5, 0., 1e-8, .2, pi), h in (-1e-10, 0., 1e-10)
        delta = anchor+h
        exact = setprecision(256) do
            d, a = BigFloat(delta), BigFloat(anchor)
            Float64((d-sin(2d)/2)-(a-sin(2a)/2))
        end
        @test r.weight_difference(delta, anchor) ≈ exact rtol=1e-6 atol=1e-28
    end
    @test_throws ArgumentError r.finite_window_integral(identity, 0., 1., .7)
    @test_throws ArgumentError r.finite_window_integral(identity, .1, .1, .7)
    @test_throws ArgumentError r.shell(r.reference(q0_synthetic_profile()), 1.; omega_lower_inv_fm=0.)
    @test_throws ArgumentError r.shell(r.reference(q0_synthetic_profile()), 1.; endpoint_policy="ignore")
end

@testset "Safe legacy parity and retained root/geometry safeguards" begin
    r = Q0_REFERENCE_TEST
    for shift in (-.3, 0., .2), q in (.1, 1.3, 3.)
        ref = r.reference(q0_synthetic_profile(;shift=shift))
        finite = r.shell(ref, q; nodes=96)
        strict = r.shell(ref, q; nodes=96, endpoint_policy="strict_zero_limit")
        @test finite.density ≈ strict.density rtol=1e-9 atol=1e-12
        @test finite.bound == strict.bound
        @test !finite.endpoint_warning
    end
    ref = r.reference(q0_synthetic_profile(shift=.4))
    bad_roots = r.Reference(ref.profile, (Float64[], [.1]))
    @test_throws ArgumentError r.shell(bad_roots, .1)
    extra_roots = r.Reference(ref.profile, (Float64[], [.6, 1.]))
    @test_throws ArgumentError r.shell(extra_roots, .1)
end
