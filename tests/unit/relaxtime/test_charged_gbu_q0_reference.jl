using Test

isdefined(Main, :ChargedGBUContourScan) || include(joinpath(@__DIR__, "..", "..", "..",
    "scripts", "analysis", "relaxtime", "run_charged_gbu_contour_scan.jl"))
const Q0_REFERENCE_TEST = Main.ChargedGBUContourScan.Reference
const Q0_REFERENCE_P = Main.ChargedGBUResearchWorkflow.P
const Q0_REFERENCE_Y = Main.ChargedGBUResearchWorkflow.Y

struct Q0SyntheticProfile{K,T}
    kernel::K
    threshold::T
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
