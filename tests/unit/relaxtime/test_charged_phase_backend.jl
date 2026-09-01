"""Pure synthetic-path tests for the strict charged phase/BU backend."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _RELAXTIME_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.RelaxTime.ChargedPhaseBackend: StrictChargedPhaseSpec,
                                           strict_retarded_phase,
                                           strict_phase_profile,
                                           strict_phase_gate,
                                           strict_mott_gate,
                                           strict_charged_bu_density,
                                           strict_charged_rpa_bu_density,
                                           strict_density_convergence_gate

@testset "Strict charged phase algebra" begin
    δ = 0.7
    inverse = cis(-δ)
    @test strict_retarded_phase(inverse) ≈ δ atol=1e-14
    @test strict_retarded_phase(conj(inverse); phase_object=:propagator, phase_sign=1) ≈ δ atol=1e-14
    @test_throws ArgumentError StrictChargedPhaseSpec(phase_sign=0)
    @test_throws ArgumentError strict_retarded_phase(inverse; phase_object=:unknown)
end

@testset "Strict phase profile has an explicit endpoint contract" begin
    ω = collect(range(0.1, 5.0; length=64))
    rise = @. 0.5 * (1.0 + tanh((ω - 1.3) / 0.12))
    fall = @. 0.5 * (1.0 - tanh((ω - 3.0) / 0.25))
    δ = π .* rise .* fall
    inverse = cis.(-δ)
    profile = strict_phase_profile(ω, inverse; spec=StrictChargedPhaseSpec(tail_points=8))
    @test profile.anchored_phase[end] ≈ 0.0 atol=1e-14
    @test profile.tail_stable
    @test maximum(abs.(profile.anchored_phase .- δ)) < 1e-6

    unstable = strict_phase_profile(
        ω,
        cis.(-δ .- 0.25 .* (ω ./ ω[end]));
        spec=StrictChargedPhaseSpec(tail_points=8, tail_tolerance=1e-3),
    )
    @test !unstable.tail_stable
end

@testset "Strict BU density uses single-charge domega/pi" begin
    inverse_fn(ω, q) = cis(-π * (0.5 * (1.0 + tanh((ω - (1.2 + 0.02q)) / 0.12)) * 0.5 * (1.0 - tanh((ω - 3.0) / 0.25))))
    strict = strict_charged_bu_density(
        inverse_fn,
        0.6,
        0.2;
        μ=0.0,
        qmax=3.0,
        q_nodes=8,
        omega_min=0.1,
        omega_max=5.0,
        omega_nodes=48,
        require_levinson=false,
    )
    legacy = strict_charged_bu_density(
        inverse_fn,
        0.6,
        0.2;
        μ=0.0,
        qmax=3.0,
        q_nodes=8,
        omega_min=0.1,
        omega_max=5.0,
        omega_nodes=48,
        omega_measure=:legacy_domega_over_2pi,
        require_levinson=false,
    )
    @test strict.accepted
    @test strict.status === :ok
    @test strict.density > 0.0
    @test strict.omega_measure_factor ≈ 1 / π
    @test strict.density ≈ 2.0 * legacy.density rtol=1e-12

    changed = strict_charged_bu_density(
        inverse_fn,
        0.6,
        0.2;
        μ=0.0,
        qmax=3.0,
        q_nodes=12,
        omega_min=0.1,
        omega_max=6.0,
        omega_nodes=72,
        require_levinson=false,
    )
    convergence = strict_density_convergence_gate(strict, changed; rtol=0.05)
    @test convergence.passed
end

@testset "Strict BU gate failure is explicit" begin
    inverse_fn(ω, q) = cis(-π * (0.5 * (1.0 + tanh((ω - 1.0) / 0.2)) * 0.5 * (1.0 - tanh((ω - 2.5) / 0.25))))
    result = strict_charged_bu_density(
        inverse_fn,
        0.4,
        0.2;
        μ=0.0,
        qmax=2.0,
        q_nodes=4,
        omega_min=0.1,
        omega_max=3.0,
        omega_nodes=24,
        threshold=0.8,
        bound_state_count=1,
    )
    @test result.status === :gate_failed
    @test !result.accepted
    @test result.failed_q_count > 0
end

@testset "Strict backend exposes Levinson/Mott and density-sign gates" begin
    ω = collect(range(0.1, 4.0; length=80))
    before_phase = [x <= 2.0 ? π : (x <= 3.0 ? π * (3.0 - x) : 0.0) for x in ω]
    after_phase = zeros(length(ω))
    before = strict_phase_profile(ω, cis.(-before_phase))
    after = strict_phase_profile(ω, cis.(-after_phase))
    before_gate = strict_phase_gate(before; threshold=2.0, bound_state_count=1)
    after_gate = strict_phase_gate(after; threshold=2.0, bound_state_count=0)
    @test before_gate.levinson.phase_passed
    @test after_gate.levinson.phase_passed
    mott = strict_mott_gate(
        before,
        after;
        before_threshold=2.0,
        after_threshold=2.0,
        before_bound_state_count=1,
        after_bound_state_count=0,
    )
    @test mott.passed

    # A falling phase has a negative derivative.  The value is retained for
    # diagnosis, but the strict backend must not accept it as a density.
    falling(ω, q) = cis(-0.8π * (1.0 - tanh((ω - 1.0) / 0.15)) / 2.0)
    negative = strict_charged_bu_density(
        falling,
        0.4,
        0.2;
        qmax=2.0,
        q_nodes=4,
        omega_min=0.1,
        omega_max=3.0,
        omega_nodes=32,
        require_levinson=false,
    )
    @test negative.density < 0.0
    @test negative.status === :invalid_density
    @test !negative.accepted
end

@testset "Charged RPA adapter composes the strict phase backend" begin
    spec = Main.RelaxTime.ChargedRPAKernel.charged_rpa_spec(:K_plus)
    K_a = 0.25
    inverse_fn(ω, q) = cis(-0.7π * exp(-((ω - (1.0 + 0.05q)) / 0.35)^2))
    polarization_fn(ω, q) = (1.0 - inverse_fn(ω, q)) / (4.0 * K_a)
    result = strict_charged_rpa_bu_density(
        spec,
        K_a,
        polarization_fn,
        0.6,
        0.2;
        qmax=2.0,
        q_nodes=4,
        omega_min=0.1,
        omega_max=3.0,
        omega_nodes=32,
        require_levinson=false,
    )
    @test result.accepted
    @test result.status === :ok
    @test result.density > 0.0
    @test result.threshold_policy === :none
end
