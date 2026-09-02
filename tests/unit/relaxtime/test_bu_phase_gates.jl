"""Pure unit coverage for the strict real-axis BU measure and phase gates."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _RELAXTIME_PATH_BU_GATES = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH_BU_GATES)
end

using Main.RelaxTime.BUPhaseGates: STRICT_SINGLE_CHARGE_OMEGA_MEASURE,
                                        LEGACY_POSITIVE_ENERGY_OMEGA_MEASURE,
                                        bu_omega_measure,
                                        bu_omega_measure_factor,
                                        anchor_phase_high_energy,
                                        count_subthreshold_roots,
                                        count_bound_states,
                                        continue_bound_state_counts,
                                        levinson_phase_gate,
                                        mott_phase_gate,
                                        bose_support_gate,
                                        convergence_gate,
                                        joint_convergence_gate,
                                        four_density_algorithm_labels

@testset "BU omega measure contract" begin
    @test bu_omega_measure(:strict) === STRICT_SINGLE_CHARGE_OMEGA_MEASURE
    @test bu_omega_measure(:legacy) === LEGACY_POSITIVE_ENERGY_OMEGA_MEASURE
    @test bu_omega_measure_factor(:single_charge_domega_over_pi) ≈ inv(pi)
    @test bu_omega_measure_factor(:legacy_domega_over_2pi) ≈ inv(2pi)
    @test bu_omega_measure_factor(:strict) ≈ 2 * bu_omega_measure_factor(:legacy)
    @test_throws ArgumentError bu_omega_measure(:unknown)
end

@testset "Independent bound-state counting and q continuation" begin
    inverse_fn(omega, q) = (omega - (0.4 + 0.1q)) * (omega - (0.9 + 0.1q))
    one_q = count_bound_states(inverse_fn, 0.5, 1.5; omega_nodes=80)
    @test one_q.independent
    @test one_q.status === :ok
    @test one_q.count == 2
    @test length(one_q.brackets) == 2

    continued = continue_bound_state_counts(
        inverse_fn,
        [0.0, 0.5, 1.0],
        q -> 1.5 + 0.1q;
        omega_nodes=80,
    )
    @test length(continued) == 3
    @test all(row.count == 2 for row in continued)
    @test ismissing(continued[1].continuation_delta)
    @test all(row.continuation_delta == 0 for row in continued[2:end])

    complex = count_bound_states(
        (omega, q) -> inverse_fn(omega, q) + 1.0e-3im,
        0.0,
        1.5;
        omega_nodes=80,
        imag_tolerance=1.0e-6,
    )
    @test complex.status === :complex_subthreshold
    @test !complex.passed
end

@testset "High-energy phase anchoring unwraps from the tail" begin
    omega = collect(0.0:1.0:5.0)
    desired = [pi, pi, 0.8pi, 0.3pi, 0.01, 0.0]
    raw = mod.(desired .+ 0.37 .+ pi, 2pi) .- pi
    anchored = anchor_phase_high_energy(omega, raw; tail_points=2)
    @test anchored.high_energy_phase_after_anchor ≈ 0.0 atol=1e-14
    @test anchored.applied_shift ≈ -0.37 atol=1e-12
    @test anchored.anchored_phase ≈ desired atol=1e-12
    @test anchored.tail_span ≈ 0.01 atol=1e-12
    @test_throws ArgumentError anchor_phase_high_energy([0.0, 0.0], [0.0, 0.0])
end

@testset "Subthreshold sign-change root count" begin
    omega = [0.0, 0.5, 1.5, 2.5, 3.0]
    inverse_values = ComplexF64[(w - 1.0) * (w - 2.0) for w in omega]
    roots = count_subthreshold_roots(omega, inverse_values, 2.8)
    @test roots.passed
    @test roots.count == 2
    @test length(roots.brackets) == 2

    complex_values = inverse_values .+ 1.0e-3im
    complex_roots = count_subthreshold_roots(
        omega,
        complex_values,
        2.8;
        imag_tolerance=1.0e-6,
    )
    @test !complex_roots.passed
    @test complex_roots.status === :complex_subthreshold
end

@testset "Levinson and Mott gates" begin
    omega = collect(0.0:1.0:9.0)
    offset_before = 0.41
    phase_before = [pi, pi, pi, pi, 0.7pi, 0.3pi, 0.05, 0.01, 0.0, 0.0] .+ offset_before
    before = levinson_phase_gate(
        omega,
        phase_before,
        3.0;
        bound_state_count=1,
        tail_points=3,
    )
    @test before.passed
    @test before.threshold_phase ≈ pi atol=1e-12
    @test abs(before.levinson_residual) < 1e-12

    phase_after = fill(-0.23, length(omega))
    after = levinson_phase_gate(
        omega,
        phase_after,
        3.0;
        bound_state_count=0,
        tail_points=3,
    )
    @test after.passed
    transition = mott_phase_gate(before, after)
    @test transition.passed
    @test transition.bound_state_count_drop == 1
    @test transition.threshold_phase_drop ≈ pi atol=1e-12

    wrong_after = merge(after, (bound_state_count=1,))
    @test !mott_phase_gate(before, wrong_after).passed
end

@testset "Bose support and convergence gates" begin
    safe = bose_support_gate(0.7, 0.2; omega_min=0.3, omega_max=3.0)
    @test safe.passed
    @test safe.status === :safe_normal_domain
    unsafe = bose_support_gate(0.7, 0.8; omega_min=0.3, omega_max=3.0)
    @test !unsafe.passed
    @test unsafe.status === :unsafe_bose_domain
    @test_throws ArgumentError bose_support_gate(0.7, 0.2; omega_min=3.0, omega_max=0.3)

    close = convergence_gate(1.0, 1.0005; rtol=1e-3)
    @test close.passed
    far = convergence_gate(1.0, 1.01; rtol=1e-3)
    @test !far.passed
    joint = joint_convergence_gate([
        (density=1.0, accepted=true, tail_stable=true, eta_inv_fm=0.01, q_nodes=8),
        (density=1.0005, accepted=true, tail_stable=true, eta_inv_fm=0.005, q_nodes=16),
    ]; rtol=1e-3)
    @test joint.passed
    @test joint.sample_count == 2
    rejected = joint_convergence_gate([
        (density=1.0, accepted=true, tail_stable=true),
        (density=1.0005, accepted=false, tail_stable=true),
    ]; rtol=1e-3)
    @test !rejected.passed
    @test four_density_algorithm_labels() == (
        :stable_particle_limit,
        :reduced_strict_bw,
        :q_pole_strict_bw,
        :phase_shift_bu,
    )
end
