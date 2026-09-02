using Test

const _CHARGED_PHASE_BACKEND_SCRIPT = joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
    "audit_charged_phase_backend.jl",
)

@testset "strict charged phase backend physical diagnostic contract" begin
    @test isfile(_CHARGED_PHASE_BACKEND_SCRIPT)
    source = read(_CHARGED_PHASE_BACKEND_SCRIPT, String)
    @test Meta.parseall(source) isa Expr
    @test occursin("strict_charged_rpa_bu_density", source)
    @test occursin("strict_density_convergence_gate", source)
    @test occursin("joint_convergence_gate", source)
    @test occursin("charged_polarization", source)
    @test occursin("solve_meson_point_from_equilibrium", source)
    @test occursin("prescription=prescription", source)
    @test occursin("ordered_pv_cut", source) || occursin("B0_pv_cut", source)
    @test occursin("CHARGED_PHASE_PRESCRIPTION", source)
    @test occursin("CHARGED_PHASE_OMEGA_MEASURE", source)
    @test occursin("threshold_fn", source)
    @test occursin("count_bound_states", source)
    @test occursin("bound_state_nodes", source)
    @test occursin("omega_measure=_omega_measure()", source)
    @test occursin("production_candidate_status=\"not_authorized\"", source)
    @test occursin("CHARGED_PHASE_\" * name * suffix", source)
    @test occursin("Q_NODES", source)
    @test occursin("OMEGA_NODES", source)
    @test occursin("convergence_passed", source)
    @test occursin("joint_convergence_passed", source)
    @test occursin("tail_failed_q_count", source)
    @test occursin("meson_mass_inv_fm", source)
    @test occursin("levinson_residual_max", source)
    @test occursin("threshold_phase_min", source)
end
