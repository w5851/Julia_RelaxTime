using Test

const _CHARGED_RPA_CONVERGENCE_SCRIPT = joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
    "audit_charged_rpa_bu_convergence.jl",
)

@testset "charged RPA/BU convergence script contract" begin
    @test isfile(_CHARGED_RPA_CONVERGENCE_SCRIPT)
    source = read(_CHARGED_RPA_CONVERGENCE_SCRIPT, String)
    @test Meta.parseall(source) isa Expr
    @test occursin("build_freezeout_scan_points", source)
    @test occursin("baseline_freezeout", source)
    @test occursin("traversal=:sqrts_descending", source)
    @test occursin("bose_support_gate", source)
    @test occursin("convergence_gate", source)
    @test occursin("four_density_algorithm_labels", source)
    for algorithm in ("stable_particle_limit", "reduced_strict_bw", "q_pole_strict_bw", "phase_shift_bu")
        @test occursin(algorithm, source)
    end
    @test occursin("for algorithm in algorithms, charge in (:plus, :minus)", source)
    @test occursin("ratio=result.ratio", source)
    @test occursin("m_pi_MeV", source)
    @test occursin("mu_pi_MeV", source)
    @test occursin("mu_K_MeV", source)
    @test occursin("rho_B_fm3", source)
    @test occursin("production_candidate_status=\"not_authorized\"", source)
    @test !occursin("kwargs..., continuation_state=nothing", source)
end
