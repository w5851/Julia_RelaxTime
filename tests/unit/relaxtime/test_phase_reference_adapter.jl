using Test

const _PHASE_REFERENCE_ADAPTER_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "relaxtime", "phase_reference_adapter.jl",
))
if !isdefined(Main, :PhaseReferenceAdapter)
    Base.include(Main, _PHASE_REFERENCE_ADAPTER_PATH)
end

const PRA = Main.PhaseReferenceAdapter

function _write_candidate_fixture(root::AbstractString; unresolved::Bool=false, duplicate::Bool=false, nonfinite::Bool=false)
    mkpath(joinpath(root, "strict", "tables"))
    write(joinpath(root, "manifest.json"), """
    {"schema_version":"pnjl_issue130_phase_reference_import_v1","reference_status":"candidate","runtime_consumption":false}
    """)
    write(joinpath(root, "strict", "manifest.json"), "{\"layer\":\"strict_reference_v1\"}")
    grid_status = unresolved ? "rho:rho_geometry_not_converged" : "ok"
    area = nonfinite ? "NaN" : "1e-6"
    boundary = """surface,xi,T_MeV,mu_MeV,rho_hadron,rho_quark,area_residual,grid_status,grid_unresolved,layer,status,geometry_converged,finite_and_converged
    maxwell,0.0,100.0,300.0,1.0,2.0,$area,$grid_status,$(unresolved ? "True" : "False"),strict_reference_v1,$(unresolved ? "native_unresolved" : "native"),$(unresolved ? "False" : "True"),True
    """
    duplicate && (boundary *= "maxwell,0.0,100.0,300.0,1.0,2.0,$area,$grid_status,$(unresolved ? "True" : "False"),strict_reference_v1,native,True,True\n")
    write(joinpath(root, "strict", "tables", "maxwell_surface_strict_reference_v1.csv"), boundary)
    write(joinpath(root, "strict", "tables", "crossover_surface_strict_reference_v1.csv"), """surface,xi,mu_MeV,T_MeV,rho,mu_CEP_proxy_MeV,layer,status,physical_region
    crossover,0.0,100.0,160.0,1.0,300.0,strict_reference_v1,native,crossover_below_cep
    """)
    write(joinpath(root, "strict", "tables", "cep_boundary_strict_reference_v1.csv"), """surface,xi,mu_CEP_proxy_MeV,T_low_MeV,T_high_MeV,T_midpoint_MeV,layer,status,boundary_mode
    cep_boundary,0.0,300.0,120.0,120.0625,120.03125,strict_reference_v1,bracket_preserved,estimated_midpoint
    """)
    write(joinpath(root, "strict", "tables", "spinodal_surface_strict_reference_v1.csv"), """surface,xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,layer,status
    spinodal,0.0,100.0,320.0,280.0,strict_reference_v1,native
    """)
    return root
end

@testset "phase-reference candidate adapter contract" begin
    root = _write_candidate_fixture(mktempdir())
    source = PRA.load_phase_reference(root)
    @test PRA.source_kind(source) === :candidate
    @test PRA.source_layer(source) === :strict
    @test !PRA.source_runtime_enabled(source)
    @test source.diagnostics.row_counts["boundary"] == 1
    @test source.diagnostics.uncertified_rows == 0
    @test PRA.available_xi(source, :boundary) == [0.0]

    @test_throws PRA.PhaseReferenceAdapterError PRA.boundary_data(source, 0.0)
    # A fully certified fixture can be explicitly gated for a solver-free
    # consumer parity query.
    gated = PRA.load_phase_reference(root; allow_runtime=true)
    boundary = PRA.boundary_data(gated, 0.0)
    @test boundary.T_values == [100.0]
    @test boundary.mu_values == [300.0]
    @test boundary.muB_CEP == 900.0
    @test PRA.crossover_rows(gated, 0.0)[1].muB_MeV == 300.0
end

@testset "phase-reference candidate rejects unresolved/nonfinite/duplicate rows" begin
    @test PRA.load_phase_reference(_write_candidate_fixture(mktempdir(); unresolved=true)).diagnostics.uncertified_rows == 1
    @test_throws PRA.PhaseReferenceAdapterError PRA.load_phase_reference_runtime(
        _write_candidate_fixture(mktempdir(); unresolved=true),
    )
    @test_throws PRA.PhaseReferenceAdapterError _write_candidate_fixture(mktempdir(); nonfinite=true) |> PRA.load_phase_reference
    @test_throws PRA.PhaseReferenceAdapterError _write_candidate_fixture(mktempdir(); duplicate=true) |> PRA.load_phase_reference
end

@testset "legacy and candidate views preserve muq to muB parity" begin
    candidate_root = _write_candidate_fixture(mktempdir())
    candidate = PRA.load_phase_reference_runtime(candidate_root)
    candidate_boundary = PRA.boundary_data(candidate, 0.0)

    legacy_root = mktempdir()
    boundary = joinpath(legacy_root, "boundary.csv")
    cep = joinpath(legacy_root, "cep.csv")
    crossover = joinpath(legacy_root, "crossover.csv")
    spinodals = joinpath(legacy_root, "spinodals.csv")
    write(boundary, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark\n0.0,100.0,300.0,1.0,2.0\n")
    write(cep, "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV\n0.0,120.03125,300.0,900.0\n")
    write(crossover, "xi,mu_MeV,T_crossover_MeV,rho\n0.0,100.0,160.0,1.0\n")
    write(spinodals, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV\n0.0,100.0,320.0,280.0\n")
    legacy = PRA.load_legacy_phase_reference(
        boundary_path=boundary, cep_path=cep, crossover_path=crossover, spinodals_path=spinodals,
    )
    legacy_boundary = PRA.boundary_data(legacy, 0.0)
    @test legacy_boundary.T_values == candidate_boundary.T_values
    @test legacy_boundary.mu_values == candidate_boundary.mu_values
    @test legacy_boundary.muB_CEP == candidate_boundary.muB_CEP
end
