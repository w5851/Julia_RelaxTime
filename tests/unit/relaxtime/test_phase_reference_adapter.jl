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

function _write_v2_candidate_fixture(root::AbstractString;
    accepted_noncertified::Bool=false,
    accepted_authorized::Bool=false,
    unresolved::Bool=false,
)
    _write_candidate_fixture(root; unresolved=unresolved)
    manifest_path = joinpath(root, "manifest.json")
    manifest = replace(read(manifest_path, String),
        "pnjl_issue130_phase_reference_import_v1" => "pnjl_issue130_phase_reference_v2")
    manifest = replace(manifest,
        "\"reference_status\":\"candidate\"" =>
            "\"reference_status\":\"candidate\",\"promotion_status\":\"accepted_for_downstream\",\"downstream_default_layer\":\"accepted\"")
    write(manifest_path, manifest)

    render_names = Dict(
        :boundary => "maxwell_surface_render.csv",
        :crossover => "crossover_surface_render.csv",
        :cep => "cep_boundary_render.csv",
        :spinodals => "spinodal_surface_render.csv",
    )
    strict_names = Dict(
        :boundary => "maxwell_surface_strict_reference_v1.csv",
        :crossover => "crossover_surface_strict_reference_v1.csv",
        :cep => "cep_boundary_strict_reference_v1.csv",
        :spinodals => "spinodal_surface_strict_reference_v1.csv",
    )
    accepted_names = Dict(
        :boundary => "maxwell_surface_accepted_phase_map_v1.csv",
        :crossover => "crossover_surface_accepted_phase_map_v1.csv",
        :cep => "cep_boundary_accepted_phase_map_v1.csv",
        :spinodals => "spinodal_surface_accepted_phase_map_v1.csv",
    )
    additions = ",source_status,acceptance_status,extrapolation,coverage_status,acceptance_scope"
    for layer in ("render", "accepted")
        mkpath(joinpath(root, layer, "tables"))
        write(joinpath(root, layer, "manifest.json"), "{\"layer\":\"$(layer)\",\"runtime_consumption\":false}")
    end
    for table in keys(strict_names)
        source = joinpath(root, "strict", "tables", strict_names[table])
        render = joinpath(root, "render", "tables", render_names[table])
        cp(source, render; force=true)
        lines = split(chomp(read(source, String)), '\n')
        header, row = first(lines), last(lines)
        acceptance_status = accepted_authorized ?
            "author_accepted_for_downstream" : "candidate_pending_author_review"
        accepted_row_status = accepted_noncertified ?
            "interpolated_noncertified,$acceptance_status,False,interpolated_common_support,downstream_phase_map_candidate" :
            "strict_certified,$acceptance_status,False,native_support,downstream_phase_map_candidate"
        accepted = join((header * additions, row * "," * accepted_row_status), '\n') * "\n"
        write(joinpath(root, "accepted", "tables", accepted_names[table]), accepted)
    end
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
    @test length(source.diagnostics.candidate_manifest_sha256) == 64

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
    strict_runtime = PRA.load_phase_reference_runtime(
        _write_candidate_fixture(mktempdir(); unresolved=true),
    )
    @test PRA.source_layer(strict_runtime) === :strict
    @test strict_runtime.diagnostics.runtime_view == "strict_certified_only"
    @test isempty(strict_runtime.tables[:boundary])
    @test_throws PRA.PhaseReferenceAdapterError _write_candidate_fixture(mktempdir(); nonfinite=true) |> PRA.load_phase_reference
    @test_throws PRA.PhaseReferenceAdapterError _write_candidate_fixture(mktempdir(); duplicate=true) |> PRA.load_phase_reference
end

@testset "v2 exposes strict/render/accepted and keeps downstream layers out of runtime" begin
    root = _write_v2_candidate_fixture(mktempdir(); accepted_noncertified=true)
    render = PRA.load_phase_reference(root; layer=:render)
    @test PRA.source_layer(render) === :render
    @test render.diagnostics.candidate_schema_version == "pnjl_issue130_phase_reference_v2"
    @test length(render.tables[:spinodals]) == 1

    accepted = PRA.load_phase_reference(root; layer=:accepted)
    @test PRA.source_layer(accepted) === :accepted
    @test accepted.diagnostics.uncertified_rows == 4
    @test_throws PRA.PhaseReferenceAdapterError PRA.load_phase_reference_runtime(root; layer=:render)
    @test_throws PRA.PhaseReferenceAdapterError PRA.load_phase_reference_runtime(root; layer=:accepted)
    @test_throws PRA.PhaseReferenceAdapterError PRA.load_phase_reference(root; layer=:derived)
    @test !isdefined(PRA, :load_phase_reference_runtime_with_fallback)
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

@testset "legacy fallback and rollback entrypoint are retired" begin
    candidate_root = _write_candidate_fixture(mktempdir(); unresolved=true)
    legacy_root = mktempdir()
    boundary = joinpath(legacy_root, "boundary.csv")
    cep = joinpath(legacy_root, "cep.csv")
    crossover = joinpath(legacy_root, "crossover.csv")
    spinodals = joinpath(legacy_root, "spinodals.csv")
    write(boundary, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark\n0.0,90.0,301.0,1.0,2.0\n")
    write(cep, "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV\n0.0,121.0,299.0,897.0\n")
    write(crossover, "xi,mu_MeV,T_crossover_MeV,rho\n0.0,100.0,160.0,1.0\n")
    write(spinodals, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV\n0.0,100.0,320.0,280.0\n")

    @test !isdefined(PRA, :load_phase_reference_runtime_with_fallback)

    rollback = PRA.load_legacy_phase_reference(
        boundary_path=boundary, cep_path=cep, crossover_path=crossover, spinodals_path=spinodals,
    )
    @test PRA.source_kind(rollback) === :legacy
    @test PRA.source_summary(rollback).runtime_view == "legacy"
    @test PRA.boundary_data(rollback, 0.0).mu_values == [301.0]
end

@testset "accepted non-certified rows are the primary runtime source" begin
    root = _write_v2_candidate_fixture(
        mktempdir(); accepted_noncertified=true, accepted_authorized=true,
    )
    legacy_root = mktempdir()
    boundary = joinpath(legacy_root, "boundary.csv")
    cep = joinpath(legacy_root, "cep.csv")
    crossover = joinpath(legacy_root, "crossover.csv")
    spinodals = joinpath(legacy_root, "spinodals.csv")
    write(boundary, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark\n0.0,90.0,301.0,1.0,2.0\n")
    write(cep, "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV\n0.0,121.0,299.0,897.0\n")
    write(crossover, "xi,mu_MeV,T_crossover_MeV,rho\n0.0,250.0,160.0,1.0\n0.0,350.0,150.0,1.1\n")
    write(spinodals, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV\n0.0,100.0,320.0,280.0\n")

    runtime = PRA.load_phase_reference_runtime(root; layer=:accepted)
    summary = PRA.source_summary(runtime)
    @test summary.runtime_view == "accepted_primary"
    @test summary.primary_layer == "accepted"
    @test !summary.fallback_enabled
    @test summary.fallback_order == "accepted_primary"
    @test summary.legacy_fallback_row_counts["boundary"] == 0
    @test summary.legacy_fallback_row_counts["crossover"] == 0
    accepted_row = only(runtime.tables[:boundary])
    @test !accepted_row.certified
    @test accepted_row.runtime_eligible
    @test accepted_row.runtime_source_layer == "accepted_primary"
    @test accepted_row.source_status == "interpolated_noncertified"
    @test length(summary.candidate_manifest_sha256) == 64
    @test length(summary.candidate_layer_manifest_sha256) == 64
    @test PRA.boundary_data(runtime, 0.0).mu_values == [300.0]
end

@testset "repository default uses accepted primary and strict is explicit" begin
    project_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    runtime = PRA.load_default_phase_reference_runtime(project_root=project_root)
    summary = PRA.source_summary(runtime)
    @test PRA.source_kind(runtime) === :candidate
    @test PRA.source_layer(runtime) === :accepted
    @test summary.runtime_view == "accepted_primary"
    @test !summary.fallback_enabled
    @test summary.fallback_order == "accepted_primary"
    @test sum(values(summary.accepted_primary_row_counts)) > 0
    @test sum(values(summary.legacy_fallback_row_counts)) == 0

    strict = PRA.load_default_phase_reference_runtime(
        project_root=project_root, source=:strict, layer=:strict,
    )
    @test PRA.source_kind(strict) === :candidate
    @test PRA.source_layer(strict) === :strict
    @test PRA.source_summary(strict).runtime_view == "strict_certified_only"
    @test !PRA.source_summary(strict).fallback_enabled
    @test_throws PRA.PhaseReferenceAdapterError PRA.load_default_phase_reference_runtime(
        project_root=project_root, source=:strict, layer=:accepted,
    )
end
