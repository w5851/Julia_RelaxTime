using Test

const _PHASE_REFERENCE_ADAPTER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "phase_reference_adapter.jl"))
if !isdefined(Main, :PhaseReferenceAdapter)
    Base.include(Main, _PHASE_REFERENCE_ADAPTER_PATH)
end
const PRA = Main.PhaseReferenceAdapter

const _PHASE_EQUILIBRIUM_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "gap_transport_scan_phase_equilibrium.jl"))

function _write_phase_csv(path::String, header::String, rows::Vector{String})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, row)
        end
    end
end

function _write_candidate_phase_bundle(root::String)
    tables = joinpath(root, "strict", "tables")
    mkpath(tables)
    write(joinpath(root, "manifest.json"), "{\"schema_version\":\"pnjl_issue130_phase_reference_import_v1\",\"reference_status\":\"candidate\",\"runtime_consumption\":false}")
    write(joinpath(root, "strict", "manifest.json"), "{\"layer\":\"strict_reference_v1\"}")
    write(joinpath(tables, "maxwell_surface_strict_reference_v1.csv"), "surface,xi,T_MeV,mu_MeV,rho_hadron,rho_quark,area_residual,grid_unresolved,layer,status,geometry_converged,finite_and_converged\nmaxwell,0.0,120.0,300.0,1.0,2.0,1e-6,false,strict_reference_v1,native,true,true\n")
    write(joinpath(tables, "crossover_surface_strict_reference_v1.csv"), "surface,xi,mu_MeV,T_MeV,rho,mu_CEP_proxy_MeV,layer,status,physical_region\ncrossover,0.0,90.0,140.0,1.0,300.0,strict_reference_v1,native,crossover_below_cep\n")
    write(joinpath(tables, "cep_boundary_strict_reference_v1.csv"), "surface,xi,mu_CEP_proxy_MeV,T_low_MeV,T_high_MeV,T_midpoint_MeV,layer,status,boundary_mode\ncep_boundary,0.0,300.0,120.0,120.0625,120.03125,strict_reference_v1,bracket_preserved,estimated_midpoint\n")
    write(joinpath(tables, "spinodal_surface_strict_reference_v1.csv"), "surface,xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,layer,status\nspinodal,0.0,100.0,280.0,320.0,strict_reference_v1,native\n")
    return root
end

@testset "phase CEP reference mu contract" begin
    @test isfile(_PHASE_EQUILIBRIUM_PATH)
    if !isdefined(Main, :GapTransportScanPhaseEquilibrium)
        Base.include(Main, _PHASE_EQUILIBRIUM_PATH)
    end

    tmp = mktempdir()
    boundary_csv = joinpath(tmp, "boundary.csv")
    _write_phase_csv(
        boundary_csv,
        "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark",
        ["0.0,120.0,300.0,1.0,2.0"],
    )

    explicit_cep_csv = joinpath(tmp, "cep_explicit.csv")
    _write_phase_csv(
        explicit_cep_csv,
        "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV",
        ["0.0,130.0,292.0,876.0"],
    )

    data = Main.GapTransportScanPhaseEquilibrium.load_phase_boundary_data(
        0.0;
        boundary_path=boundary_csv,
        cep_path=explicit_cep_csv,
    )
    @test data.T_CEP == 130.0
    @test data.muq_CEP == 292.0
    @test data.muB_CEP == 876.0
    @test data.mu_CEP == data.muq_CEP

    legacy_cep_csv = joinpath(tmp, "cep_legacy.csv")
    _write_phase_csv(
        legacy_cep_csv,
        "xi,T_CEP_MeV,mu_CEP_MeV",
        ["0.0,130.0,292.0"],
    )

    legacy = Main.GapTransportScanPhaseEquilibrium.load_phase_boundary_data(
        0.0;
        boundary_path=boundary_csv,
        cep_path=legacy_cep_csv,
    )
    @test legacy.muq_CEP == 292.0
    @test legacy.muB_CEP == 876.0
    @test legacy.mu_CEP == legacy.muq_CEP
end

@testset "transport consumer preserves legacy/candidate phase parity" begin
    candidate_root = _write_candidate_phase_bundle(mktempdir())
    candidate = PRA.load_phase_reference_runtime(candidate_root)

    candidate_data = Main.GapTransportScanPhaseEquilibrium.load_phase_boundary_data(
        0.0;
        boundary_path="",
        cep_path="",
        phase_reference=candidate,
        phase_reference_mode=:runtime,
    )
    @test candidate_data.T_values == [120.0]
    @test candidate_data.mu_values == [300.0]
    @test candidate_data.muq_CEP == 300.0
    @test candidate_data.muB_CEP == 900.0

    candidate_cross, candidate_xi = Main.GapTransportScanPhaseEquilibrium.interpolate_crossover_temperature(
        0.0,
        90.0;
        phase_reference=candidate,
        phase_reference_mode=:runtime,
    )
    @test candidate_cross == 140.0
    @test candidate_xi == 0.0
end
