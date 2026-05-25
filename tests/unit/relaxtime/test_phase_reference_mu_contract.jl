using Test

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
