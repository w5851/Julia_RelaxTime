using Test

const MAG_SCAN_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(MAG_SCAN_PROJECT_ROOT, "src", "models", "Models.jl"))
end

const _MS = Models.MagneticScan

@testset "MagneticScan contract" begin
    @test !isempty(_MS.DEFAULT_T_VALUES)
    @test all(>(0.0), _MS.DEFAULT_T_VALUES)
    @test !isempty(_MS.DEFAULT_EB_VALUES)
    @test isdefined(Models, :run_magnetic_scan)
    @test Models.run_magnetic_scan isa Function

    @test _MS._validate_magnetic_scan_inputs(
        [150.0], [0.0], [0.0, 2.0e4], [0.0],
        :PNJLMagnetic, :fixed_mu, 4, 1, 1,
    ) === nothing
    @test_throws ArgumentError _MS._validate_magnetic_scan_inputs(
        [0.0], [0.0], [0.0], [0.0],
        :PNJLMagnetic, :fixed_mu, 4, 1, 1,
    )
    @test_throws ArgumentError _MS._validate_magnetic_scan_inputs(
        [150.0], [0.0], [0.0], [0.1],
        :PNJLMagnetic, :fixed_mu, 4, 1, 1,
    )
    @test_throws ArgumentError _MS._validate_magnetic_scan_inputs(
        [150.0], [0.0], [0.0], [0.0],
        :PNJL, :fixed_mu, 4, 1, 1,
    )
    @test_throws ArgumentError _MS._validate_magnetic_scan_inputs(
        [150.0], [0.0], [0.0], [0.0],
        :PNJLMagnetic, :fixed_rho, 4, 1, 1,
    )

    @test _MS._scan_key(150.0, 60.0, 2.0e4, 0.0) ==
        (150.0, 60.0, 20000.0, 0.0)

    tmp = mktempdir()
    selected = joinpath(tmp, "selected.csv")
    write(selected, "T_MeV,mu_MeV,eB_MeV2,xi,message\n150,60,20000,0,ok\n")
    @test (150.0, 60.0, 20000.0, 0.0) in _MS._load_completed_keys(selected)

    @test occursin("_candidates.csv", _MS._derived_candidates_output(selected))
    @test occursin("selected_candidate_index", _MS.SELECTED_HEADER)
    @test occursin("candidate_index", _MS.CANDIDATE_HEADER)

    @test_throws ArgumentError Models.run_tmu_scan(
        T_values=[150.0], mu_values=[0.0], xi_values=[0.0],
        model_kind=:PNJLMagnetic, output_path=joinpath(tmp, "tmu.csv"),
    )
    @test_throws ArgumentError Models.run_trho_scan(
        T_values=[150.0], rho_values=[0.2], xi_values=[0.0],
        model_kind=:PNJLMagnetic, output_path=joinpath(tmp, "trho.csv"),
    )
end
