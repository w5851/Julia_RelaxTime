using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _synthetic_phase_classification(; shift=0.0, status=:valid, reason="ok", area=5e-5)
    sres = Models.SShapeResult(
        true,
        310.0 + shift,
        290.0 + shift,
        0.8 + 0.1 * shift,
        2.0 + 0.1 * shift,
        2,
    )
    return (
        status=status,
        reason=reason,
        mu_transition=status == :valid ? 300.0 + shift : nothing,
        rho_hadron=status == :valid ? 1.0 + 0.1 * shift : nothing,
        rho_quark=status == :valid ? 1.8 + 0.1 * shift : nothing,
        area_residual=status == :valid ? area : nothing,
        sres=sres,
    )
end

function _synthetic_phase_result(xi; nonlinear=0.0)
    offset = 2.0 * xi + nonlinear
    boundary = [(
        T_MeV=100.0,
        mu_transition_MeV=300.0 + offset,
        rho_hadron=1.0 + 0.1 * offset,
        rho_quark=2.0 + 0.1 * offset,
        area_residual=5e-5,
        converged=true,
    )]
    spinodal = [(
        T_MeV=100.0,
        mu_spinodal_hadron_MeV=310.0 + offset,
        mu_spinodal_quark_MeV=290.0 + offset,
        rho_spinodal_hadron=0.8 + 0.1 * offset,
        rho_spinodal_quark=2.2 + 0.1 * offset,
    )]
    crossover = [(
        mu_MeV=0.0,
        T_crossover_MeV=180.0 + offset,
        rho=0.3 + 0.01 * offset,
        method="peak",
        converged=true,
        derivative=4.0 + 0.1 * offset,
        variable="phi_u",
    )]
    return Models.PhasePipelineResult(
        xi=xi,
        cep=Models.CEPResult(found=true, T_cep_MeV=130.0 + offset, mu_cep_MeV=295.0 + offset),
        first_order_boundary=boundary,
        spinodal=spinodal,
        crossover_line=crossover,
    )
end

@testset "Phase grid convergence" begin
    tol = Models.PhaseGeometryTolerances(
        position_MeV=0.05,
        density=0.01,
        maxwell_area=1e-4,
        response_rtol=0.05,
    )

    @testset "rho coarse/fine geometry uses separate gates" begin
        coarse = _synthetic_phase_classification(shift=0.0)
        fine = _synthetic_phase_classification(shift=0.02)
        err = Models._compare_phase_geometry(coarse, fine, tol)
        @test err.comparable
        @test err.converged
        @test isapprox(err.position_MeV, 0.02; atol=1e-12)
        @test isapprox(err.density, 0.002; atol=1e-12)
        @test err.maxwell_area == 5e-5

        area_bad = Models._compare_phase_geometry(
            coarse,
            _synthetic_phase_classification(shift=0.02, area=2e-4),
            tol,
        )
        @test !area_bad.converged
        @test area_bad.maxwell_area == 2e-4
    end

    @testset "stable no-S classification converges without invented geometry" begin
        invalid = _synthetic_phase_classification(status=:invalid, reason="no_s_shape")
        err = Models._compare_phase_geometry(invalid, invalid, tol)
        @test err.comparable
        @test err.converged
        @test err.reason == "stable_no_s_shape"
    end

    @testset "temperature midpoint interpolation detects curvature" begin
        left = _synthetic_phase_classification(shift=0.0)
        right = _synthetic_phase_classification(shift=0.04)
        midpoint = _synthetic_phase_classification(shift=0.02)
        linear = Models._phase_geometry_midpoint_error(left, midpoint, right, tol)
        @test linear.converged

        curved = Models._phase_geometry_midpoint_error(
            left,
            _synthetic_phase_classification(shift=0.20),
            right,
            tol,
        )
        @test curved.comparable
        @test !curved.converged
        @test curved.position_MeV > tol.position_MeV
    end

    @testset "xi midpoint audit includes boundary, spinodal, crossover and CEP" begin
        left = _synthetic_phase_result(-0.1)
        midpoint = _synthetic_phase_result(0.0)
        right = _synthetic_phase_result(0.1)
        linear = Models._phase_result_midpoint_error(left, midpoint, right, tol)
        @test linear.comparable
        @test linear.converged

        curved = Models._phase_result_midpoint_error(
            left,
            _synthetic_phase_result(0.0; nonlinear=0.2),
            right,
            tol,
        )
        @test curved.comparable
        @test !curved.converged
        @test curved.position_MeV > tol.position_MeV
    end

    @test_throws ArgumentError Models._validate_phase_geometry_tolerances(
        Models.PhaseGeometryTolerances(position_MeV=0.0),
    )
end
