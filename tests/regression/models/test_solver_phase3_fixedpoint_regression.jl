using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const HBARC_MEV_FM = 197.327

@inline _to_fm_inv(x_mev::Real) = Float64(x_mev) / HBARC_MEV_FM

function _solve_point(model, mode, T_MeV::Real; mu_MeV::Real=0.0, seed)
    T_fm = _to_fm_inv(T_MeV)
    kwargs = (
        seed_guess=seed,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        iterations=120,
    )
    if mode isa Models.FixedMu
        return Models.solve_constraint(model, mode, T_fm; μ_fm=_to_fm_inv(mu_MeV), kwargs...)
    end
    return Models.solve_constraint(model, mode, T_fm; kwargs...)
end

@testset "solver phase3 fixedpoint regression" begin
    model = Models.create_model(:PNJL)
    P = Models.pnjl_module()

    cases = [
        (
            name=:fixedmu_v1,
            mode=Models.FixedMu(),
            T_MeV=120.0,
            mu_MeV=200.0,
            seed=copy(P.HADRON_SEED_5),
            expect=(
                converged=true,
                pressure=21.61145770689436,
                rho_norm=0.0047488846370136945,
                entropy=0.06254708485157146,
                energy=-21.558981456992434,
                residual_norm=1.1783690026873215e-11,
            ),
        ),
        (
            name=:fixedrho_v1,
            mode=Models.FixedRho(0.2),
            T_MeV=120.0,
            mu_MeV=0.0,
            seed=copy(P.HADRON_SEED_8),
            expect=(
                converged=true,
                pressure=21.626374914491763,
                rho_norm=0.20000000099770002,
                entropy=0.2902034886921015,
                energy=-21.318354182416222,
                residual_norm=1.206180825306215e-9,
            ),
        ),
        (
            name=:fixedasymrho_v1,
            mode=Models.FixedAsymmetricRho(0.05, 1.0, 0.0),
            T_MeV=110.0,
            mu_MeV=0.0,
            seed=copy(P.HADRON_SEED_8),
            expect=(
                converged=false,
                pressure=-Inf,
                rho_norm=NaN,
                entropy=NaN,
                energy=NaN,
                residual_norm=Inf,
                selection_reason=:no_candidate_passed_constraints,
                failed_constraints=Symbol[:solver_failed],
            ),
        ),
    ]

    for case in cases
        got = _solve_point(model, case.mode, case.T_MeV; mu_MeV=case.mu_MeV, seed=case.seed)
        @test got.converged == case.expect.converged

        if case.expect.converged
            @test isapprox(got.pressure, case.expect.pressure; rtol=1e-8, atol=1e-10)
            @test isapprox(got.rho_norm, case.expect.rho_norm; rtol=1e-8, atol=1e-10)
            @test isapprox(got.entropy, case.expect.entropy; rtol=1e-8, atol=1e-10)
            @test isapprox(got.energy, case.expect.energy; rtol=1e-8, atol=1e-10)
            @test isapprox(got.residual_norm, case.expect.residual_norm; rtol=1e-6, atol=1e-12)
        else
            @test got.selection_reason == case.expect.selection_reason
            @test Symbol.(get(got, :failed_constraints, Symbol[])) == case.expect.failed_constraints
            @test !isfinite(got.pressure)
            @test !isfinite(got.rho_norm)
            @test !isfinite(got.entropy)
            @test !isfinite(got.energy)
            @test !isfinite(got.residual_norm)
        end
    end
end
