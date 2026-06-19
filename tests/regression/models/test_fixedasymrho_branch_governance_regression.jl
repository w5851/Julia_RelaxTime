using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const RUN_FIXEDASYMRHO_BRANCH_GOVERNANCE_REGRESSION =
    get(ENV, "RUN_FIXEDASYMRHO_BRANCH_GOVERNANCE_REGRESSION", "0") in ("1", "true", "TRUE", "yes", "YES")

function _fixedasym_case(label::Symbol, T_MeV::Real, rho_target::Real, low_seed_MeV::AbstractVector{<:Real})
    return (
        label=label,
        T_fm=Float64(T_MeV) / Main.Constants_PNJL.ħc_MeV_fm,
        mode=Main.Models.FixedAsymmetricRho(Float64(rho_target), 0.876, 0.0),
        low_seed=Float64[
            Float64(low_seed_MeV[1]),
            Float64(low_seed_MeV[2]),
            Float64(low_seed_MeV[3]),
            Float64(low_seed_MeV[4]),
            Float64(low_seed_MeV[5]),
            Float64(low_seed_MeV[6]) / Main.Constants_PNJL.ħc_MeV_fm,
            Float64(low_seed_MeV[7]) / Main.Constants_PNJL.ħc_MeV_fm,
            Float64(low_seed_MeV[8]) / Main.Constants_PNJL.ħc_MeV_fm,
        ],
    )
end

@testset "FixedAsymmetricRho branch governance regression" begin
    if !RUN_FIXEDASYMRHO_BRANCH_GOVERNANCE_REGRESSION
        @test_skip "Set RUN_FIXEDASYMRHO_BRANCH_GOVERNANCE_REGRESSION=1 to run slow multiroot branch-governance regression"
    else
        M = Main.Models
        C = Main.Constants_PNJL
        model = M.create_model(:PNJL)
        cases = [
            _fixedasym_case(:t120_rho035, 120.0, 0.35, [
                -5.25687661161,
                -5.24870896285,
                1.90800605097,
                0.0968992840914,
                0.158964020531,
                520.069332664,
                529.922361410,
                26.7743615559,
            ]),
            _fixedasym_case(:t130_rho080, 130.0, 0.80, [
                -5.20811219203,
                -5.19003678068,
                2.04098684988,
                0.186031265548,
                0.236845343269,
                515.596772744,
                526.298601248,
                13.5599246618,
            ]),
        ]

        for case in cases
            @testset "$(case.label)" begin
                low = M.solve(
                    model,
                    case.mode,
                    case.T_fm;
                    seed_guess=case.low_seed,
                    xi=0.0,
                    p_num=8,
                    t_num=4,
                    iterations=80,
                    residual_norm_max=1e-3,
                )
                @test low.converged
                @test low.residual_norm <= 1e-3
                @test low.pressure < 19.0
                @test low.mu_vec[1] * C.ħc_MeV_fm > 500.0

                multiseeds = M.get_all_seeds(M.MultiSeed(), [case.T_fm], case.mode)
                seed_pool = M.build_seed_pool(case.mode;
                    primary_seed=case.low_seed,
                    extra_seed_pool=multiseeds,
                    seed_extend=(seed, _) -> Float64.(seed),
                )
                selected = M.solve_multi(
                    model,
                    case.mode,
                    case.T_fm;
                    seeds=[entry.seed for entry in seed_pool],
                    xi=0.0,
                    p_num=8,
                    t_num=4,
                    iterations=80,
                    residual_norm_max=1e-3,
                    semantic_mode=:ground_state,
                    evaluate_all_attempts=true,
                )

                @test selected.converged
                @test selected.residual_norm <= 1e-3
                @test selected.pressure > low.pressure + 3.0
                @test selected.pressure > 21.0
                @test selected.mu_vec[1] * C.ħc_MeV_fm < 350.0
            end
        end
    end
end
