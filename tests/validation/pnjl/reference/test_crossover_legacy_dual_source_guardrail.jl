using Test

if !isdefined(Main, :_load_crossover_reference_targets)
    Base.include(Main, joinpath(@__DIR__, "..", "pnjl_validation_helpers.jl"))
end

const PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_HBARC_MEV_FM = 197.327
const PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_TARGETS_PATH = validation_targets_path(
    "pnjl",
    "reference",
    "pnjl_crossover_legacy_dual_source_targets_v1.csv",
)

using .Models

@testset "PNJL crossover legacy dual-source guardrail" begin
    targets = _load_crossover_reference_targets(PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_TARGETS_PATH)
    @test !isempty(targets)

    T_range = (
        120.0 / PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_HBARC_MEV_FM,
        220.0 / PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_HBARC_MEV_FM,
    )

    for row in targets
        result = Models.detect_crossover(
            row.mu_MeV / PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_HBARC_MEV_FM,
            T_range;
            method=Symbol(row.method),
            variable=Symbol(row.variable),
            xi=row.xi,
            p_num=8,
            t_num=4,
            n_scan=20,
            max_iter=20,
            solver_backend=:models,
        )

        @test result.found
        @test result.T_crossover !== nothing

        actual_T_MeV = Float64(result.T_crossover) * PNJL_CROSSOVER_LEGACY_DUAL_SOURCE_HBARC_MEV_FM
        @test actual_T_MeV >= row.lower_T_MeV
        @test actual_T_MeV <= row.upper_T_MeV
    end
end
