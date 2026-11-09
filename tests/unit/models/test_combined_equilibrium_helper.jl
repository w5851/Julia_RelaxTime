using Test

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

module CombinedEquilibriumHelperFixture
const Models = Main.Models
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "combined_meson_density_equilibrium.jl"))
struct RecordingModel <: Models.AbstractQCDModel end
function Models.solve_multi(::RecordingModel, mode, T; kwargs...)
    return (; mode, T, options=(; kwargs...))
end
end

@testset "combined equilibrium helper preserves seed and pressure policy" begin
    H = CombinedEquilibriumHelperFixture
    T = 120.0 / 197.327
    mode = Models.FixedAsymmetricRho(0.35, 0.876, 0.0)
    catalog = Models.get_all_seeds(Models.MultiSeed(), [T], mode)
    seeds = H._trho_asymmetric_seed_candidates(T, mode, nothing)
    @test seeds == Models._unique_path_seeds(catalog)
    previous = copy(first(catalog))
    previous[1] += 0.1
    warm = H._trho_asymmetric_seed_candidates(T, mode, previous)
    @test first(warm) == previous
    @test warm[2:end] == seeds
    @test H._trho_asymmetric_seed_candidates(T, mode, first(catalog)) == seeds
    result = H._solve_trho_asymmetric_equilibrium(H.RecordingModel(), mode, T;
        seeds=warm, xi=0.0, p_num=8, t_num=4, iterations=300, residual_norm_max=1e-6)
    @test result.mode === mode
    @test result.T == T
    @test result.options.seeds === warm
    @test result.options.semantic_mode === :ground_state
    @test result.options.evaluate_all_attempts
    @test result.options.iterations == 300
    @test result.options.residual_norm_max == 1e-6
end
