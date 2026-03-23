using Test

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_SCRIPT = joinpath(_PROJECT_ROOT, "src", "models", "Models.jl")
const _CONSTANTS_SCRIPT = joinpath(_PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl")

isfile(_MODELS_SCRIPT) || error("Models script missing: $(_MODELS_SCRIPT)")
isfile(_CONSTANTS_SCRIPT) || error("Constants script missing: $(_CONSTANTS_SCRIPT)")

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_SCRIPT)
end

if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_SCRIPT)
end

const _MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

@testset "Mixed meson root quality with continuation" begin
    muB_MeV = 0.0
    muq_fm = (muB_MeV / 3.0) / _HBARC_MEV_FM

    seed_state = nothing
    for T_MeV in (200.0, 220.0, 240.0)
        T_fm = T_MeV / _HBARC_MEV_FM
        res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=40,),
            meson_seed_state=seed_state,
            meson_root_policy=(residual_norm_max=1e-6, require_converged=true, use_trust_region_fallback=true),
        )

        seed_state = res.meson_seed_state
        @test hasproperty(res.meson_results[:eta], :root_quality)
        @test hasproperty(res.meson_results[:eta_prime], :root_quality)
        @test hasproperty(res.meson_results[:eta], :root_diagnostics)
        @test hasproperty(res.meson_results[:eta_prime], :root_diagnostics)
        @test hasproperty(res.meson_results[:eta].root_diagnostics, :selected_method)
        @test hasproperty(res.meson_results[:eta].root_diagnostics, :attempts)
        @test !isempty(res.meson_results[:eta].root_diagnostics.attempts)

        if res.meson_results[:eta].root_quality != :good
            @test !res.meson_results[:eta].converged || res.meson_results[:eta].residual > 1e-6
        end
        if res.meson_results[:eta_prime].root_quality != :good
            @test !res.meson_results[:eta_prime].converged || res.meson_results[:eta_prime].residual > 1e-6
        end
    end
end
