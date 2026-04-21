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

const _MESON_WORKFLOW_MOD = Models.meson_workflow_module()
const _HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm

@testset "Fortran-track experiment path remains non-trivial" begin
    muB_MeV = 0.0
    muq_fm = (muB_MeV / 3.0) / _HBARC_MEV_FM

    seed_default = nothing
    seed_track = nothing

    baseline_bad = Dict{Float64,Float64}()
    track_resid = Dict{Float64,Float64}()

    for T_MeV in (200.0, 220.0, 240.0)
        T_fm = T_MeV / _HBARC_MEV_FM

        res_default = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=40,),
            meson_seed_state=seed_default,
            meson_root_policy=(residual_norm_max=1e-6, require_converged=true, use_trust_region_fallback=true),
        )
        seed_default = res_default.meson_seed_state

        res_track = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=40,),
            meson_seed_state=seed_track,
            meson_root_policy=(residual_norm_max=1e-6, require_converged=true, use_trust_region_fallback=false),
            meson_experiment_mode=:fortran_track,
        )
        seed_track = res_track.meson_seed_state

        if T_MeV in (220.0, 240.0)
            baseline_bad[T_MeV] = res_default.meson_results[:eta_prime].residual
            track_resid[T_MeV] = res_track.meson_results[:eta_prime].residual
            @test isfinite(track_resid[T_MeV])
        end
    end

    @test baseline_bad[220.0] > 1e-3
    @test baseline_bad[240.0] > 1e-3
    @test track_resid[220.0] > 1e-3
    @test track_resid[240.0] > 1e-3
end
