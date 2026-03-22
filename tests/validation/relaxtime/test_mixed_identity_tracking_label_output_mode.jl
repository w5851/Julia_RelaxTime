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

@testset "Identity tracking keeps label-facing output" begin
    muq_fm = (600.0 / 3.0) / _HBARC_MEV_FM
    T_fm = 220.0 / _HBARC_MEV_FM

    res = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        muq_fm;
        xi=0.0,
        mesons=(:eta, :eta_prime),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=40,),
        mixed_branch_align=:identity_track_label_output,
    )

    @test hasproperty(res, :mixed_seed_tracking)
    @test res.meson_results[:eta].branch_sign == -1
    @test res.meson_results[:eta_prime].branch_sign == 1
end
