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

@testset "Identity alignment lowers mixed-branch selected score" begin
    muq_fm = (600.0 / 3.0) / _HBARC_MEV_FM

    label_scores = Float64[]
    id_scores = Float64[]

    seed_state_label = nothing
    seed_state_id = nothing

    for T_MeV in (220.0, 240.0)
        T_fm = T_MeV / _HBARC_MEV_FM

        res_label = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=40,),
            meson_seed_state=seed_state_label,
            mixed_branch_align=:label,
        )
        seed_state_label = res_label.meson_seed_state

        res_id = _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:eta, :eta_prime),
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=40,),
            meson_seed_state=seed_state_id,
            mixed_branch_align=:identity,
        )
        seed_state_id = res_id.meson_seed_state

        push!(label_scores, res_label.meson_results[:eta].branch_score_selected)
        push!(label_scores, res_label.meson_results[:eta_prime].branch_score_selected)
        push!(id_scores, res_id.meson_results[:eta].branch_score_selected)
        push!(id_scores, res_id.meson_results[:eta_prime].branch_score_selected)
    end

    @test sum(id_scores) <= sum(label_scores)
end

@testset "Identity-track-label-output keeps label branch semantics" begin
    muq_fm = (600.0 / 3.0) / _HBARC_MEV_FM

    seed_state = nothing
    tracking_state = nothing
    for T_MeV in (220.0, 240.0)
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
            mixed_seed_tracking_state=tracking_state,
            mixed_branch_align=:identity_track_label_output,
        )

        seed_state = res.meson_seed_state
        tracking_state = res.mixed_seed_tracking

        @test res.meson_results[:eta].branch_sign == -1
        @test res.meson_results[:eta_prime].branch_sign == 1
        @test tracking_state !== nothing
        @test haskey(tracking_state, :eta_minus)
        @test haskey(tracking_state, :eta_plus)
    end
end

@testset "strict_sign_binding keeps fixed mixed sign mapping" begin
    muq_fm = (600.0 / 3.0) / _HBARC_MEV_FM

    seed_state = nothing
    tracking_state = nothing
    for T_MeV in (220.0, 240.0)
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
            mixed_seed_tracking_state=tracking_state,
            mixed_branch_align=:strict_sign_binding,
        )

        seed_state = res.meson_seed_state
        tracking_state = res.mixed_seed_tracking

        @test res.meson_results[:eta].branch_sign == -1
        @test res.meson_results[:eta_prime].branch_sign == 1
        @test tracking_state !== nothing
        @test haskey(tracking_state, :eta_minus)
        @test haskey(tracking_state, :eta_plus)
    end
end
