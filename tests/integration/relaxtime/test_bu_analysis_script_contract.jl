using Test

const _BU_ANALYSIS_ROOT = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
))

@testset "BU analysis gate script contracts" begin
    e5_path = joinpath(_BU_ANALYSIS_ROOT, "compare_bu_derivative_vs_byparts_e5.jl")
    charged_path = joinpath(_BU_ANALYSIS_ROOT, "audit_bu_meson_density_literature_alignment.jl")
    helper_path = joinpath(_BU_ANALYSIS_ROOT, "bu_kernel_gate_utils.jl")
    @test all(isfile, (e5_path, charged_path, helper_path))

    e5_source = read(e5_path, String)
    @test occursin("omega_all = vcat(omega_min", e5_source)
    @test occursin("finite_window_bu_identity", e5_source)
    for field in (
        "density_current_boundary_lower",
        "density_current_boundary_upper",
        "density_current_byparts_total",
        "current_closure_abs",
        "density_gbu_boundary_lower",
        "density_gbu_boundary_upper",
        "density_gbu_byparts_total",
        "gbu_closure_abs",
        "bose_omega_min",
        "bose_omega_max",
    )
        @test occursin(field, e5_source)
    end

    charged_source = read(charged_path, String)
    @test occursin("bu2020_mu_s_0p2", charged_source)
    @test occursin("friesen2019_mu_s_0p55", charged_source)
    @test occursin("mesons=(:pi_plus, :K_plus, :pi_minus, :K_minus)", charged_source)
    @test occursin("One and only one upstream state solve per flavor profile", charged_source)
    for field in (
        "equilibrium_group",
        "upstream_state_id",
        "mu_u_MeV",
        "mu_d_MeV",
        "mu_s_MeV",
        "mu_formula",
        "m_pi_minus_mu_pi_MeV",
        "m_K_minus_mu_K_MeV",
        "n_pi_plus",
        "n_K_plus",
        "n_pi_minus",
        "n_K_minus",
        "channel_status",
        "strict_requested_window_status",
    )
        @test occursin(field, charged_source)
    end

    @test Meta.parseall(read(helper_path, String)) !== nothing
    @test Meta.parseall(e5_source) !== nothing
    @test Meta.parseall(charged_source) !== nothing
end
