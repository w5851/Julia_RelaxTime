using Test

const _FIXED_MUB_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(_FIXED_MUB_PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "FixedMuBConservedCharges algebra and residual contract" begin
    mode = Models.FixedMuBConservedCharges(1.2; charge_to_baryon_ratio=0.4)
    @test mode.muB_fm == 1.2
    @test mode.charge_to_baryon_ratio == 0.4
    @test mode.strangeness_density_target == 0.0
    @test Models.state_dim(mode) == 8
    @test Models.state_var_dim(mode) == 5
    @test Models.mu_var_dim(mode) == 3
    @test Models.solution_dim(mode) == 8
    @test Models.param_dim(mode) == 1
    @test occursin("rho_Q/rho_B=0.4", Models.constraint_description(mode))
    @test occursin("FixedMuBConservedCharges", sprint(show, mode))

    @test_throws ArgumentError Models.FixedMuBConservedCharges(NaN)
    @test_throws ArgumentError Models.FixedMuBConservedCharges(1.0, Inf, 0.0)
    @test_throws ArgumentError Models.FixedMuBConservedCharges(1.0, 0.4, NaN)

    conserved = (mu_B=1.2, mu_Q=-0.08, mu_S=0.31)
    flavor = Models.flavor_mu_from_bqs(conserved...)
    roundtrip = Models.conserved_mu_from_flavor(flavor...)
    @test roundtrip.mu_B ≈ conserved.mu_B atol=1e-14
    @test roundtrip.mu_Q ≈ conserved.mu_Q atol=1e-14
    @test roundtrip.mu_S ≈ conserved.mu_S atol=1e-14
    scaled = Models.flavor_mu_from_bqs(2 .* collect(conserved)...)
    @test collect(scaled) ≈ 2 .* collect(flavor) atol=1e-14

    densities = Models.conserved_densities_from_flavor([3.0, 2.0, 0.5])
    @test densities.rho_B ≈ 5.5 / 3
    @test densities.rho_Q ≈ 3.5 / 3
    @test densities.rho_S == -0.5
    @test_throws ArgumentError Models.conserved_densities_from_flavor([1.0, 2.0])

    components = Models.build_constraint_components(mode)
    @test Models.constraint_name.(components) == [
        :stationarity,
        :fixed_muB,
        :conserved_charge_density,
    ]
    @test Models.constraint_total_dim(components) == Models.state_dim(mode)

    T_fm = 0.7
    params = Models.GapParams(
        T_fm,
        Models.cached_nodes(8, 4),
        0.0;
        p_num=8,
        t_num=4,
        model_kind=:PNJL,
    )
    x = [-1.5, -1.5, -2.1, 0.2, 0.2, flavor.mu_u, flavor.mu_d, flavor.mu_s]
    conditions = Models.build_conditions(mode, params)
    expected = conditions([T_fm], x)
    residual! = Models.build_residual!(mode, params)
    actual = zeros(8)
    residual!(actual, x)
    explicit = Models.explicit_residual(mode, x, [T_fm], params)

    @test actual ≈ expected rtol=1e-12 atol=1e-12
    @test explicit ≈ expected rtol=1e-12 atol=1e-12
    @test abs(actual[6]) < 1e-14

    spec = Models.build_problem_spec(mode)
    @test spec.mode === mode
    @test spec.x_dim == 8
    @test spec.theta_dim == 1

    model = Models.create_model(:PNJL)
    synthetic_spec = Models.ProblemSpec(
        mode;
        x_dim=8,
        theta_dim=1,
        forward_solve=(model_arg, T_arg; kwargs...) -> (
            converged=true,
            solution=copy(x),
            x_state=copy(x[1:5]),
            mu_vec=copy(x[6:8]),
            omega=-1.0,
            pressure=1.0,
            rho_norm=0.2,
            entropy=0.3,
            energy=0.4,
            masses=[1.0, 1.0, 2.0],
            iterations=1,
            residual_norm=0.0,
        ),
    )
    vec_result = Models.solve_vec(model, mode, [T_fm]; problem_spec=synthetic_spec)
    named_result = Models.solve_named(model, mode, (T_fm=T_fm,); problem_spec=synthetic_spec)
    @test vec_result isa Models.SolverResult
    @test named_result isa Models.SolverResult
    @test vec_result.solution == named_result.solution == x
    @test_throws ArgumentError Models.solve_vec(model, mode, [T_fm, mode.muB_fm]; problem_spec=synthetic_spec)
    @test_throws ArgumentError Models.solve_named(model, mode, (temperature=T_fm,); problem_spec=synthetic_spec)

    seed = copy(Models.pnjl_module().HADRON_SEED_8)
    @test_throws ArgumentError spec.forward_solve(
        model,
        T_fm;
        seed_guess=seed,
        rho0=0.0,
        p_num=8,
        t_num=4,
        nlsolve_method=:trust_region,
    )
    @test_throws ArgumentError spec.forward_solve(
        model,
        T_fm;
        seed_guess=seed,
        rho0=0.16,
        p_num=8,
        t_num=4,
        nlsolve_method=:unsupported,
        trust_region_fallback=false,
    )
end
