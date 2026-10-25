using Test, TOML
if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end
const TSP = Models.TransportSamplingProfiles

@testset "transport sampling profile effective numerical context" begin
    profile = TSP.load_profile("thermal_weighted_v1")
    runtime = TSP.build_runtime(profile)
    baseline = Models.create_model(:PNJL)
    @test runtime.model.params.thermal_p_max_inv_fm == 25.0
    @test all(k -> k === :thermal_p_max_inv_fm || isequal(getfield(runtime.model.params,k),
        getfield(baseline.params,k)), fieldnames(typeof(baseline.params)))
    @test runtime.workflow_cache.model_cache[:PNJL] === runtime.model
    @test length(runtime.p_grid) == 48
    @test length(runtime.cos_grid) == length(runtime.phi_grid) == 36
    @test sum(runtime.cos_w) ≈ 2.0
    @test sum(runtime.cos_w .* runtime.cos_grid.^2) ≈ 2/3
    @test sum(runtime.phi_w) ≈ 2pi
    @test runtime.tau_kwargs.n_sigma_points == 24
    @test runtime.tau_kwargs.sampling_config.rtol == 1e-3
    @test runtime.tau_kwargs.sampling_config.relative_velocity === :cm
    @test TSP.effective_config(profile).sampling.relative_velocity === :cm
    @test TSP.effective_config(profile).a_polarization == "full_RS"
    @test TSP.effective_config(profile).a_coupling == "full_RS"
    @test TSP.effective_config(profile).internal_density.cos_nodes == 36
    p,w = TSP.AS._cached_semi_infinite_momentum_grid(profile.terminal_integral.p_nodes,
        profile.terminal_integral.scale_inv_fm)
    @test sum(w .* exp.(-p)) ≈ 1.0 rtol=1e-6
    @test sum(w .* p.^2 .* exp.(-p)) ≈ 2.0 rtol=1e-5
end

@testset "profile uses full-RS A and coupling conventions" begin
    profile = TSP.load_profile("thermal_weighted_v1")
    eq = (x_state=[-0.1,-0.1,-0.2,0.4,0.3],mu_vec=[0.2,0.2,0.2],masses=[1.0,1.0,2.0])
    matched = TSP.point_state(profile,eq,1.0,0.2;xi=0.37)
    isotropic = TSP.point_state(profile,eq,1.0,0.2;xi=0.37,propagator_xi_policy=:isotropic)
    @test matched.thermo_params == isotropic.thermo_params
    @test matched.quark_params == isotropic.quark_params
    full_k = Main.RelaxTime.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,
        Main.RelaxTime.EffectiveCouplings.calculate_G_from_A(matched.quark_params.A.u,eq.masses[1]),
        Main.RelaxTime.EffectiveCouplings.calculate_G_from_A(matched.quark_params.A.s,eq.masses[3]))
    isotropic_k = Main.RelaxTime.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,
        Main.RelaxTime.EffectiveCouplings.calculate_G_from_A(isotropic.propagator_quark_params.A.u,eq.masses[1]),
        Main.RelaxTime.EffectiveCouplings.calculate_G_from_A(isotropic.propagator_quark_params.A.s,eq.masses[3]))
    @test matched.K_coeffs == full_k
    @test matched.K_coeffs != isotropic_k
    @test matched.propagator_quark_params.A != isotropic.propagator_quark_params.A
    zero_a = TSP.point_state(profile,eq,1.0,0.2;xi=0.0)
    zero_b = TSP.point_state(profile,eq,1.0,0.2;xi=0.0,propagator_xi_policy=:isotropic)
    @test all(isapprox(zero_a.propagator_quark_params.A[k],zero_b.propagator_quark_params.A[k];
        rtol=1e-12) for k in (:u,:d,:s))
    @test all(isapprox(zero_a.K_coeffs[k],zero_b.K_coeffs[k]; rtol=1e-12) for k in keys(zero_a.K_coeffs))
    @test_throws ArgumentError TSP.point_state(profile,merge(eq,(masses=[NaN,1.0,2.0],)),1.0,0.2)
end

@testset "invalid sampling profiles fail before a numerical run" begin
    source = TOML.parsefile(TSP.load_profile("thermal_weighted_v1").source_path)
    mktempdir() do dir
        for (section,key,value) in (("tau_integral","p_nodes",0),("background","thermal_p_max_inv_fm",Inf),
            ("a_integral","ignored_knob",3),("sampling","max_evaluations",0))
            data=deepcopy(source)
            data[section][key]=value
            path=joinpath(dir,"invalid.toml")
            open(path,"w") do io
                TOML.print(io,data)
            end
            @test_throws ArgumentError TSP.load_profile(path)
        end
    end
end
