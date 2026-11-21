using Test, ForwardDiff, TOML
if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "Ward A and direction-averaged linear loop" begin
    rt = Main.RelaxTime
    cor = rt.OneLoopIntegralsCorrection
    gl = rt.GaussLegendre
    p, pw = gl.gauleg(0.0, 25.0, 96)
    c, cw = gl.gauleg(-1.0, 1.0, 24)
    args = (0.7, 0.2, 0.9, 0.4, 0.55)
    @test cor.A_aniso_ward(args..., 0.0, p,pw,c,cw) ≈
        rt.OneLoopIntegrals.A(args...,p,pw) rtol=2e-14
    # Direct derivative of the thermal log integrand has M/Eξ.
    for xi in (-0.5, 0.3, 0.5)
        value = cor.A_aniso_ward(args...,xi,p,pw,c,cw)
        reference = 0.0
        for (v,w) in zip(p,pw), (x,wx) in zip(c,cw)
            E = sqrt(v*v*(1+xi*x*x)+args[1]^2)
            fp = rt.OneLoopIntegrals.distribution_value(:pnjl,:plus,E,args[2:5]...)
            fm = rt.OneLoopIntegrals.distribution_value(:pnjl,:minus,E,args[2:5]...)
            reference += w*wx*v*v*(fp+fm)/E
        end
        @test value ≈ 2reference-4rt.OneLoopIntegrals.const_integral_term_A(args[1]) rtol=1e-14
    end
    @test_throws ArgumentError cor.A_aniso_ward(args...,-1.0,p,pw,c,cw)
    for (a,b) in ((0.8,1.4),(0.8,-1.4),(1e-7,1.4))
        numerical = sum(w/(3*(a*x+b)) for (x,w) in zip(c,cw))
        @test cor._real_direction_kernel(a,b,:sphere_average) ≈ numerical rtol=2e-13
        # At fixed x, the azimuthal average is x² cos²β+(1-x²)sin²β/2.
        for x in (-0.8,0.0,0.3,1.0)
            direction_mean = sum(w*(x*x*y*y+(1-x*x)*(1-y*y)/2)/2 for (y,w) in zip(c,cw))
            @test direction_mean ≈ 1/3 rtol=2e-14
        end
    end
    @test cor._imag_direction_kernel(0.8,0.3,:sphere_average) ≈ pi/(3*0.8)
    bargs = (1.4,0.4,0.7,1.1,0.2,0.2,0.9,0.4,0.55,0.0)
    @test cor.B0_correction(bargs...;b0_direction=:aligned) ==
        cor.B0_correction(bargs...;b0_direction=:sphere_average)
    for xi in (-0.5,0.5)
        rest = (1.4,0.0,0.7,1.1,0.2,0.2,0.9,0.4,0.55,xi)
        @test cor.B0_correction(rest...;b0_direction=:aligned) ==
            cor.B0_correction(rest...;b0_direction=:sphere_average)
        base = cor.real_integrand_k_zero(:quark,2.0,0.7,1.1,1.4,xi,0.9,0.2,0.4,0.55)
        near = cor.real_integrand_k_positive(:quark,2.0,1e-5,0.7,1.1,1.4,xi,0.9,0.2,0.4,0.55;
            b0_direction=:sphere_average)
        @test near ≈ base rtol=1e-9
    end
    @test_throws ArgumentError cor.B0_correction(bargs...;b0_direction=:unknown)
    # The new angular kernel remains differentiable in the real external energy.
    derivative = ForwardDiff.derivative(b->cor._real_direction_kernel(0.8,b,:sphere_average),1.4)
    h = 1e-5
    finite_difference = (cor._real_direction_kernel(0.8,1.4+h,:sphere_average)-
        cor._real_direction_kernel(0.8,1.4-h,:sphere_average))/(2h)
    @test derivative ≈ finite_difference rtol=1e-8
end

@testset "prescription identities and actual propagation" begin
    sp = Models.TransportSamplingProfiles
    rt = Main.RelaxTime
    aligned = sp.load_profile("ward_moller_aligned_v1")
    sphere = sp.load_profile("ward_moller_sphere_v1")
    legacy = sp.load_profile("thermal_weighted_v1")
    @test aligned.sampling.relative_velocity === sphere.sampling.relative_velocity === :medium_moller
    @test aligned.background == sphere.background == legacy.background
    @test aligned.a_integral == sphere.a_integral == legacy.a_integral
    @test aligned.tau_integral == sphere.tau_integral == legacy.tau_integral
    @test aligned.terminal_integral == sphere.terminal_integral == legacy.terminal_integral
    @test sp.effective_config(sphere).a_coupling == "ward_matched_RS"
    @test sp.effective_config(sphere).b0_direction == "sphere_average"
    qp = (m=(u=0.7,d=0.7,s=1.1),μ=(u=0.2,d=0.2,s=0.2))
    tp = (T=0.9,Φ=0.4,Φbar=0.55,ξ=0.3)
    st1, st2 = sp.fixed_point_state(aligned,qp,tp), sp.fixed_point_state(sphere,qp,tp)
    @test st1.quark_params == st2.quark_params
    @test st1.K_coeffs == st2.K_coeffs
    pc = rt.PolarizationCache
    pc.reset_cache!()
    pargs = (:P,1.4,0.4,0.7,0.7,0.2,0.2,0.9,0.4,0.55,0.3,
        st1.quark_params.A.u,st1.quark_params.A.u,0)
    a = pc.polarization_aniso_cached(pargs...;b0_direction=:aligned)
    b = pc.polarization_aniso_cached(pargs...;b0_direction=:sphere_average)
    @test a != b
    @test pc.get_cache_stats().cache_misses == 2
    @test pc.polarization_aniso_cached(pargs...;b0_direction=:aligned) == a
    @test pc.get_cache_stats().cache_hits == 1
    options = (p_nodes=4,angle_nodes=3,phi_nodes=5,n_sigma_points=4)
    as = rt.AverageScatteringRate
    prob1 = as._prepare_thermal_cache_problem(:us_to_us,st1.quark_params,st1.thermo_params,st1.K_coeffs;options...)
    prob2 = as._prepare_thermal_cache_problem(:us_to_us,st2.quark_params,st2.thermo_params,st2.K_coeffs;options...)
    @test prob1.identity != prob2.identity
    # Check that the thermo field reaches the public scattering entrypoint.
    s = (qp.m.u+qp.m.s+0.5)^2
    cs1 = rt.TotalCrossSection.total_cross_section(:us_to_us,s,st1.quark_params,st1.thermo_params,st1.K_coeffs;n_points=8)
    cs2 = rt.TotalCrossSection.total_cross_section(:us_to_us,s,st2.quark_params,st2.thermo_params,st2.K_coeffs;n_points=8)
    @test isfinite(cs1) && isfinite(cs2) && cs1 > 0 && cs2 > 0
    @test !isapprox(cs1,cs2;rtol=1e-6)
    # A frozen endpoint from the published high-density background. The soft
    # pion relation includes the explicit current-quark mass m0/Mu.
    dense_q = (m=(u=0.10385714343595118,d=0.1038571434359508,s=2.1679484929487445),
        μ=(u=1.5203192153038187,d=1.5203192153038187,s=1.5203192153038187))
    dense_t = (T=0.6371961908588305,Φ=0.3818517764595931,Φbar=0.40758787988004697,ξ=-0.5)
    dense = sp.fixed_point_state(aligned,dense_q,dense_t)
    constants = Main.Constants_PNJL
    residual = 1 + constants.N_color/(2pi^2)*dense.K_coeffs.K123_plus*
        (dense.quark_params.A.u+dense.quark_params.A.d) - constants.m_ud0_inv_fm/dense_q.m.u
    @test abs(residual) < 1e-8
    stock = sp.fixed_point_state(legacy,dense_q,dense_t)
    stock_residual = 1 + constants.N_color/(2pi^2)*stock.K_coeffs.K123_plus*
        (stock.quark_params.A.u+stock.quark_params.A.d) - constants.m_ud0_inv_fm/dense_q.m.u
    @test stock_residual ≈ -0.0736026 atol=1e-6
end
