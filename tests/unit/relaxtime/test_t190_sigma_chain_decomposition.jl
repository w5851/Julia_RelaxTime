using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "analysis", "relaxtime", "t190_sigma_chain_decomposition_lib.jl"))

@testset "T190 sigma amplitude decomposition" begin
    st = build_state_point(190.0, 0.0, -0.10)
    th = process_threshold_info(:uubar_to_ddbar, st.quark_params)
    s = th.s_th + 0.2
    tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
    t = 0.5 * (tb.t_min + tb.t_max)

    parts = decompose_qqbar_amplitude_terms(:uubar_to_ddbar, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
    @test isfinite(parts.M2_total)
    @test isapprox(parts.M2_total, parts.M_s_sq + parts.M_t_sq + parts.M_interf; rtol=1e-9, atol=1e-11)
    @test isapprox(parts.M_s_sq, parts.M_s_S + parts.M_s_P; rtol=1e-9, atol=1e-11)
    @test isapprox(parts.M_t_sq, parts.M_t_S + parts.M_t_P; rtol=1e-9, atol=1e-11)
end

@testset "T190 s-channel factor decomposition" begin
    st = build_state_point(190.0, 0.0, -0.10)
    th = process_threshold_info(:uubar_to_uubar, st.quark_params)
    s = th.s_th + 0.25
    tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
    t = 0.5 * (tb.t_min + tb.t_max)

    fac = decompose_qqbar_s_channel_factors(:uubar_to_uubar, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
    @test isfinite(fac.M_s_S)
    @test isfinite(fac.M_s_P)
    @test isapprox(fac.M_s_S, fac.abs_D_s_S_sq * fac.kin_s_S; rtol=1e-9, atol=1e-11)
    @test isapprox(fac.M_s_P, fac.abs_D_s_P_sq * fac.kin_s_P; rtol=1e-9, atol=1e-11)
end

@testset "T190 P-channel propagator absolute strength decomposition" begin
    st = build_state_point(190.0, 0.0, -0.10)
    for process in (:uubar_to_ddbar, :uubar_to_uubar)
        th = process_threshold_info(process, st.quark_params)
        s = th.s_th + 0.2
        tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
        t = 0.5 * (tb.t_min + tb.t_max)

        p = decompose_p_channel_propagator_strength(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)

        @test isfinite(p.abs_D_s_P_sq_total)
        @test p.abs_D_s_P_sq_total >= 0.0
        @test isapprox(p.D_s_P_total, p.D_s_P_simple + p.D_s_P_mixed; rtol=1e-9, atol=1e-11)
        @test isapprox(p.abs_D_s_P_sq_total, abs2(p.D_s_P_total); rtol=1e-9, atol=1e-11)
        @test isapprox(p.M_s_P_total, p.abs_D_s_P_sq_total * p.kin_s_P; rtol=1e-9, atol=1e-11)
    end
end

@testset "T190 mixed_P propagator chain decomposition" begin
    st = build_state_point(190.0, 0.0, -0.10)
    for process in (:uubar_to_ddbar, :uubar_to_uubar)
        th = process_threshold_info(process, st.quark_params)
        s = th.s_th + 0.2
        tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
        t = 0.5 * (tb.t_min + tb.t_max)

        d = decompose_mixed_p_propagator_chain(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
        @test isfinite(d.detK_plus)
        @test isfinite(d.abs_detM_sq)
        @test d.abs_detM_sq > 0.0
        @test isapprox(d.M00, d.M00_from_K0 + d.M00_from_Piuu + d.M00_from_Piss; rtol=1e-9, atol=1e-11)
        @test isapprox(d.M08, d.M08_from_K08 + d.M08_from_Piuu + d.M08_from_Piss; rtol=1e-9, atol=1e-11)
        @test isapprox(d.M88, d.M88_from_K8 + d.M88_from_Piuu + d.M88_from_Piss; rtol=1e-9, atol=1e-11)
        @test isapprox(d.M00_from_K0, st.K_coeffs.K0_plus; rtol=1e-9, atol=1e-11)
        @test isapprox(d.M08_from_K08, st.K_coeffs.K08_plus; rtol=1e-9, atol=1e-11)
        @test isapprox(d.M88_from_K8, st.K_coeffs.K8_plus; rtol=1e-9, atol=1e-11)
        @test isapprox(d.detM_term_M00M88, d.M00 * d.M88; rtol=1e-9, atol=1e-11)
        @test isapprox(d.detM_term_M08sq, d.M08 * d.M08; rtol=1e-9, atol=1e-11)
        @test isapprox(d.detM, d.detM_term_M00M88 - d.detM_term_M08sq; rtol=1e-9, atol=1e-11)
        @test isapprox(d.detM_cross_term, -2.0 * real(d.detM_term_M00M88 * conj(d.detM_term_M08sq)); rtol=1e-9, atol=1e-11)
        @test isapprox(d.abs_detM_sq, d.abs_detM_term_M00M88_sq + d.abs_detM_term_M08sq_sq + d.detM_cross_term; rtol=1e-9, atol=1e-11)
        @test isapprox(d.D_mixed_P, d.prefactor * d.JMJ; rtol=1e-9, atol=1e-11)
        @test isapprox(d.abs_D_mixed_P_sq, abs2(d.D_mixed_P); rtol=1e-9, atol=1e-11)
        @test isapprox(d.abs_D_mixed_P_sq, 4.0 * d.detK_plus^2 * d.abs_JMJ_sq / d.abs_detM_sq; rtol=1e-9, atol=1e-11)
    end
end
