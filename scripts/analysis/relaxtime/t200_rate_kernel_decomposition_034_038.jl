#!/usr/bin/env julia

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function main()
    process = :uubar_to_ddbar
    xis = [0.34, 0.36, 0.38]
    T_MeV = 200.0
    muB_MeV = 0.0

    p_nodes = 20
    angle_nodes = 4
    phi_nodes = 8
    n_sigma_points = 6
    sigma_grid_n = Main.AverageScatteringRate.DEFAULT_SIGMA_GRID_N
    design_p_nodes = Main.AverageScatteringRate.DEFAULT_W0CDF_P_NODES
    design_angle_nodes = Main.AverageScatteringRate.DEFAULT_W0CDF_ANGLE_NODES
    design_phi_nodes = Main.AverageScatteringRate.DEFAULT_W0CDF_PHI_NODES

    p_grid, p_w = Main.GaussLegendre.gauleg(0.0, 15.0, p_nodes)
    cos_grid, cos_w = Main.GaussLegendre.gauleg(-1.0, 1.0, angle_nodes)
    phi_grid, phi_w = Main.GaussLegendre.gauleg(0.0, 2 * pi, phi_nodes)

    out_csv = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_rate_kernel_decomposition_034_038.csv"
    ensure_parent_dir(out_csv)

    io = open(out_csv, "w")
    try
        println(io, "xi,rate_direct,rate_func,rate_recomputed,prefactor,rho_i,rho_j,omega_total,omega_over_rho2,F_sum,F_mean_v,FV_mean_sigma,FV_mean_s,FV_std_s,FV_near_frac_ds_0p02,FV_near_frac_ds_0p05,FV_near_frac_ds_0p10")

        rows = collect(CSV.File(raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_channel_diag.csv"; comment="#"))
        direct_rate = Dict{Float64, Float64}()
        for r in rows
            if Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && String(r.species) == "u" && String(r.channel) == "uubar_to_ddbar"
                direct_rate[Float64(r.xi)] = Float64(r.rate)
            end
        end

        for xi in xis
            st = build_state_point(T_MeV, muB_MeV, xi)
            pi_sym, pj_sym, pc_sym, pd_sym = Main.TotalCrossSection.parse_particles_from_process(process)

            mi = Main.AverageScatteringRate.get_mass(pi_sym, st.quark_params)
            mj = Main.AverageScatteringRate.get_mass(pj_sym, st.quark_params)
            mc = Main.AverageScatteringRate.get_mass(pc_sym, st.quark_params)
            md = Main.AverageScatteringRate.get_mass(pd_sym, st.quark_params)

            μi = Main.AverageScatteringRate.get_mu(pi_sym, st.quark_params)
            μj = Main.AverageScatteringRate.get_mu(pj_sym, st.quark_params)

            cache = Main.AverageScatteringRate.build_w0cdf_pchip_cache(
                process,
                st.quark_params,
                st.thermo_params,
                st.K_coeffs;
                N=sigma_grid_n,
                design_p_nodes=design_p_nodes,
                design_angle_nodes=design_angle_nodes,
                design_phi_nodes=design_phi_nodes,
                p_cutoff=Main.Constants_PNJL.Λ_inv_fm,
                n_sigma_points=n_sigma_points,
                threshold_subtraction=true,
                asym_window=0.6,
                asym_fit_min_points=8,
                asym_extra_points=10,
            )

            # 与扫描一致：finite_15 + apply_s_domain_cut
            _, s_bo, s_up = Main.AverageScatteringRate._resolve_sigma_support_bounds(
                mi,
                mj,
                mc,
                md,
                Main.Constants_PNJL.Λ_inv_fm,
                p_grid,
                p_w,
            )

            ρ_i = Main.AverageScatteringRate.number_density(
                pi_sym,
                mi,
                μi,
                st.thermo_params.T,
                st.thermo_params.Φ,
                st.thermo_params.Φbar,
                st.thermo_params.ξ;
                p_grid=p_grid,
                p_w=p_w,
                angle_nodes=angle_nodes,
            )
            ρ_j = Main.AverageScatteringRate.number_density(
                pj_sym,
                mj,
                μj,
                st.thermo_params.T,
                st.thermo_params.Φ,
                st.thermo_params.Φbar,
                st.thermo_params.ξ;
                p_grid=p_grid,
                p_w=p_w,
                angle_nodes=angle_nodes,
            )

            F_sum = 0.0
            FV_sum = 0.0
            FVsigma_sum = 0.0
            FV_s_sum = 0.0
            FV_s2_sum = 0.0
            FV_near_002 = 0.0
            FV_near_005 = 0.0
            FV_near_010 = 0.0

            s_th = max((mi + mj)^2, (mc + md)^2)

            for (p_i, w_pi) in zip(p_grid, p_w)
                Ei = Main.AverageScatteringRate.energy_from_p(p_i, mi)
                for (p_j, w_pj) in zip(p_grid, p_w)
                    Ej = Main.AverageScatteringRate.energy_from_p(p_j, mj)
                    for (cθi, w_cθi) in zip(cos_grid, cos_w)
                        sθi = sqrt(max(1.0 - cθi * cθi, 0.0))
                        f_i = Main.AverageScatteringRate.distribution_with_anisotropy(
                            pi_sym,
                            p_i,
                            mi,
                            μi,
                            st.thermo_params.T,
                            st.thermo_params.Φ,
                            st.thermo_params.Φbar,
                            st.thermo_params.ξ,
                            cθi,
                        )
                        f_i == 0.0 && continue

                        for (cθj, w_cθj) in zip(cos_grid, cos_w)
                            sθj = sqrt(max(1.0 - cθj * cθj, 0.0))
                            f_j = Main.AverageScatteringRate.distribution_with_anisotropy(
                                pj_sym,
                                p_j,
                                mj,
                                μj,
                                st.thermo_params.T,
                                st.thermo_params.Φ,
                                st.thermo_params.Φbar,
                                st.thermo_params.ξ,
                                cθj,
                            )
                            f_j == 0.0 && continue

                            for (φ, wφ) in zip(phi_grid, phi_w)
                                cosΘ = cθi * cθj + sθi * sθj * cos(φ)
                                s = mi^2 + mj^2 + 2.0 * (Ei * Ej - p_i * p_j * cosΘ)
                                if (s <= s_bo) || (s >= s_up)
                                    continue
                                end

                                s_rt = sqrt(s)
                                Ei_cm = (s + mi^2 - mj^2) / (2.0 * s_rt)
                                Ej_cm = (s - mi^2 + mj^2) / (2.0 * s_rt)
                                pi_cm = sqrt(max(0.0, (s - (mi + mj)^2) * (s - (mi - mj)^2))) / (2.0 * s_rt)
                                pj_cm = pi_cm
                                v_rel_num = (Ei_cm * Ej_cm + pi_cm * pj_cm)^2 - (mi * mj)^2
                                v_rel = v_rel_num > 0.0 ? sqrt(v_rel_num) / (Ei_cm * Ej_cm) : 0.0
                                (v_rel == 0.0 || v_rel > 2.0) && continue

                                σ = Main.AverageScatteringRate._get_sigma_core(
                                    cache,
                                    s,
                                    st.quark_params,
                                    st.thermo_params,
                                    st.K_coeffs;
                                    n_points=n_sigma_points,
                                    interpolation_mode=:linear,
                                )

                                w = w_pi * w_pj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2)
                                F = w * f_i * f_j
                                FV = F * v_rel
                                FVs = FV * σ

                                F_sum += F
                                FV_sum += FV
                                FVsigma_sum += FVs
                                FV_s_sum += FV * s
                                FV_s2_sum += FV * s^2

                                ds = s - s_th
                                if ds <= 0.02
                                    FV_near_002 += FVs
                                end
                                if ds <= 0.05
                                    FV_near_005 += FVs
                                end
                                if ds <= 0.10
                                    FV_near_010 += FVs
                                end
                            end
                        end
                    end
                end
            end

            prefactor = (Main.AverageScatteringRate.DQ^2) / (32.0 * pi^5 * ρ_i * ρ_j)
            rate_recomputed = prefactor * FVsigma_sum

            rate_func = Main.AverageScatteringRate.average_scattering_rate(
                process,
                st.quark_params,
                st.thermo_params,
                st.K_coeffs;
                p_nodes=p_nodes,
                angle_nodes=angle_nodes,
                phi_nodes=phi_nodes,
                p_grid=p_grid,
                p_w=p_w,
                cos_grid=cos_grid,
                cos_w=cos_w,
                phi_grid=phi_grid,
                phi_w=phi_w,
                cs_cache=cache,
                n_sigma_points=n_sigma_points,
                sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
                threshold_subtraction=true,
                asym_window=0.6,
                asym_fit_min_points=8,
                asym_extra_points=10,
                interpolation_mode=:linear,
            )

            F_mean_v = FV_sum / F_sum
            FV_mean_sigma = FVsigma_sum / FV_sum
            FV_mean_s = FV_s_sum / FV_sum
            FV_var_s = max(FV_s2_sum / FV_sum - FV_mean_s^2, 0.0)
            FV_std_s = sqrt(FV_var_s)

            println(io, join((
                xi,
                get(direct_rate, xi, NaN),
                rate_func,
                rate_recomputed,
                prefactor,
                ρ_i,
                ρ_j,
                FVsigma_sum,
                FVsigma_sum / (ρ_i * ρ_j),
                F_sum,
                F_mean_v,
                FV_mean_sigma,
                FV_mean_s,
                FV_std_s,
                FV_near_002 / FVsigma_sum,
                FV_near_005 / FVsigma_sum,
                FV_near_010 / FVsigma_sum,
            ), ','))
        end
    finally
        close(io)
    end
end

main()
