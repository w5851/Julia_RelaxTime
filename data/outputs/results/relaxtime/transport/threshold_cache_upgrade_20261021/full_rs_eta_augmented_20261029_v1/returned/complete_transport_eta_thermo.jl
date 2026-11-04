#!/usr/bin/env julia
# Replay only missing EOS columns from saved equilibria, on quote-remote.
# Arguments: FROZEN_SOURCE_ROOT AUGMENTED_CASE ETA_CASE
using CSV, JSON3, SHA, Dates, StaticArrays, LinearAlgebra, Sockets

length(ARGS) == 3 || error("expected FROZEN_SOURCE_ROOT AUGMENTED_CASE ETA_CASE")
const SOURCE_ROOT, CASE, ETA_CASE = abspath.(ARGS)
const OUTPUT = joinpath(CASE, "thermo")
ispath(OUTPUT) && error("refusing to overwrite thermodynamic replay")
digest(path) = bytes2hex(sha256(read(path)))
const EXPECTED_SHA = "2e5c411c716d5cbc524c3fc6cfb335fe05b39444"
const INPUT = JSON3.read(read(joinpath(CASE, "input_manifest.json"), String))
String(INPUT.calculation_sha) == EXPECTED_SHA || error("calculation source mismatch")
for item in INPUT.files
    path = joinpath(CASE, String(item.path))
    filesize(path) == item.bytes && digest(path) == String(item.sha256) || error("changed input: $path")
end
const SOURCE_HASHES = JSON3.read(read(joinpath(ETA_CASE, "expected_source_hashes.json"), String))
for (relative, expected) in pairs(SOURCE_HASHES)
    path = joinpath(SOURCE_ROOT, String(relative))
    digest(path) == String(expected) || error("changed frozen source: $relative")
end
include(joinpath(SOURCE_ROOT, "src", "models", "Models.jl"))
include(joinpath(SOURCE_ROOT, "scripts", "relaxtime", "scan_quality.jl"))
const SP = Models.TransportSamplingProfiles
const RT = Main.RelaxTime.RelaxationTime
const SPECIES = (:u, :d, :s, :ubar, :dbar, :sbar)
const HC = Main.Constants_PNJL.ħc_MeV_fm
LinearAlgebra.BLAS.set_num_threads(1)

function checked_relative(actual, expected, name; tolerance=1e-10)
    isfinite(actual) && isfinite(expected) && expected != 0 || error("invalid $name")
    difference = abs(actual / expected - 1)
    difference <= tolerance || error("$name mismatch: relative difference $difference")
    difference
end

function main()
    selection = JSON3.read(read(joinpath(CASE, "selection", "selection.json"), String))
    point_path, rate_path = (joinpath(CASE, "selection", f) for f in ("points.csv", "channels.csv"))
    digest(point_path) == String(selection.points_sha256) || error("selected points changed")
    digest(rate_path) == String(selection.channels_sha256) || error("selected rates changed")
    points, rates = collect(CSV.File(point_path)), collect(CSV.File(rate_path))
    length(points) == 57 && length(unique(p.point_id for p in points)) == 57 || error("invalid selected points")
    profile = SP.load_profile(joinpath(ETA_CASE, "inputs_v2", "fine.toml"))
    model = SP.profile_model(profile)
    p_nodes, cos_nodes = profile.background.p_nodes, profile.background.cos_nodes
    rows, checks = NamedTuple[], NamedTuple[]
    for point in points
        x = SVector{5}(Float64(point.phi_u), Float64(point.phi_d), Float64(point.phi_s),
                       Float64(point.Phi), Float64(point.Phibar))
        T, mu, xi = Float64(point.T_fm), Float64(point.muq_fm), Float64(point.xi)
        mu_vec = SVector{3}(mu, mu, mu)
        abs(T * HC - Float64(point.T_MeV)) < 1e-8 || error("temperature units mismatch")
        abs(3mu * HC - Float64(point.muB_MeV)) < 1e-8 || error("chemical potential units mismatch")
        masses = Models.calculate_mass_vec(model, SVector{3}(x[1], x[2], x[3]))
        mass_error = maximum(checked_relative(masses[i], Float64(getproperty(point, field)), "mass")
                             for (i, field) in enumerate((:m_u, :m_d, :m_s)))
        pressure, _, entropy, epsilon = Models.model_thermo(model, x, mu_vec, T;
            p_num=p_nodes, t_num=cos_nodes, xi=xi)
        omega = Models.omega(model, x, T, mu_vec; p_num=p_nodes, t_num=cos_nodes, xi=xi)
        entropy_error = checked_relative(entropy, Float64(point.s_fm3inv), "entropy")
        checked_relative(entropy, Float64(point.entropy), "second entropy column")
        omega_error = checked_relative(omega, Float64(point.omega), "omega")
        all(isfinite, (pressure, entropy, epsilon, omega)) || error("nonfinite thermodynamics")
        densities = (; (q => Float64(getproperty(point, Symbol("n_", q))) for q in SPECIES)...)
        tau = (; (q => Float64(getproperty(point, Symbol("tau_", q))) for q in SPECIES)...)
        selected_rates = filter(r -> String(r.point_id) == String(point.point_id), rates)
        length(selected_rates) == 17 || error("incomplete selected channel rates")
        values = Dict(Symbol(r.channel) => Float64(r.rate) for r in selected_rates)
        Set(keys(values)) == Set(RT.REQUIRED_RATE_KEYS_FOR_TAU) || error("wrong channel identities")
        for r in selected_rates
            String(r.variant) == String(point.variant) && String(r.source_run) == String(point.source_run) ||
                error("selected rates use a different precision/source")
            all(Float64(getproperty(r, k)) == Float64(getproperty(point, k)) for k in (:T_MeV, :muB_MeV, :xi)) ||
                error("selected rates use a different coordinate")
        end
        reconstructed = RT.relaxation_rates(densities, values)
        rate_error = maximum(checked_relative(getproperty(reconstructed, q),
            Float64(getproperty(point, Symbol("tauinv_", q))), "inverse relaxation time") for q in SPECIES)
        tau_error = maximum(checked_relative(1 / getproperty(reconstructed, q), getproperty(tau, q),
                                             "relaxation time") for q in SPECIES)
        ratio_error = maximum((
            checked_relative(Float64(point.eta) / entropy, Float64(point.eta_over_s), "eta/s"),
            checked_relative(Float64(point.zeta) / entropy, Float64(point.zeta_over_s), "zeta/s"),
            checked_relative(Float64(point.sigma) / T, Float64(point.sigma_over_T), "sigma/T")))
        flag, reason, metric = ScanQuality.assess_point_quality(tau, Float64(point.eta_over_s), Float64(point.sigma_over_T))
        rho = ((densities.u - densities.ubar) + (densities.d - densities.dbar) + (densities.s - densities.sbar)) / 3
        push!(rows, (point_id=String(point.point_id), xi=xi, T_MeV=Float64(point.T_MeV), muB_MeV=Float64(point.muB_MeV),
            rho_baryon=rho, rho_norm=rho / Main.Constants_PNJL.ρ0_inv_fm3,
            omega_fm4inv=omega, P_fm4inv=pressure, epsilon_fm4inv=epsilon,
            omega_MeV_fm3=omega * HC, P_MeV_fm3=pressure * HC, epsilon_MeV_fm3=epsilon * HC,
            eps_minus_3P_over_T4=(epsilon - 3pressure) / T^4,
            sigma_over_T_over_eta_over_s=Float64(point.sigma_over_T) / Float64(point.eta_over_s),
            zeta_over_s_over_eta_over_s=Float64(point.zeta_over_s) / Float64(point.eta_over_s),
            quality_flag=flag, quality_reason=reason, quality_metric=metric))
        push!(checks, (point_id=String(point.point_id), mass_relative_error=mass_error,
            entropy_relative_error=entropy_error, omega_relative_error=omega_error,
            rates_relative_error=rate_error, tau_relative_error=tau_error, ratios_relative_error=ratio_error))
    end
    mkpath(OUTPUT)
    CSV.write(joinpath(OUTPUT, "rows.csv"), rows)
    CSV.write(joinpath(OUTPUT, "replay_checks.csv"), checks)
    maxima = Dict(String(k) => maximum(getproperty(r, k) for r in checks) for k in keys(first(checks)) if k != :point_id)
    summary = (status="passed", points=length(rows), channels=length(rates),
        rows_sha256=digest(joinpath(OUTPUT, "rows.csv")), checks_sha256=digest(joinpath(OUTPUT, "replay_checks.csv")),
        selection_sha256=digest(joinpath(CASE, "selection", "selection.json")),
        calculation_sha=EXPECTED_SHA, source_root=SOURCE_ROOT, source_files_checked=length(SOURCE_HASHES),
        source_manifest_sha256=digest(joinpath(ETA_CASE, "expected_source_hashes.json")),
        script_sha256=digest(@__FILE__), julia_version=string(VERSION), hostname=gethostname(),
        completed_utc=string(now(UTC)), tolerance=1e-10, maximum_relative_errors=maxima,
        effective_thermo=(p_nodes=p_nodes, cos_nodes=cos_nodes,
            thermal_p_max_inv_fm=profile.background.thermal_p_max_inv_fm),
        original_solver=(residual_norm_max=1e-10, max_iter=80, stable_multiseed=true),
        plotted_values_recomputed=false, gap_solver_called=false, scattering_integration_called=false,
        scope="EOS replay from saved condensates and Polyakov fields; mass/entropy/omega and rate-to-tau/ratio reconstruction checks")
    open(joinpath(OUTPUT, "summary.json"), "w") do io
        JSON3.pretty(io, summary)
        println(io)
    end
    println(JSON3.write(summary))
end
main()
