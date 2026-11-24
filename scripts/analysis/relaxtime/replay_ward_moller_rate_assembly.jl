#!/usr/bin/env julia
# Replay only the species-rate assembly and terminal kernels.
# julia --project=SOURCE replay_ward_moller_rate_assembly.jl SOURCE GRID FIG7 SCAN OUT
# GRID and FIG7 are the frozen pre-fix result folders; SCAN contains the original
# source-indexed nodes and threshold centers listed by the Fig. 7 manifest.
include(joinpath(@__DIR__, "run_ward_moller_mott_response.jl"))

const WR_RT = MM_RT.RelaxationTime
const WR_TC = MM_RT.TransportCoefficients
const WR_WORKFLOW = Models.TransportWorkflow
const WR_SPECIES = (:u, :d, :s, :ubar, :dbar, :sbar)
const WR_KEYS = Set(WR_RT.REQUIRED_RATE_KEYS_FOR_TAU)
const WR_GRID_FILES = ("mode_a_grid.csv", "mode_b_grid.csv",
    "mechanism_transport.csv", "midwindow_joint_transport.csv")

wr_read(path) = JSON3.read(read(path, String))
wr_float(row, f) = parse(Float64, getproperty(row, f))
wr_species(row, prefix) = NamedTuple{WR_SPECIES}(Tuple(
    wr_float(row, Symbol(prefix, s)) for s in WR_SPECIES))
wr_rel(calc, ref) = abs(calc/ref - 1)

function wr_verify(path, digest)
    mm_digest(path) == digest || error("source hash mismatch: $path")
    path
end

function wr_validate_rates(rates)
    Set(keys(rates)) == WR_KEYS || error("all 17 independent channel rates are required")
    all(x -> isfinite(x) && x >= 0, values(rates)) || error("invalid cached rate")
    rates
end

function wr_read_rates(path)
    points = Dict{String,Dict{Symbol,Float64}}()
    for row in CSV.File(path)
        rates = get!(() -> Dict{Symbol,Float64}(), points, String(row.point_id))
        channel = Symbol(row.channel)
        haskey(rates, channel) && error("duplicate rate $(row.point_id)/$channel")
        rates[channel] = Float64(row.rate)
    end
    foreach(wr_validate_rates, values(points))
    points
end

function wr_assemble(densities, rates, old_tau, id; rtol)
    inverse = WR_RT.relaxation_rates(densities, rates)
    # Check all species against the archived output and isolate the one changed term.
    shift = 2densities.ubar*(rates[:subar_to_subar] - rates[:usbar_to_usbar])
    legacy = merge(inverse, (s=inverse.s-shift,))
    defect = maximum(wr_rel(legacy[s], inv(old_tau[s])) for s in WR_SPECIES)
    defect <= rtol || error("legacy rate closure at $id: $defect")
    all(x -> isfinite(x) && x > 0, inverse) || error("invalid assembled rate at $id")
    return inverse, map(inv, inverse), defect
end

function wr_grid(file, points, config; check_same_ids=true)
    rows, checks = NamedTuple[], NamedTuple[]
    for row in CSV.File(file; types=String)
        id = String(row.point_id)
        rates = points[id]
        density, old_tau = wr_species(row, "n_"), wr_species(row, "tau_")
        inv_tau, tau, rate_defect = wr_assemble(density, rates, old_tau, id; rtol=1e-11)
        masses = (u=wr_float(row,:m_u), d=wr_float(row,:m_d), s=wr_float(row,:m_s))
        mu = wr_float(row,:muq_fm)
        qp = (m=masses, μ=(u=mu,d=mu,s=mu))
        tp = (T=wr_float(row,:T_fm), Φ=wr_float(row,:Phi),
            Φbar=wr_float(row,:Phibar), ξ=wr_float(row,:xi))
        s = wr_float(row,:s_fm3inv)
        J, Z = wr_species(row,"J_"), wr_species(row,"Z_")
        eta_old = WR_TC.shear_viscosity(qp,tp; tau=old_tau,config)
        sigma_T_old = sum(old_tau[f]*J[f] for f in WR_SPECIES)
        zeta_old = sum(old_tau[f]*Z[f] for f in WR_SPECIES)
        kernel_defect = maximum((wr_rel(eta_old,wr_float(row,:eta)),
            wr_rel(sigma_T_old,wr_float(row,:sigma_over_T)),wr_rel(zeta_old,wr_float(row,:zeta))))
        kernel_defect < 1e-10 || error("legacy terminal kernel closure at $id: $kernel_defect")
        # Other species retain the archived times so the non-strange columns are exact.
        tau_used = merge(old_tau, (s=tau.s,))
        eta = WR_TC.shear_viscosity(qp,tp; tau=tau_used,config)
        zeta = sum(tau_used[f]*Z[f] for f in WR_SPECIES)
        sigma_T = sum(tau_used[f]*J[f] for f in WR_SPECIES)
        isapprox(sigma_T, WR_TC.electric_conductivity(qp,tp; tau=tau_used,config)/tp.T;
            rtol=1e-11) || error("corrected conductivity at $id")
        eta_s, zeta_s = eta/s, zeta/s
        corrected = (tau_s=tau.s, tauinv_s=inv_tau.s, eta=eta, zeta=zeta, sigma=sigma_T*tp.T,
            eta_over_s=eta_s, zeta_over_s=zeta_s, sigma_over_T=sigma_T,
            sigma_over_T_over_eta_over_s=sigma_T/eta_s, zeta_over_s_over_eta_over_s=zeta_s/eta_s,
            su_sd_fraction=2density.u*rates[:us_to_us]/inv_tau.s)
        push!(rows, merge(NamedTuple(row), corrected))
        push!(checks, (; rate_defect,kernel_defect))
    end
    ids = Set(r.point_id for r in rows)
    length(ids) == length(rows) || error("duplicate grid point ID")
    check_same_ids && ids != Set(keys(points)) && error("grid and rate IDs differ")
    return rows, (points=length(rows), legacy_rate_defect=maximum(r.rate_defect for r in checks),
        legacy_kernel_defect=maximum(r.kernel_defect for r in checks))
end

function wr_window(runtime, scan, manifest)
    centers = Dict(t.meson => wr_read(joinpath(scan,"center_$(t.meson).json")).location for t in MM_TARGETS)
    b = runtime.profile.background
    rows = NamedTuple[]
    nodes = map(collect(pairs(manifest.transport_artifact_sha256))) do (relative, digest)
        path = wr_verify(joinpath(scan, String(relative)), digest)
        (; path, digest, node=wr_read(path))
    end
    # Continue outward from each published center. This also keeps platform
    # differences in the Newton iteration from changing the selected branch.
    sort!(nodes; by=r -> abs(r.node.delta_xi))
    seeds = Dict((m, side) => collect(c.center_state) for (m,c) in centers for side in (-1,1))
    for (; path, digest, node) in nodes
        node.metadata.core_commit == MM_SHA || error("unexpected kernel at $path")
        target = mm_target(Symbol(node.meson))
        center = centers[target.meson]
        Float64(center.xi) == node.xi_mott || error("incorrect threshold center at $path")
        branch = (target.meson, node.delta_xi < 0 ? -1 : 1)
        bg = mm_equilibrium(runtime, target, node.xi; seed=seeds[branch])
        seeds[branch] = collect(bg.eq.x_state)
        inputs = WR_WORKFLOW._transport_inputs_from_equilibrium(runtime.workflow_cache,
            bg.eq, bg.T, bg.mu; xi=bg.xi, p_num=b.p_nodes, t_num=b.cos_nodes)
        therm = WR_WORKFLOW._thermo_background_from_equilibrium(runtime.workflow_cache,
            bg.eq,bg.T; xi=bg.xi,p_num=b.p_nodes,t_num=b.cos_nodes)
        rates = wr_validate_rates(Dict{Symbol,Float64}(pairs(node.rates)))
        old_tau = NamedTuple{WR_SPECIES}(Tuple(Float64(getproperty(node.tau,s)) for s in WR_SPECIES))
        _, tau, rate_defect = wr_assemble(inputs.densities, rates, old_tau, path; rtol=1e-9)
        shear(t) = WR_TC.shear_viscosity(bg.state.quark_params,bg.state.thermo_params;
            tau=t,config=runtime.transport_config)/therm.entropy
        old_eta = shear(old_tau)
        kernel_defect = wr_rel(old_eta, node.eta_over_s)
        kernel_defect < 1e-9 || error("legacy window shear closure at $path: $kernel_defect")
        # Add only the change from the corrected species weight to the archived value.
        unit_s = NamedTuple{WR_SPECIES}(Tuple(s === :s ? 1.0 : 0.0 for s in WR_SPECIES))
        eta = node.eta_over_s + (tau.s-old_tau.s)*shear(unit_s)
        push!(rows, (meson=String(node.meson),tier=String(node.tier),index=Int(node.index),
            old_eta_over_s=node.eta_over_s,eta_over_s=eta,old_tau_s=old_tau.s,tau_s=tau.s,
            n_ubar_over_u=inputs.densities.ubar/inputs.densities.u,
            rate_defect=rate_defect,kernel_defect=kernel_defect,original_node_sha256=String(digest)))
    end
    length(rows) == 132 || error("incomplete local-window input")
    sort!(rows; by=r -> (r.meson,r.tier,r.index))
end

function wr_main(args)
    length(args) == 5 || error("expected SOURCE GRID FIG7 SCAN OUT")
    _, grid, fig7, scan, out = abspath.(args)
    ispath(out) && error("use a new output directory")
    BLAS.set_num_threads(1)
    grid_manifest = wr_read(joinpath(grid,"numerical_manifest.json"))
    fig_manifest = wr_read(joinpath(fig7,"numerical_manifest.json"))
    fig_manifest.core_commit == MM_SHA || error("unexpected window kernel")
    input(name) = wr_verify(joinpath(grid,name),getproperty(grid_manifest.files,Symbol(name)))
    profile = mm_profile()
    std = MM_SP.build_runtime(profile)
    joint = MM_SP.build_runtime(MM_SP.load_profile(input("midwindow_joint_profile.toml")))
    all_rates = wr_read_rates(input("channel_rates.csv"))
    @assert length(all_rates) == 1867 && sum(length,values(all_rates)) == 31739
    mkpath(out)
    checks = Dict{String,Any}()
    for file in WR_GRID_FILES
        lookup = file == "mechanism_transport.csv" ? wr_read_rates(input("mechanism_rates.csv")) :
            file == "midwindow_joint_transport.csv" ? wr_read_rates(input("midwindow_joint_rates.csv")) : all_rates
        config = (file == "midwindow_joint_transport.csv" ? joint : std).transport_config
        corrected, check = wr_grid(input(file), lookup, config; check_same_ids=lookup !== all_rates)
        CSV.write(joinpath(out,file), corrected)
        checks[file] = check
        println("REPLAY ",file," ",check)
    end
    windows = wr_window(std, scan, fig_manifest)
    CSV.write(joinpath(out,"window_replay.csv"),windows)
    println("REPLAY Fig. 7 points=",length(windows))
    files = Dict(basename(p) => mm_digest(p) for p in readdir(out; join=true))
    mm_json(joinpath(out,"rate_assembly_replay.json"),(
        host=get(ENV,"COMPUTERNAME",Sys.MACHINE),julia=string(VERSION),
        assembly_source_sha256=mm_digest(joinpath(MM_SOURCE,"src","relaxtime","RelaxationTime.jl")),
        driver_sha256=mm_digest(@__FILE__), channel_rates_sha256=mm_digest(input("channel_rates.csv")),
        original_grid_manifest_sha256=mm_digest(joinpath(grid,"numerical_manifest.json")),
        original_window_manifest_sha256=mm_digest(joinpath(fig7,"numerical_manifest.json")),
        profile_sha256=profile.source_sha256,joint_profile_sha256=joint.profile.source_sha256,
        grid_channels=31739,window_points=length(windows),
        max_window_legacy_rate_defect=maximum(r.rate_defect for r in windows),
        max_window_legacy_kernel_defect=maximum(r.kernel_defect for r in windows),
        checks=checks,files=files))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && wr_main(ARGS)
