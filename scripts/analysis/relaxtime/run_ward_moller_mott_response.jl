#!/usr/bin/env julia
# Recompute local Mott windows on the pinned Ward/Møller production kernel.
# julia --project=SOURCE run_ward_moller_mott_response.jl SOURCE OUT locate
# julia --project=SOURCE run_ward_moller_mott_response.jl SOURCE OUT scan MESON TIER INDEX COUNT
using CSV, JSON3, SHA, LinearAlgebra

const MM_SOURCE = abspath(get(ARGS, 1, joinpath(@__DIR__, "..", "..", "..")))
if !isdefined(Main, :Models)
    include(joinpath(MM_SOURCE, "src", "models", "Models.jl"))
end
include(joinpath(MM_SOURCE, "scripts", "relaxtime", "transport_sampling_scan_config.jl"))
include(joinpath(MM_SOURCE, "scripts", "relaxtime", "gap_transport_scan_phase_equilibrium.jl"))
const MM_RT = Main.RelaxTime
const MM_AS = MM_RT.AverageScatteringRate
const MM_SP = Models.TransportSamplingProfiles
const MM_PHASE = GapTransportScanPhaseEquilibrium
const MM_MP = MM_RT.MesonMass
const MM_EC = MM_RT.EffectiveCouplings
const MM_PAIR = MM_RT.PairLoopQuadrature
const MM_HC = Main.Constants_PNJL.ħc_MeV_fm
const MM_SHA = "096199b96d34045bf8e2d9611c2d786d6f92fb8e"
const MM_PROFILE_SHA = "a26d470cf09730ba79d4a440bf9e3d7bc17c202f59323591edafe4a47d249313"

const MM_TARGETS = (
    (meson=:pi, T_MeV=138.3096002910106, muB_MeV=900.0, search=(0.17, 0.19), scout_nodes=5),
    (meson=:K, T_MeV=200.0898076157829, muB_MeV=0.0, search=(-0.20, -0.16), scout_nodes=5),
    (meson=:sigma, T_MeV=125.7360002645551, muB_MeV=900.0, search=(-0.24, 0.24), scout_nodes=25),
)
mm_target(meson) = only(t for t in MM_TARGETS if t.meson === meson)
mm_digest(p) = bytes2hex(sha256(read(p)))

function mm_json(path, record)
    open(path * ".writing", "w") do io
        JSON3.pretty(io, record)
        println(io)
    end
    mv(path * ".writing", path; force=true)
end

function mm_profile()
    path = joinpath(MM_SOURCE, "config", "models", "pnjl", "transport_sampling",
        "ward_moller_sphere_v1.toml")
    mm_digest(path) == MM_PROFILE_SHA || error("Unexpected production profile hash")
    profile = MM_SP.load_profile(path)
    @assert profile.physics == (a_method=:ward_matched, b0_direction=:sphere_average)
    @assert profile.sampling.relative_velocity === :medium_moller
    @assert (profile.background.p_nodes, profile.background.cos_nodes,
        profile.background.thermal_p_max_inv_fm) == (96, 24, 25.0)
    profile
end

function mm_runtime(profile, tier)
    runtime = MM_SP.build_runtime(profile)
    tier == "production" && return runtime
    rtol = tier == "fine" ? 2e-6 : tier == "reference" ? 2e-7 :
        throw(ArgumentError("unknown sampling tier $tier"))
    cfg = MM_AS.ThermalSamplingConfig(; merge(MM_AS._thermal_sampling_identity(profile.sampling),
        (rtol=rtol, atol=1e-13, max_evaluations=160_000, max_nodes=48_000, max_rounds=70))...)
    merge(runtime, (tau_kwargs=merge(runtime.tau_kwargs, (sampling_config=cfg,)),))
end

function mm_equilibrium(runtime, target, xi; seed=nothing)
    p = runtime.profile.background
    T, mu = target.T_MeV/MM_HC, target.muB_MeV/(3MM_HC)
    opts = (model=runtime.model, p_num=p.p_nodes, t_num=p.cos_nodes,
        residual_norm_max=1e-11, max_iter=100)
    eq = if seed === nothing
        MM_PHASE.solve_with_multiseed_governance(T, mu, xi, opts)
    else
        continuation = Float64.(collect(seed))
        raw = Models.solve_constraint(runtime.model, Models.FixedMu(), T;
            μ_fm=mu, seed_guess=continuation, seed_candidates=[continuation],
            xi, p_num=p.p_nodes, t_num=p.cos_nodes,
            solver=Models.NLsolveGapSolver(method=:newton,jacobian=:forward,xtol=1e-14,ftol=2e-14),
            residual_norm_max=1e-12, iterations=100, physicality_check=Models.is_physical_solution)
        MM_PHASE._normalize_equilibrium_result(raw;solver_backend=:models_strict_local)
    end
    eq !== nothing && eq.residual_norm <= opts.residual_norm_max ||
        error("mean-field solution failed for $(target.meson) at xi=$xi")
    (; eq, T, mu, xi, state=MM_SP.point_state(runtime.profile, eq, T, mu; xi))
end

# Exact k=0 endpoints of the same on-shell pair loops. Pair subtraction avoids
# loss of precision at a pole touching the integration endpoint. Common mu
# and equal light masses are enforced here and used by the manuscript scan.
function mm_endpoint(bg, meson; nodes=128)
    st = bg.state
    qp, tp, kc = st.quark_params, st.thermo_params, st.K_coeffs
    mu = qp.μ.u
    @assert mu == qp.μ.s && isapprox(qp.m.u, qp.m.d; rtol=1e-12)
    m, ms = qp.m.u, qp.m.s
    pu = MM_PAIR.PairLoopPlan(m, mu, tp; n_nodes=nodes)
    N = Main.Constants_PNJL.N_color/(8pi^2)
    if meson === :pi
        pol = -N*(qp.A.u+qp.A.d-8m^2*pu.threshold_integral)
        return complex(1-4MM_MP._meson_coupling(:pi, kc)*pol)
    end
    ps = MM_PAIR.PairLoopPlan(ms, mu, tp; n_nodes=nodes)
    if meson === :K
        b0 = 2/(m+ms)*(m*pu.threshold_integral+ms*ps.threshold_integral)
        pol = -N*(qp.A.u+qp.A.s-4m*ms*b0)
        return complex(1-4MM_MP._meson_coupling(:K, kc)*pol)
    elseif meson === :sigma
        polu = complex(-2N*qp.A.u)
        pols = -N*(2qp.A.s + 4(ms^2-m^2)*MM_PAIR.pair_b0_pv(ps, 4m^2))
        return MM_MP._mixed_inverse(:S, :light, MM_EC.mixing_matrix_elements(polu, pols, kc, :S))
    end
    throw(ArgumentError("unknown meson $meson"))
end

function mm_native_inverse(bg, meson, x)
    st = bg.state
    qp = st.quark_params
    w = meson === :K ? qp.m.u+qp.m.s : 2qp.m.u
    MM_MP.meson_mass_equation(meson, sqrt(w^2+x), 0.0, 0.0, qp, st.thermo_params, st.K_coeffs)
end

function mm_bisect(f, lo, hi; atol=2e-13)
    fl, fh = f(lo), f(hi)
    all(isfinite, (fl, fh)) && signbit(fl) != signbit(fh) || error("invalid threshold bracket")
    while hi-lo > atol
        x = (lo+hi)/2
        fx = f(x)
        isfinite(fx) || error("nonfinite threshold inverse")
        if signbit(fx) == signbit(fl)
            lo, fl = x, fx
        else
            hi, fh = x, fx
        end
    end
    (xi=(lo+hi)/2, bracket=(lo, hi), endpoint_values=(fl, fh))
end

function mm_location(runtime, target, outdir)
    function endpoint(b, n)
        v = mm_endpoint(b, target.meson; nodes=n)
        abs(imag(v)) < 1e-9 || error("open cut at the candidate threshold")
        real(v)
    end
    xs = collect(range(target.search...; length=target.scout_nodes))
    coarse_bgs = [mm_equilibrium(runtime, target, x) for x in xs]
    scout = [(xi=b.xi, inverse=endpoint(b, 128), light_MeV=MM_HC*b.eq.masses[1]) for b in coarse_bgs]
    CSV.write(joinpath(outdir, "$(target.meson)_scout.csv"), scout)
    brackets = [(l.xi, r.xi) for (l, r) in zip(scout[1:end-1], scout[2:end])
        if signbit(l.inverse) != signbit(r.inverse) && abs(l.light_MeV-r.light_MeV) < 20.0]
    length(brackets) == 1 || error("Expected one smooth $(target.meson) Mott threshold; found $brackets")
    lo, hi = only(brackets)
    seed = only(b.eq.x_state for b in coarse_bgs if b.xi==lo)
    bgs = Dict{Float64,NamedTuple}()
    bg(x) = get!(bgs, x) do
        mm_equilibrium(runtime, target, x;seed)
    end
    coarse = mm_bisect(x->endpoint(bg(x), 128), lo, hi)
    fine = mm_bisect(x->endpoint(bg(x), 256), lo, hi)
    abs(coarse.xi-fine.xi)<1e-8 || error("pair quadrature did not resolve $(target.meson) threshold")
    center = bg(fine.xi)
    native = [let v=mm_native_inverse(center, target.meson, x)
        (x_inv_fm2=x, re=real(v), im=imag(v))
    end for x in (-1e-4, -1e-6, 1e-6, 1e-4)]
    guards = [bg(fine.xi+d) for d in (-0.008, -0.004, 0.0, 0.004, 0.008)]
    for b in guards
        independent = mm_equilibrium(runtime, target, b.xi)
        maximum(abs, b.eq.masses-independent.eq.masses)<1e-6 ||
            error("continued branch is not the minimum of Omega")
    end
    masses = [b.eq.masses[1]*MM_HC for b in guards]
    maximum(abs, diff(masses)) < 20.0 || error("discontinuous background in threshold window")
    println("MOTT ", target.meson, " xi=", fine.xi, " M_u [MeV]=", masses)
    (target..., xi=fine.xi, bracket=fine.bracket, endpoint_values=fine.endpoint_values,
        root_128_shift=coarse.xi-fine.xi, center_inverse=endpoint(center, 256),
        native_offsets=native, branch_masses_MeV=masses,
        center_masses_inv_fm=center.eq.masses, center_state=center.eq.x_state)
end

mm_offsets() = sort!(unique(vcat(collect(-32:4:32), collect(-8:8),
    [-1.5,-0.5,0.5,1.5])) .* 0.00025)
mm_reference_offsets() = 0.00025 .* [-32,-28,-4,-1,-0.5,0,0.5,1,4,28,32]

function mm_metadata(runtime)
    (core_commit=MM_SHA, run_host=get(ENV, "COMPUTERNAME", get(ENV, "HOSTNAME", Sys.MACHINE)),
        julia_version=string(VERSION),
        driver_sha256=mm_digest(@__FILE__), profile=MM_SP.effective_config(runtime.profile),
        runtime_sampling=MM_AS._thermal_sampling_identity(runtime.tau_kwargs.sampling_config))
end

function mm_scan(runtime, target, location, rowdir, offsets, index, count, tier)
    1 <= index <= count || error("invalid shard index")
    mkpath(rowdir)
    for i in index:count:length(offsets)
        path = joinpath(rowdir, lpad(i, 3, '0')*".json")
        isfile(path) && error("refusing to overwrite completed node $path")
        u = offsets[i]
        bg = mm_equilibrium(runtime, target, Float64(location.xi)+u; seed=location.center_state)
        started = time()
        result = MM_SP.solve_profile_transport(runtime, bg.eq, bg.T, bg.mu;
            xi=bg.xi, compute_bulk=false)
        eta_s = result.transport.eta/result.thermo_background.entropy
        row = (; metadata=mm_metadata(runtime), meson=target.meson, target.T_MeV, target.muB_MeV,
            tier=tier, index=i, xi=bg.xi, xi_mott=Float64(location.xi), delta_xi=u, eta_over_s=eta_s,
            m_u_MeV=bg.eq.masses[1]*MM_HC, m_s_MeV=bg.eq.masses[3]*MM_HC,
            omega=bg.eq.omega, residual=bg.eq.residual_norm,
            tau=result.tau, rates=result.rates, rate_diagnostics=result.sampling_diagnostics,
            seconds=time()-started)
        mm_json(path, row)
        println("NODE $tier $(target.meson) $i/$(length(offsets)) u=$u eta/s=$eta_s seconds=$(round(row.seconds; digits=1))")
        flush(stdout)
    end
end

function mm_main(args)
    length(args) in (3, 4, 7) || error("SOURCE OUT locate [MESON] | SOURCE OUT scan MESON TIER INDEX COUNT")
    outdir, stage = abspath(args[2]), args[3]
    mkpath(outdir)
    LinearAlgebra.BLAS.set_num_threads(1)
    profile = mm_profile()
    if stage == "locate" && length(args) in (3,4)
        runtime = mm_runtime(profile, "production")
        targets = length(args)==3 ? MM_TARGETS : (mm_target(Symbol(args[4])),)
        for t in targets
            path = joinpath(outdir,"center_$(t.meson).json")
            isfile(path) && error("center already exists: $path")
            mm_json(path,(metadata=mm_metadata(runtime),location=mm_location(runtime,t,outdir)))
        end
    elseif stage == "scan" && length(args) == 7
        meson, tier = Symbol(args[4]), args[5]
        target = mm_target(meson)
        record = JSON3.read(read(joinpath(outdir, "center_$meson.json"), String))
        record.metadata.core_commit == MM_SHA || error("unexpected threshold provenance")
        location = record.location
        Symbol(location.meson) === meson || error("unexpected threshold selection")
        index, count = parse.(Int, args[6:7])
        for precision in (tier=="both" ? ("fine","reference") : (tier,))
            runtime = mm_runtime(profile, precision)
            offsets = precision == "reference" ? mm_reference_offsets() : mm_offsets()
            mm_scan(runtime, target, location, joinpath(outdir, precision, String(meson)),
                offsets, index, count, precision)
        end
    else
        error("unknown or malformed stage")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    mm_main(ARGS)
end
