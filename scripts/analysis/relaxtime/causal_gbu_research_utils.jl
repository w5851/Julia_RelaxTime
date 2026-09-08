"""Analysis-only GBU contracts, independent counting and frozen-background evaluation."""
module CausalGBUResearch
using TOML, CSV, JSON3, SHA
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const METHOD = joinpath(ROOT, "config", "models", "pnjl", "charged_gbu_research_v1.toml")
isdefined(Main, :RelaxTime) || Base.include(Main, joinpath(ROOT, "src", "relaxtime", "RelaxTime.jl"))
using Main.RelaxTime.CausalSpectralBubble
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec
using Main.RelaxTime.BUPhaseGates: certify_gap_roots
using Main.RelaxTime.ChargedPhaseBackend: bu_phase_integral_parts
using Main.RelaxTime.GaussLegendre: gauleg
const CHANNELS = (:pi_plus, :pi_minus, :K_plus, :K_minus)
hashfile(p) = bytes2hex(sha256(read(p)))

function method_contract(path=METHOD)
    return validate_method(TOML.parsefile(path))
end

function validate_method(d)
    d["schema"] == "charged_gbu_research_v1" || error("unknown research method")
    for key in ("production_authorized", "meson_feedback", "additional_pion_fugacity",
                "fold_for_density", "delete_landau", "clip_density")
        d[key] === false || error("method v1 forbids $(key)")
    end
    expected = (
        observable="fixed_quark_only_gbu_partial_yield", background="FixedMuBConservedCharges",
        phase="minus_arg_retarded_inverse", frequency="external_k0",
        internal_frequency="lambda=k0+mu1-mu2", bose_weight="g(k0;mu=0)",
        weight="delta-sin(2delta)/2", counting="independent_analytic_gap_roots_plus_signed_continua",
        regulator="two_line_vacuum_thermal", vacuum_cutoff_source="model_Lambda",
        q0_reference="lambda_invariant_timelike_extrapolation",
        charge_to_baryon_ratio=0.4, strangeness_density_fm3=0.0)
    for (key,value) in pairs(expected)
        d[String(key)] == value || error("method v1 mismatch: $(key)")
    end
    for (key,value) in pairs((phase_tolerance=0.005,root_residual_tolerance=1e-8,
            static_imag_tolerance=1e-8,relative_density_target=0.01,mott_bracket_width_MeV=0.25,
            complete_curve_requires_sparse_acceptance=true))
        d["gates"][String(key)] == value || error("method v1 gate mismatch: $(key)")
    end
    for c in CHANNELS
        spec = charged_rpa_spec(c)
        d["channels"][String(c)] == String.([spec.pair..., spec.kernel_pair]) || error("channel mismatch")
    end
    return d
end

Base.@kwdef struct Settings
    mesh::Int = 128
    np::Int = 128
    nx::Int = 64
    ne::Int = 64
    nw::Int = 1200
    nr::Int = 128
    nq::Int = 12
    qmax::Float64 = 8.0
    thermal::Float64 = 20.0
    lower::Float64 = 1e-4
    upper::Float64 = 48.0
    phase_tol::Float64 = 0.005
end

function validate(s::Settings)
    all(n -> n >= 8, (s.mesh,s.np,s.nx,s.ne,s.nw,s.nr,s.nq)) || throw(ArgumentError("nodes must be >=8"))
    all(isfinite, (s.qmax,s.thermal,s.lower,s.upper,s.phase_tol)) &&
        s.qmax > 0 && s.thermal > 0 && 0 < s.lower < s.upper && s.phase_tol > 0 ||
        throw(ArgumentError("invalid integration limits"))
    return s
end

function sparse_settings()
    d=method_contract()["numerics"]
    env_int(name,key)=parse(Int,get(ENV,"GBU_SPARSE_"*name,string(d[key])))
    env_float(name,key)=parse(Float64,get(ENV,"GBU_SPARSE_"*name,string(d[key])))
    return validate(Settings(
        mesh=env_int("MESH","segment_nodes"), np=env_int("P_NODES","momentum_nodes"),
        nx=env_int("X_NODES","angle_nodes"), ne=env_int("ENERGY_NODES","energy_nodes"),
        nw=env_int("OMEGA_NODES","omega_nodes"),
        nr=env_int("ROOT_NODES","root_nodes"), nq=env_int("Q_NODES","q_nodes"),
        qmax=env_float("QMAX","qmax_inv_fm"), thermal=env_float("THERMAL_MAX","thermal_max_inv_fm"),
        lower=env_float("OMEGA_MIN","omega_min_inv_fm"), upper=env_float("OMEGA_MAX","omega_max_inv_fm")))
end

"""Argument-principle count on a CCW rectangle. This does not unwrap a density phase.

The caller must certify analyticity and absence of poles in the rectangle.
Two contour resolutions must agree, with no near-zero boundary or large step.
"""
function contour_count(f, left, right, bottom, top; nodes=128, boundary_tol=1e-9,max_nodes=2048)
    all(isfinite, (left,right,bottom,top)) && left < right && bottom < top && 8<=nodes<=max_nodes ||
        throw(ArgumentError("invalid counting rectangle"))
    function one(n)
        corners = ComplexF64[left+bottom*im,right+bottom*im,right+top*im,left+top*im]
        zs = ComplexF64[]
        for j in 1:4
            a,b = corners[j],corners[mod1(j+1,4)]
            # Cluster near the real-axis crossing without shrinking the counting domain.
            knots = real(a)==real(b) && imag(a)*imag(b)<0 ? [a,complex(real(a),0.0),b] : [a,b]
            for k in 1:length(knots)-1
                m = length(knots)==2 ? n : k==1 ? div(n,2) : n-div(n,2)
                start,stop = knots[k],knots[k+1]
                append!(zs,[start+(stop-start)*(1-cospi(t/m))/2 for t in 0:m-1])
            end
        end
        values = f.(zs)
        all(isfinite, values) || return (winding=NaN, max_step=Inf, minimum=0.0)
        increments = [angle(values[mod1(j+1,length(values))]/values[j]) for j in eachindex(values)]
        return (winding=sum(increments)/(2pi),max_step=maximum(abs,increments),minimum=minimum(abs,values))
    end
    n=nodes
    a,b = one(n),one(2n)
    function verdict(a,b)
        safe = all(isfinite, (a.winding,b.winding)) && min(a.minimum,b.minimum)>boundary_tol &&
            max(a.max_step,b.max_step)<pi/2
        count = safe ? round(Int,b.winding) : -1
        passed = safe && count >= 0 && abs(a.winding-count)<1e-7 && abs(b.winding-count)<1e-7
        return count,passed
    end
    count,passed=verdict(a,b)
    while !passed && n<max_nodes && min(a.minimum,b.minimum)>boundary_tol
        n*=2
        a,b=b,one(2n)
        count,passed=verdict(a,b)
    end
    return (count=count,passed=passed,coarse_winding=a.winding,fine_winding=b.winding,
        max_step=max(a.max_step,b.max_step),boundary_minimum=min(a.minimum,b.minimum),nodes_per_edge=2n)
end

"""Invariant q=0 reference in the internal lambda coordinate; spacelike is outside this model."""
function q0_reference_coordinate(k0, q, shift)
    q >= 0 || throw(ArgumentError("q must be nonnegative"))
    lambda = k0+shift
    lambda >= q || return nothing
    return sqrt(max(0.0,lambda^2-q^2))-shift
end

function frozen_background(input)
    hashes = Dict(f=>hashfile(joinpath(input,f)) for f in ("plot_manifest.json","charged_phase_profile_detail.csv"))
    m = JSON3.read(read(joinpath(input,"plot_manifest.json"),String))
    rows = collect(CSV.File(joinpath(input,"charged_phase_profile_detail.csv")))
    length(rows) == m.rows || error("frozen row count mismatch")
    triple(d) = (u=Float64(d.u),d=Float64(d.d),s=Float64(d.s))
    couplings = Dict(c=>Float64(first(filter(r->r.channel==String(c) && r.variant=="pv_cut",rows)).coupling_fm2) for c in CHANNELS)
    all(hashfile(joinpath(input,f))==h for (f,h) in hashes) || error("background input changed during read")
    return (m=triple(m.masses_inv_fm),mu=triple(m.chemical_potentials_inv_fm),
        T=Float64(m.thermo.T_inv_fm),Phi=Float64(m.thermo.Phi),PhiBar=Float64(m.thermo.PhiBar),
        coupling=couplings,vacuum=Float64(Main.Constants_PNJL.Λ_inv_fm),
        label="retained_T$(m.background.T_MeV)_muB$(m.background.muB_MeV)",
        T_MeV=Float64(m.background.T_MeV),muB_MeV=Float64(m.background.muB_MeV),
        residual=Float64(m.background.gap_residual_norm),solver_called=false,
        input_directory=abspath(input),input_hashes=hashes)
end

function saved_background(input, T_MeV, pn)
    path=joinpath(input,"backgrounds.csv")
    hash=hashfile(path)
    r=only(filter(r->r.T_MeV==T_MeV && r.p_nodes==pn,collect(CSV.File(path))))
    hashfile(path)==hash || error("saved background changed during read")
    return (m=(u=r.m_u,d=r.m_d,s=r.m_s),mu=(u=r.mu_u,d=r.mu_d,s=r.mu_s),
        T=r.T_MeV/Main.Constants_PNJL.ħc_MeV_fm,Phi=r.Phi,PhiBar=r.PhiBar,
        coupling=Dict(:pi_plus=>r.K12,:pi_minus=>r.K12,:K_plus=>r.K45,:K_minus=>r.K45),
        vacuum=Float64(Main.Constants_PNJL.Λ_inv_fm),T_MeV=r.T_MeV,muB_MeV=r.muB_MeV,
        residual=r.residual,label="saved_T$(T_MeV)_p$(pn)",solver_called=false,
        input_directory=abspath(input),input_hashes=Dict("backgrounds.csv"=>hash))
end

function bubble_at(bg, channel, q, s::Settings)
    validate(s)
    a,b = charged_rpa_spec(channel).pair
    g = build_spectral_bubble(q,bg.m[a],bg.mu[a],bg.m[b],bg.mu[b],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=s.thermal,momentum_nodes=s.np,angle_nodes=s.nx)
    p = build_bubble_dispersion(g;segment_nodes=s.mesh,energy_nodes=s.ne)
    shift = bg.mu[a]-bg.mu[b]
    inverse(z) = 1-4bg.coupling[channel]*cauchy_transform(p,z+shift)
    return (grid=g,profile=p,inverse=inverse,shift=shift,
        landau=hypot(q,bg.m[a]-bg.m[b])-shift,threshold=hypot(q,bg.m[a]+bg.m[b])-shift)
end

function gap_audit(b, s::Settings; count_contour=true)
    left,right = max(0.0,b.landau),b.threshold
    right-left > 2e-6 || return (roots=NamedTuple[],count=0,passed=false,status="no_normal_gap",
        contour_count=-1,contour_passed=false,root_nodes_stable=false,static_passed=false)
    gaps = [(left,right)]
    lo = certify_gap_roots((w,_)->b.inverse(w),b.grid.q_inv_fm,gaps;physical_sheet=true,real_axis=true,omega_nodes=s.nr)
    hi = certify_gap_roots((w,_)->b.inverse(w),b.grid.q_inv_fm,gaps;physical_sheet=true,real_axis=true,omega_nodes=2s.nr)
    stable = lo.passed && hi.passed && lo.count==hi.count &&
        all(abs(lo.roots[i].omega_inv_fm-hi.roots[i].omega_inv_fm)<1e-7 for i in eachindex(lo.roots))
    # Schwarz reflection extends the same real spectral Cauchy function through this gap.
    f(z) = imag(z)>=0 ? b.inverse(z) : conj(b.inverse(conj(z)))
    h = min(0.02,(right-left)/20)
    c = count_contour ? contour_count(f,left+1e-6,right-1e-6,-h,h) :
        (count=-1,passed=false)
    static = b.inverse(0.0)
    static_ok = real(static)>0 && abs(imag(static))<1e-8
    passed = stable && static_ok && (!count_contour || (c.passed && c.count==hi.count))
    return (roots=hi.roots,count=hi.count,passed=passed,status=passed ? "conditional_gap_checked" : "gap_or_static_check_failed",
        contour_count=c.count,contour_passed=c.passed,root_nodes_stable=stable,static_passed=static_ok)
end

"""One-sided phase-limit check, independent of the expected number of roots.

The last two probes must agree. An unresolved threshold zero is not rounded
into a bound-state count; the finite offsets and all sampled phases are retained.
"""
function threshold_phase_limit(inverse, threshold; offsets=(1e-7,1e-8,1e-9,1e-10),
                               phase_tol=0.005,zero_tol=1e-8)
    isfinite(threshold) && length(offsets)>=2 && all(isfinite,offsets) &&
        all(>(0),offsets) && all(<(0),diff(collect(offsets))) && phase_tol>0 && zero_tol>0 ||
        throw(ArgumentError("invalid threshold-limit probes"))
    energies=threshold .+ collect(offsets)
    all(>(threshold),energies) && all(<(0),diff(energies)) ||
        throw(ArgumentError("threshold offsets are unresolved at machine precision"))
    values=inverse.(energies)
    at_threshold=inverse(threshold)
    safe=all(isfinite,values) && isfinite(at_threshold) &&
        abs(at_threshold)>zero_tol && all(v->abs(v)>zero_tol,values)
    phases=safe ? -angle.(values) : fill(NaN,length(values))
    change=abs(last(phases)-phases[end-1])/pi
    passed=safe && maximum(abs,diff(phases))<pi && change<phase_tol/10
    return (phase=last(phases),offset=last(offsets),phase_change_over_pi=change,
        offsets=collect(offsets),phases=phases,threshold_inverse_abs=abs(at_threshold),passed=passed)
end

function shell(bg,channel,q,s::Settings; reference=nothing,count_contour=true)
    b = bubble_at(bg,channel,q,s)
    ga = gap_audit(b,s;count_contour=count_contour)
    inverse = b.inverse
    roots = [r.omega_inv_fm for r in ga.roots]
    reference_boundary_ok = true
    if reference !== nothing
        inverse = w -> begin
            v = q0_reference_coordinate(w,q,b.shift)
            v === nothing ? complex(1.0,0.0) : reference.bubble.inverse(v)
        end
        roots = [hypot(r.omega_inv_fm+b.shift,q)-b.shift for r in reference.gap.roots]
        reference_boundary_ok = abs(angle(reference.bubble.inverse(-b.shift))) < s.phase_tol
    end
    s.upper > b.threshold+0.01 || error("upper endpoint is below threshold")
    limit=threshold_phase_limit(inverse,b.threshold;phase_tol=s.phase_tol)
    # Resolve the one-sided threshold limit separately from the bulk energy mesh.
    uw = sort!(unique(vcat(b.threshold .+ [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3],
        collect(range(b.threshold+0.002,s.upper;length=s.nw)))))
    ld = b.landau
    lw = ld > s.lower+1e-8 ? exp.(range(log(s.lower),log(ld-1e-8);length=s.nw)) : Float64[]
    if reference !== nothing && !isempty(lw)
        onset = q-b.shift
        s.lower < onset < ld && (lw = sort!(unique(vcat(lw,[onset,onset+1e-8]))))
    end
    up,lp = [-angle(inverse(w)) for w in uw],[-angle(inverse(w)) for w in lw]
    branches = maximum(abs,diff(up))<pi && (length(lp)<2 || maximum(abs,diff(lp))<pi)
    u = bu_phase_integral_parts(uw,up,bg.T;weight=:gbu)
    l = isempty(lw) ? nothing : bu_phase_integral_parts(lw,lp,bg.T;weight=:gbu)
    bose_ok = all(w->w>0,roots)
    bound = bose_ok ? sum((1/expm1(w/bg.T) for w in roots);init=0.0) : NaN
    landau = l === nothing ? 0.0 : l.derivative
    active_gap = reference===nothing ? ga : reference.gap
    count_ok = active_gap.passed
    lr = limit.phase/pi-length(roots)
    gates=(gap=count_ok,bose=bose_ok,branch=branches,levinson=abs(lr)<s.phase_tol,
        threshold_limit=limit.passed,upper=abs(last(up))<s.phase_tol,reference_boundary=reference_boundary_ok)
    failed_gates=join([String(k) for (k,v) in pairs(gates) if !v],";")
    pass = isempty(failed_gates)
    return (channel=String(channel),q_inv_fm=Float64(q),route=reference===nothing ? "direct_finite_q" : "q0_lambda_reference",
        bound=bound,unitary=u.derivative,landau=landau,shell_inv_fm2=q^2/(2pi^2)*(bound+u.derivative+landau),
        gap_count=length(roots),contour_count=active_gap.contour_count,contour_passed=active_gap.contour_passed,
        static_passed=active_gap.static_passed,root_nodes_stable=active_gap.root_nodes_stable,levinson_residual=lr,
        threshold_phase_over_pi=limit.phase/pi,threshold_offset=limit.offset,
        threshold_phase_at_1e7_over_pi=first(limit.phases)/pi,
        threshold_limit_change=limit.phase_change_over_pi,threshold_limit_passed=limit.passed,
        upper_phase=last(up),branch_passed=branches,
        reference_boundary_passed=reference_boundary_ok,failed_gates=failed_gates,
        passed=pass,status=pass ? "conditional_window" : "gate_failed",
        production_authorized=false)
end

function density(bg,channel,s; reference=false,count_contour=true)
    q,weights = gauleg(0.0,s.qmax,s.nq)
    ref = if reference
        b = bubble_at(bg,channel,0.0,s)
        (bubble=b,gap=gap_audit(b,s;count_contour=count_contour))
    else
        nothing
    end
    rows = [shell(bg,channel,v,s;reference=ref,count_contour=count_contour) for v in q]
    value = sum(weights[i]*rows[i].shell_inv_fm2 for i in eachindex(q))
    return (density=value,passed=all(r.passed for r in rows) && isfinite(value) && value>=0,rows=rows)
end

"""New directory only; source snapshots bind diagnostics to loaded research code."""
function start_output(output)
    ispath(output) && error("refusing to overwrite $(output)")
    mkpath(output)
    paths = String[]
    for dir in (joinpath(ROOT,"src"),joinpath(ROOT,"config"))
        for (root,_,files) in walkdir(dir), file in files
            endswith(file,".jl") || endswith(file,".toml") || continue
            push!(paths,joinpath(root,file))
        end
    end
    append!(paths,[joinpath(@__DIR__,f) for f in readdir(@__DIR__) if startswith(f,"causal_gbu_") || startswith(f,"audit_causal_gbu_")])
    hashes = Dict(relpath(p,ROOT)=>hashfile(p) for p in paths)
    for p in paths
        target = joinpath(output,"source_snapshot",relpath(p,ROOT))
        mkpath(dirname(target)); cp(p,target)
    end
    return hashes
end

function finish_output(output, hashes, record)
    all(hashfile(joinpath(ROOT,p))==h for (p,h) in hashes) || error("source changed during run")
    record["source_hashes"] = hashes
    record["method_sha256"] = hashfile(METHOD)
    record["git_head"] = readchomp(`git -C $ROOT rev-parse HEAD`)
    record["production_authorized"] = false
    record["output_hashes"] = Dict(f=>hashfile(joinpath(output,f)) for f in readdir(output) if isfile(joinpath(output,f)))
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,record)
    end
end
end
