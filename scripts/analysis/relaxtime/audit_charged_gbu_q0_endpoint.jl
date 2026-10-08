#!/usr/bin/env julia

isdefined(Main, :ChargedGBUContourScan) || include("run_charged_gbu_contour_scan.jl")

"""Small, analysis-only endpoint audit on exactly three retained backgrounds.

No equilibrium solve, scan, gate relaxation or density promotion is performed.
The old finite-Lth spectrum and the infinite target are compared at the SAME
coordinates. Finite IR windows are probes, never repaired production yields.
"""
module ChargedGBUQ0EndpointAudit

using CSV, JSON3, SHA
using Main.Models

const Scan = Main.ChargedGBUContourScan
const W = Scan.Workflow
const I = W.I
const P = W.P
const R = W.R
const Ref = Scan.Reference
const ROOT = Scan.ROOT
const HISTORY = joinpath(ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_phase_backend", "fig4_like_freezeout_ratio_dense_20260905_v3")
const SOURCE_RUN = "36971423991"
const SOURCE_SHA = "e6b3de025d59213e480927ebffa5f188e3e71ebd"
const LEVELS = ((label="screening", mesh=64, cut=32, tail=32),
    (label="refined128", mesh=128, cut=64, tail=64),
    (label="refined256", mesh=256, cut=96, tail=96))
const IMAG_TOL = 1e-8 # Identical to the retained shell gate, not a new tolerance.

hashfile(p) = bytes2hex(sha256(read(p)))
readjson(p) = JSON3.read(read(p, String))
writecsv(path,rows) = CSV.write(path,rows;transform=(col,value)->value===nothing ? missing : value)

"""Separate the two original static-gate conditions; do not diagnose a phase."""
function inverse_record(f)
    finite = isfinite(f)
    return (real_inverse=real(f), imaginary_inverse=imag(f), phase=-angle(f),
        gbu_weight=W.Y.O.gbu_weight(-angle(f)), finite=finite,
        real_positive=finite && real(f)>0,
        imaginary_passed=finite && abs(imag(f))<IMAG_TOL,
        static_passed=finite && real(f)>0 && abs(imag(f))<IMAG_TOL)
end

function restore_seed(model, seed, T_MeV, muB_MeV, residual)
    length(seed)==8 && all(isfinite, seed) || throw(ArgumentError("finite eight-component saved seed required"))
    T_MeV>0 && isfinite(T_MeV) && isfinite(muB_MeV) && isfinite(residual) ||
        throw(ArgumentError("finite saved background coordinates and residual required"))
    state = Models.meanfield_state(Float64.(seed[1:5]))
    masses = Models.calculate_mass_vec(model, state.phi)
    interaction = Main.RelaxTime.MesonInteractionKernel.build_full_kmt_interaction(
        state.phi; G=model.params.G_fm2, K=model.params.K_fm5)
    coupling = Dict(ch=>Main.RelaxTime.ChargedRPAKernel.charged_rpa_coupling(
        R.charged_rpa_spec(ch), interaction) for ch in R.CHANNELS)
    tri(x) = (u=Float64(x[1]), d=Float64(x[2]), s=Float64(x[3]))
    return (m=tri(masses), mu=tri(seed[6:8]), T=T_MeV/Main.Constants_PNJL.ħc_MeV_fm,
        Phi=Float64(state.Phi), PhiBar=Float64(state.PhiBar), coupling=coupling,
        vacuum=Float64(model.params.Λ_inv_fm), T_MeV=Float64(T_MeV),
        muB_MeV=Float64(muB_MeV), residual=Float64(residual))
end

function retained_cases(input_root; history=HISTORY)
    model = Models.create_model(:PNJL) # Parameter restoration only; never solve.
    Float64(model.params.Λ_inv_fm)==Float64(Main.Constants_PNJL.Λ_inv_fm) ||
        error("historical and current vacuum cutoffs differ")
    manifest_path = joinpath(history, "manifest.json")
    old = readjson(manifest_path)
    old.schema=="charged_gbu_dense_freezeout_v3" || error("unexpected historical schema")
    inputs = Dict{String,String}(manifest_path=>hashfile(manifest_path))
    for (rel, expected) in pairs(old.output_hashes)
        path=joinpath(history, String(rel))
        hashfile(path)==expected || error("historical output hash mismatch: $(rel)")
        inputs[path]=hashfile(path)
    end
    rows = collect(CSV.File(joinpath(history, "backgrounds.csv")))
    cases = NamedTuple[]
    for (label,T,muB) in (("historical_7p7GeV",139.61211954370506,421.6498501015441),
            ("historical_3GeV",79.95694382094246,719.0764156129742))
        r = only(filter(r->isapprox(r.T_MeV,T;atol=1e-9,rtol=0) &&
            isapprox(r.muB_MeV,muB;atol=1e-9,rtol=0), rows))
        bg=(m=(u=r.m_u,d=r.m_d,s=r.m_s),mu=(u=r.mu_u,d=r.mu_d,s=r.mu_s),
            T=r.T_MeV/Main.Constants_PNJL.ħc_MeV_fm,Phi=r.Phi,PhiBar=r.PhiBar,
            coupling=Dict(:pi_plus=>r.K12,:pi_minus=>r.K12,:K_plus=>r.K45,:K_minus=>r.K45),
            vacuum=Float64(model.params.Λ_inv_fm),T_MeV=r.T_MeV,muB_MeV=r.muB_MeV,residual=r.residual)
        push!(cases,(label=label,bg=bg,seed=nothing,
            condensates=[r.phi_u,r.phi_d,r.phi_s],source="retained historical backgrounds.csv"))
    end
    matches = [joinpath(dir,f) for (dir,_,files) in walkdir(input_root) for f in files
        if f=="point_T22_mu16_145p000000_375p000000.json"]
    path = only(matches); data=readjson(path)
    path_manifest=joinpath(dirname(dirname(path)), "manifest.json"); m=readjson(path_manifest)
    data.schema=="charged_gbu_contour_point_v2" && m.schema=="charged_gbu_contour_scan_v2" &&
        data.T_MeV==145 && data.muB_MeV==375 && data.scan_identity==m.scan_identity &&
        m.git_head==SOURCE_SHA && m.density_route=="q0_lambda_reference" || error("frozen point provenance mismatch")
    for (rel,expected) in pairs(m.source_hashes)
        hashfile(joinpath(ROOT,String(rel)))==expected || error("numerical source changed: $(rel)")
    end
    inputs[path]=hashfile(path);inputs[path_manifest]=hashfile(path_manifest)
    seed=Float64.(data.background.seed)
    bg=restore_seed(model,seed,data.T_MeV,data.muB_MeV,Float64(data.background.residual))
    push!(cases,(label="grid_145_375_onset",bg=bg,seed=seed,condensates=seed[1:3],
        source="frozen Actions run $(SOURCE_RUN); algebraic restoration only"))
    return cases,inputs,old.settings
end

function q_probes(shift)
    q8,_=R.gauleg(0.,8.,8);q24,_=R.gauleg(0.,8.,24)
    extras=shift>0 ? [shift/2,shift*(1-1e-4),shift*(1+1e-4)] : Float64[]
    return sort!(unique!(vcat([0.,first(q24)],q8,extras)))
end

function coordinate_rows(bg,ch)
    i,j=R.charged_rpa_spec(ch).pair;shift=bg.mu[i]-bg.mu[j]
    landau=abs(bg.m[i]-bg.m[j]);threshold=bg.m[i]+bg.m[j]
    return [(q_inv_fm=q,lambda_inv_fm=Ref.reference_coordinate(0.,q,shift),
        in_timelike_map=Ref.reference_coordinate(0.,q,shift)!==nothing,
        in_q0_landau=begin x=Ref.reference_coordinate(0.,q,shift);x!==nothing && 0<abs(x)<landau end,
        in_q0_analytic_gap=begin x=Ref.reference_coordinate(0.,q,shift);x!==nothing && landau<abs(x)<threshold end)
        for q in q_probes(shift)]
end

"""Finite-window derivative/bulk/endpoint identity, NOT a zero-endpoint yield."""
function window_record(phase,T,q,lo,hi;nodes=4800)
    0<lo<hi && T>0 && q>=0 && nodes>=8 || throw(ArgumentError("positive finite-window bounds required"))
    xs=exp.(range(log(lo),log(hi);length=nodes));xs[1]=lo;xs[end]=hi
    parts=R.bu_phase_integral_parts(xs,phase.(xs),T;weight=:gbu)
    measure=q^2/(2pi^2)
    return (lower_inv_fm=lo,upper_inv_fm=hi,nodes=nodes,
        derivative_inv_fm2=measure*parts.derivative,bulk_inv_fm2=measure*parts.bulk,
        lower_boundary_inv_fm2=measure*parts.lower_boundary,
        upper_boundary_inv_fm2=measure*parts.upper_boundary,
        identity_error_inv_fm2=measure*(parts.derivative-parts.reconstructed))
end

function endpoint_rows(case,ch,p,level;raw=false)
    bg=case.bg;k=p.kernel;K=bg.coupling[ch]
    entries=NamedTuple[(point="q0_external_static",q_inv_fm=0.,lambda_inv_fm=k.shift),
        (point="q0_internal_origin",q_inv_fm=0.,lambda_inv_fm=0.)]
    append!(entries,[(point="boosted_external_static",q_inv_fm=c.q_inv_fm,
        lambda_inv_fm=c.lambda_inv_fm) for c in coordinate_rows(bg,ch)])
    rows=NamedTuple[]
    for e in entries
        x=e.lambda_inv_fm
        f=x===nothing ? complex(1.) : P.inverse(p,x)
        rho=x===nothing ? 0. : I.total_cut(bg,ch,0.,x;nodes=96)
        raw_imag=-4K*rho
        raw_values=raw && x!==nothing && e.q_inv_fm in (0.,first(R.gauleg(0.,8.,8)[1])) ?
            [(n,1-4K*I.polarization(k,complex(x);nodes=n)) for n in (128,256)] : []
        push!(rows,merge((case=case.label,channel=String(ch),level=level),e,inverse_record(f),
            (raw_imaginary_inverse=raw_imag,profile_imag_error=imag(f)-raw_imag,
             raw_real128=isempty(raw_values) ? nothing : real(raw_values[1][2]),
             raw_real256=isempty(raw_values) ? nothing : real(raw_values[2][2]),
             raw_imag128=isempty(raw_values) ? nothing : imag(raw_values[1][2]),
             raw_imag256=isempty(raw_values) ? nothing : imag(raw_values[2][2]),
             shift_inv_fm=k.shift,landau_edge_inv_fm=k.landau,
             timelike_map=x!==nothing,in_landau=x!==nothing && 0<abs(x)<k.landau)))
    end
    return rows
end

function restore_settings(settings)
    return R.Settings(; (Symbol(key)=>Float64(val) for (key,val) in pairs(settings)
        if String(key) in ("qmax","thermal","lower","upper","phase_tol"))...,
        mesh=Int(settings.mesh),np=Int(settings.np),nx=Int(settings.nx),ne=Int(settings.ne),
        nw=Int(settings.nw),nr=Int(settings.nr),nq=Int(settings.nq))
end

function old_rows(case,ch,settings)
    bg=case.bg;s=restore_settings(settings)
    b=R.bubble_at(bg,ch,0.,s)
    rows=NamedTuple[]
    entries=NamedTuple[(point="q0_external_static",q_inv_fm=0.,lambda_inv_fm=b.shift),
        (point="q0_internal_origin",q_inv_fm=0.,lambda_inv_fm=0.)]
    append!(entries,[(point="boosted_external_static",q_inv_fm=c.q_inv_fm,
        lambda_inv_fm=c.lambda_inv_fm) for c in coordinate_rows(bg,ch)])
    for e in entries
        f=e.lambda_inv_fm===nothing ? complex(1.) : b.inverse(e.lambda_inv_fm-b.shift)
        push!(rows,merge((case=case.label,channel=String(ch),level="old_Lth24_mesh512"),e,
            inverse_record(f),(shift_inv_fm=b.shift,thermal_max_inv_fm=s.thermal,)))
    end
    return rows,b
end

function infrared_rows(case,ch,p,old)
    bg=case.bg;k=p.kernel
    rows=NamedTuple[]
    for q in (first(R.gauleg(0.,8.,8)[1]),first(R.gauleg(0.,8.,24)[1]))
        k.shift>q || continue
        endpoint=Ref.reference_coordinate(0.,q,k.shift)
        0<endpoint<k.landau || continue
        hi=min(.02,(hypot(k.landau,q)-k.shift)/4)
        hi>1e-4 || continue
        for (name,f) in (("infinite_refined256",x->P.inverse(p,x)),
                        ("old_Lth24_mesh512",x->old.inverse(x-k.shift)))
            phase=w->-angle(f(Ref.reference_coordinate(w,q,k.shift)))
            for lo in (1e-3,1e-4,1e-5,1e-6,1e-7)
                lo<hi || continue
                push!(rows,merge((case=case.label,channel=String(ch),kernel=name,
                    q_inv_fm=q,phase_at_zero=phase(0.),weight_at_zero=W.Y.O.gbu_weight(phase(0.))),
                    window_record(phase,bg.T,q,lo,hi)))
            end
        end
    end
    return rows
end

function direct_rows(case,ch,level)
    bg=case.bg;q=first(R.gauleg(0.,8.,8)[1])
    k=I.kernel(bg,ch,q;cut_nodes=level.cut,split_inv_fm=max(36.,q+28bg.T+2.))
    p=P.profile(k;mesh=level.mesh,tail_nodes=level.tail)
    f=P.inverse(p,k.shift);raw_imag=-4bg.coupling[ch]*I.total_cut(bg,ch,q,k.shift;nodes=96)
    raw=level==last(LEVELS) ? [1-4bg.coupling[ch]*I.polarization(k,complex(k.shift);nodes=n)
        for n in (128,256)] : []
    return merge((case=case.label,channel=String(ch),level=level.label,
        q_inv_fm=q,lambda_inv_fm=k.shift,raw_imaginary_inverse=raw_imag,
        profile_imag_error=imag(f)-raw_imag,
        raw_real128=isempty(raw) ? nothing : real(raw[1]),raw_real256=isempty(raw) ? nothing : real(raw[2])),inverse_record(f))
end

function parse_args(args)
    "--help" in args && return nothing
    length(args)==4 || throw(ArgumentError("--input-root <frozen-shards> --output <new-directory> required"))
    pairs=Dict(args[j]=>args[j+1] for j in 1:2:length(args))
    Set(keys(pairs))==Set(["--input-root","--output"]) || throw(ArgumentError("unsupported audit arguments"))
    return (input_root=abspath(pairs["--input-root"]),output=abspath(pairs["--output"]))
end

function run_audit(options)
    ispath(options.output) && throw(ArgumentError("audit output exists; refusing overwrite"))
    cases,inputs,oldsettings=retained_cases(options.input_root)
    length(cases)==3 || error("exactly three approved retained backgrounds required")
    mkpath(options.output)
    hashes=W.snapshot(options.output,W.DEFAULT_CONFIG)
    W.writejson(joinpath(options.output,"inputs.json"),(source_run=SOURCE_RUN,source_sha=SOURCE_SHA,
        hashes=inputs,cases=cases,solver_called=false))
    # Retain the actual selected input bytes, not only absolute-path pointers.
    for (j,path) in enumerate(sort!(collect(keys(inputs))))
        target=joinpath(options.output,"frozen_inputs","input$(j)_"*basename(path))
        mkpath(dirname(target));cp(path,target)
    end
    endpoints=NamedTuple[];historical=NamedTuple[];direct=NamedTuple[];infrared=NamedTuple[]
    for case in cases,ch in (:K_plus,:pi_plus,:K_minus)
        println("[q0-endpoint-audit] $(case.label) $(ch)");flush(stdout)
        fine=nothing
        for level in LEVELS
            k=I.kernel(case.bg,ch,0.;cut_nodes=level.cut,split_inv_fm=max(36.,28case.bg.T+2.))
            p=P.profile(k;mesh=level.mesh,tail_nodes=level.tail)
            append!(endpoints,endpoint_rows(case,ch,p,level.label;raw=level==last(LEVELS)))
            fine=p
            ch==:K_plus && push!(direct,direct_rows(case,ch,level))
        end
        h,b=old_rows(case,ch,oldsettings);append!(historical,h)
        append!(infrared,infrared_rows(case,ch,fine,b))
        for (name,rows) in (("endpoints",endpoints),("historical",historical),("direct",direct),("infrared",infrared))
            isempty(rows) || writecsv(joinpath(options.output,name*".csv"),rows)
        end
    end
    all(hashfile(path)==hash for (path,hash) in inputs) || error("frozen inputs changed during audit")
    all(hashfile(joinpath(ROOT,path))==hash for (path,hash) in hashes) || error("source changed during audit")
    W.writejson(joinpath(options.output,"manifest.json"),(schema="charged_gbu_q0_endpoint_audit_v1",
        git_head=readchomp(`git -C $ROOT rev-parse HEAD`),source_scan_run=SOURCE_RUN,source_scan_sha=SOURCE_SHA,
        actions_run_id=get(ENV,"GITHUB_RUN_ID","local"),case_count=3,channels=["K_plus","pi_plus","K_minus"],
        levels=LEVELS,static_imag_tolerance=IMAG_TOL,source_hashes=hashes,input_hashes=inputs,
        output_hashes=W.output_hashes(options.output),solver_called=false,new_backgrounds=0,
        full_density_evaluated=false,gate_tolerances_changed=false,production_authorized=false,
        status="diagnostic_probes_complete;not_a_condensation_verdict"))
    println("[q0-endpoint-audit] complete: $(options.output)")
end

function main(args=ARGS)
    options=parse_args(args)
    options===nothing && return println("Usage: audit_charged_gbu_q0_endpoint.jl --input-root <frozen-shards> --output <new-directory>")
    run_audit(options)
end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && ChargedGBUQ0EndpointAudit.main()
