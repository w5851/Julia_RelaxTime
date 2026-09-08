"""Opt-in research-production orchestration; never loaded by the default provider.

The verified infinite-thermal kernel remains single-source in analysis during
PR310 review. This explicit adapter is not promotion of the legacy MesonDensity
defaults. All runtime dependencies are snapshotted for reproducibility.
"""
module ChargedGBUResearchWorkflow
using CSV,JSON3,SHA,TOML,Distributed
using Main.Models
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(ROOT,"scripts","analysis","relaxtime","causal_gbu_infinite_qgate.jl"))
const Q=CausalGBUInfiniteQGate
const A=Q.A
const Y=Q.Y
const P=Q.P
const I=Q.I
const R=Q.R
const DEFAULT_CONFIG=joinpath(ROOT,"config","models","pnjl","charged_gbu_infinite_v1.toml")
hashfile(p)=bytes2hex(sha256(read(p)))
jsonsafe(x::AbstractFloat)=isfinite(x) ? x : nothing
jsonsafe(x::NamedTuple)=Dict(String(k)=>jsonsafe(v) for (k,v) in pairs(x))
jsonsafe(x::AbstractDict)=Dict(String(k)=>jsonsafe(v) for (k,v) in pairs(x))
jsonsafe(x::AbstractArray)=map(jsonsafe,x)
jsonsafe(x)=x
function writejson(p,x)
    temp=p*".partial"
    open(io->JSON3.write(io,jsonsafe(x)),temp,"w")
    mv(temp,p;force=true)
    open(io->print(io,hashfile(p)),p*".sha256","w")
end
function readjson(p)
    isfile(p*".sha256") && readchomp(p*".sha256")==hashfile(p) || error("checkpoint hash mismatch: $(p)")
    return JSON3.read(read(p,String))
end

"""Hash outputs without including the manifest or its previous checksum on resume."""
function output_hashes(output)
    return Dict(relpath(joinpath(root,f),output)=>hashfile(joinpath(root,f))
        for (root,_,files) in walkdir(output) for f in files
        if !occursin("source_snapshot",root) && !(f in ("manifest.json","manifest.json.sha256")))
end

function validate_config(c)
    for (key,value) in (("schema","charged_gbu_infinite_v1"),("thermal_target","infinity"),
        ("observable","fixed_quark_only_gbu_partial_yield"),("vacuum_regulator","two_line_Lambda"),
        ("background","FixedMuBConservedCharges"),("charge_to_baryon_ratio",0.4),
        ("strangeness_density_fm3",0.0),("meson_feedback",false),("production_default",false))
        c[key]==value || throw(ArgumentError("unsupported method contract: $(key)"))
    end
    n=c["numerics"];g=c["gates"]
    n["mesh"]>=512 && n["cut_nodes"]>=96 && n["omega_nodes"]>=96 && n["background_nodes"]>=48 ||
        throw(ArgumentError("numerics below accepted method floor"))
    orders=Int.(n["q_orders"])
    length(orders)>=2 && issorted(orders) && length(unique(orders))==length(orders) && first(orders)>=8 ||
        throw(ArgumentError("at least two increasing q orders required"))
    n["q_bulk_max_inv_fm"]==8 && n["q_tail_edges_inv_fm"]==[8,12,16,24,32] && n["tail_order"]>=8 ||
        throw(ArgumentError("unsupported tail bands"))
    for (key,limit) in (("kernel_absolute_inv_fm2",1e-6),("q_relative",.005),
        ("component_relative",.01),("inner_relative",.005),("eta_relative",.01),
        ("tail_fraction",1e-6),("background_residual",1e-7))
        0<g[key]<=limit || throw(ArgumentError("unsupported gate tolerance: $(key)"))
    end
    return c
end

function energy_grid(values)
    e=Float64.(values)
    !isempty(e) && all(x->isfinite(x) && x>0,e) && length(unique(e))==length(e) ||
        throw(ArgumentError("finite positive unique energies required"))
    return sort(e;rev=true)
end

function background(model,pt,c;seed=nothing)
    h=Main.Constants_PNJL.ħc_MeV_fm
    kw=(p_num=Int(c["numerics"]["background_nodes"]),t_num=8,
        residual_norm_max=c["gates"]["background_residual"],iterations=300,xi=0.)
    seed===nothing || (kw=merge(kw,(seed_guess=seed,)))
    eq=Models.solve(model,Models.FixedMuBConservedCharges(pt.muB_MeV/h,.4,0.),pt.T_MeV/h;kw...)
    eq.converged && eq.residual_norm<=c["gates"]["background_residual"] || error("equilibrium residual not accepted")
    state=Models.meanfield_state(eq.x_state)
    kernel=Main.RelaxTime.MesonInteractionKernel.build_full_kmt_interaction(state.phi;G=model.params.G_fm2,K=model.params.K_fm5)
    coupling=Dict(ch=>Main.RelaxTime.ChargedRPAKernel.charged_rpa_coupling(R.charged_rpa_spec(ch),kernel) for ch in R.CHANNELS)
    bg=(m=(u=Float64(eq.masses[1]),d=Float64(eq.masses[2]),s=Float64(eq.masses[3])),
        mu=(u=Float64(eq.mu_vec[1]),d=Float64(eq.mu_vec[2]),s=Float64(eq.mu_vec[3])),
        T=pt.T_MeV/h,Phi=Float64(state.Phi),PhiBar=Float64(state.PhiBar),coupling=coupling,
        vacuum=Float64(model.params.Λ_inv_fm),T_MeV=pt.T_MeV,muB_MeV=pt.muB_MeV,residual=Float64(eq.residual_norm))
    return (bg=bg,seed=Float64.(eq.solution),omega=Float64(eq.omega))
end

function restore_background(d)
    b=d.bg
    tri(x)=(u=Float64(x.u),d=Float64(x.d),s=Float64(x.s))
    return (m=tri(b.m),mu=tri(b.mu),T=Float64(b.T),Phi=Float64(b.Phi),PhiBar=Float64(b.PhiBar),
        coupling=Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(b.coupling)),vacuum=Float64(b.vacuum),
        T_MeV=Float64(b.T_MeV),muB_MeV=Float64(b.muB_MeV),residual=Float64(b.residual))
end

function profile(bg,ch,q,c;refined=false)
    n=c["numerics"]
    k=I.kernel(bg,ch,q;cut_nodes=Int(n["cut_nodes"]),split_inv_fm=max(36.,q+28bg.T+2.))
    return P.profile(k;mesh=Int(n["mesh"])*(refined ? 2 : 1),tail_nodes=refined ? 192 : 96)
end

function local_gate(bg,ch,c)
    p=profile(bg,ch,1.4,c;refined=true);k=p.kernel;g=c["gates"]
    s=Y.shell(p;nodes=192); coarse=Y.shell(profile(bg,ch,1.4,c);nodes=96)
    zs=[complex(k.threshold),complex((k.landau+k.threshold)/2),2+0.8im]
    raw_error=maximum(abs(P.polarization(p,z)-I.polarization(k,z;nodes=128)) for z in zs)
    radial_error=abs(I.radial_polarization(bg,ch,1.4,2+0.8im;nodes=192)-I.polarization(k,2+0.8im;nodes=128))
    ca=A.contour_audit(p;nodes=32);cb=A.contour_audit(p;nodes=64,radius=30.,indent=3e-6)
    values=[A.eta_shell(p,eta;lower=lo,nodes=96) for eta in (.001,.0003,.0001),lo in (.001,.0003)]
    errors=abs.(values.-s.density)/abs(s.density)
    eta_pass=all(errors[j+1,l]<=errors[j,l] for j in 1:2,l in 1:2) && maximum(errors[end,:])<g["eta_relative"] &&
        abs(values[end,1]-values[end,2])/abs(s.density)<g["eta_relative"]
    inner=abs(s.density-coarse.density)/abs(s.density)
    return (passed=raw_error<g["kernel_absolute_inv_fm2"] && radial_error<g["kernel_absolute_inv_fm2"] &&
        ca.passed && cb.passed && eta_pass && inner<g["inner_relative"],raw_error=raw_error,
        radial_error=radial_error,contour_passed=ca.passed && cb.passed,eta_passed=eta_pass,
        eta_relative=maximum(errors[end,:]),inner_relative=inner,scope="q1.4_representative_per_channel")
end

function reduce_bands(lo,hi,tails,g)
    n=hi.density+sum(r.density for r in tails)
    relative=abs(hi.density-lo.density)/abs(n)
    component=maximum(abs(getproperty(hi,f)-getproperty(lo,f)) for f in (:bound,:landau,:pair))/abs(n)
    tail=abs(last(tails).density)/abs(n)
    decreasing=all(abs(tails[j+1].density)<=abs(tails[j].density) for j in 1:length(tails)-1)
    omega=(hi.omega_tail+sum(r.omega_tail for r in tails))/abs(n)
    passed=isfinite(n) && n>0 && relative<g["q_relative"] && component<g["component_relative"] &&
        tail<g["tail_fraction"] && omega<g["tail_fraction"] && decreasing
    return (density=n,passed=passed,q_relative=relative,component_relative=component,
        final_tail_fraction=tail,omega_tail_fraction=omega,tail_decreasing=decreasing)
end

function channel_density(bg,ch,c,progress;checkpoint=identity)
    localcheck=local_gate(bg,ch,c)
    localcheck.passed || return (passed=false,status="local_gate_failed",density=NaN,local_gate=localcheck,rows=NamedTuple[])
    m=Q.mott_momenta(bg,ch)
    m.passed || error("Mott q split not converged")
    edges=sort!(unique!(vcat([0.,2bg.vacuum,8.],m.roots)))
    rows=NamedTuple[];n=c["numerics"]
    function band(pieces,order,label)
        start=length(rows)
        for (a,b) in pieces
            qs,ws=R.gauleg(a,b,order)
            for (q,w) in zip(qs,ws)
                p=profile(bg,ch,q,c);s=Y.shell(p;nodes=Int(n["omega_nodes"]))
                topology=label=="bulk" ? Q.topology(p) : (passed=true,)
                topology.passed || error("gap/cut topology failed q=$(q)")
                push!(rows,(q_inv_fm=q,weight=w,band=label,order=order,density=s.density,bound=s.bound,
                    landau=s.landau,pair=s.pair,root_count=s.root_count,negative_root_count=s.negative_root_count,
                    topology_evaluated=label=="bulk",omega_tail=s.omega_tail_conditional_bound))
            end
        end
        selected=rows[start+1:end]
        checkpoint(rows)
        field(f)=sum(r.weight*getproperty(r,f) for r in selected)
        progress("$(label) n=$(order) rows=$(length(rows))")
        return (density=field(:density),bound=field(:bound),landau=field(:landau),pair=field(:pair),omega_tail=field(:omega_tail))
    end
    tail_edges=n["q_tail_edges_inv_fm"]
    tails=[band([(a,b)],Int(n["tail_order"]),"tail$(j)") for (j,(a,b)) in enumerate(zip(tail_edges[1:end-1],tail_edges[2:end]))]
    previous=nothing;decision=nothing
    for order in Int.(n["q_orders"])
        current=band(collect(zip(edges[1:end-1],edges[2:end])),order,"bulk")
        if previous!==nothing
            decision=reduce_bands(previous,current,tails,c["gates"])
            decision.passed && return merge(decision,(status="accepted",local_gate=localcheck,mott=m.checks,
                final_order=order,rows=rows,bulk=current,tails=tails))
        end
        previous=current
    end
    return merge(decision,(status="q_convergence_failed",local_gate=localcheck,mott=m.checks,rows=rows))
end

function ratio_row(e,pt,results)
    function ratio(a,b)
        haskey(results,a) && haskey(results,b) || return (NaN,false)
        x,y=results[a],results[b]
        ok=x.passed && y.passed && isfinite(x.density) && isfinite(y.density) && x.density>=0 && y.density>0
        return (ok ? x.density/y.density : NaN,ok)
    end
    plus,pg=ratio("K_plus","pi_plus");minus,mg=ratio("K_minus","pi_minus")
    return (sqrt_s_NN_GeV=e,T_MeV=pt.T_MeV,muB_MeV=pt.muB_MeV,Kplus_over_pi_plus=plus,
        Kminus_over_pi_minus=minus,plus_passed=pg,minus_passed=mg,production_default=false)
end

function snapshot(output,config)
    paths=String[]
    for dir in (joinpath(ROOT,"src"),joinpath(ROOT,"config"),joinpath(ROOT,"scripts","analysis","relaxtime"))
        for (root,_,files) in walkdir(dir),f in files
            endswith(f,".jl") || endswith(f,".toml") || continue
            push!(paths,joinpath(root,f))
        end
    end
    push!(paths,joinpath(ROOT,"scripts","relaxtime","run_charged_gbu_freezeout_scan.jl"))
    push!(paths,joinpath(ROOT,"scripts","relaxtime","workflow","charged_gbu_plot.jl"))
    hashes=Dict(relpath(p,ROOT)=>hashfile(p) for p in unique(paths))
    for p in unique(paths)
        target=joinpath(output,"source_snapshot",relpath(p,ROOT));mkpath(dirname(target));cp(p,target)
    end
    return hashes
end

"""Single workflow owns mapping, BQS solves, four channels, gates and output.

Resume requires identical source/config/energy hashes; successful and failed
channel records are retained, never silently overwritten or reused across methods.
"""
function run_scan(;output,config=DEFAULT_CONFIG,energies=nothing,resume=false,make_plot=true,workers=3)
    workers in 1:4 || throw(ArgumentError("workers must be 1..4"))
    c=validate_config(TOML.parsefile(config))
    es=energy_grid(energies===nothing ? c["energies_GeV"] : energies)
    output=abspath(output); identity=bytes2hex(sha256(JSON3.write((config=c,energies=es))))
    statepath=joinpath(output,"run.json")
    if ispath(output)
        resume || error("output exists; explicit --resume required")
        state=readjson(statepath);state.identity==identity || error("resume settings mismatch")
        all(hashfile(joinpath(ROOT,String(p)))==String(h) for (p,h) in pairs(state.source_hashes)) || error("resume source mismatch")
        hashes=Dict(String(p)=>String(h) for (p,h) in pairs(state.source_hashes))
    else
        mkpath(output);hashes=snapshot(output,config)
        writejson(statepath,(identity=identity,config=c,energies_GeV=es,source_hashes=hashes,
            git_head=readchomp(`git -C $ROOT rev-parse HEAD`),route="charged_gbu_infinite_v1",production_default=false))
    end
    points=Models.build_freezeout_scan_points(es;profile=Models.load_freezeout_profile(profile=c["freezeout_profile"]),traversal=:sqrts_descending)
    model=Models.create_model(:PNJL);jobs=[];seed=nothing
    for pt in points
        dir=joinpath(output,"energy_$(pt.sqrt_s_NN_GeV)");mkpath(dir);path=joinpath(dir,"background.json")
        try
            if !isfile(path)
                b=background(model,pt,c;seed=seed);writejson(path,b)
            end
            b=readjson(path);bg=restore_background(b);seed=Float64.(b.seed)
            for ch in R.CHANNELS;push!(jobs,(bg=bg,ch=ch,dir=dir,energy=pt.sqrt_s_NN_GeV));end
            println("[charged-gbu-background] s=$(pt.sqrt_s_NN_GeV) residual=$(bg.residual)");flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            writejson(joinpath(dir,"background_failure.json"),(reason=sprint(showerror,err),))
        end
    end
    # Process isolation is required: wide-coordinate fallback uses setprecision.
    pending=filter(job->!isfile(joinpath(job.dir,"$(job.ch).json")),jobs)
    if workers==1
        foreach(job->run_job(job,c),pending)
    elseif !isempty(pending)
        ids=addprocs(workers;exeflags=`--project=$ROOT --threads=1`,dir=ROOT)
        try
            for id in ids
                remotecall_wait(Base.include,id,Main,joinpath(ROOT,"src","models","Models.jl"))
                remotecall_wait(Base.include,id,Main,@__FILE__)
            end
            pmap(job->Main.ChargedGBUResearchWorkflow.run_job(job,c),WorkerPool(ids),pending)
        finally
            rmprocs(ids)
        end
    end
    ratios=NamedTuple[];densities=NamedTuple[]
    for pt in reverse(points)
        dir=joinpath(output,"energy_$(pt.sqrt_s_NN_GeV)");results=Dict{String,Any}()
        for ch in R.CHANNELS
            path=joinpath(dir,"$(ch).json")
            r=isfile(path) ? readjson(path) : nothing
            results[String(ch)]=(passed=r!==nothing && r.passed,density=r===nothing || r.density===nothing ? NaN : Float64(r.density))
            push!(densities,merge((sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,channel=String(ch),
                status=r===nothing ? "background_failed" : String(r.status),
                reason=r!==nothing && hasproperty(r,:reason) ? String(r.reason) : ""),results[String(ch)]))
        end
        push!(ratios,ratio_row(pt.sqrt_s_NN_GeV,pt,results))
    end
    CSV.write(joinpath(output,"ratios.csv"),ratios);CSV.write(joinpath(output,"densities.csv"),densities)
    all(hashfile(joinpath(ROOT,p))==h for (p,h) in hashes) || error("source changed during scan")
    passed=all(r.plus_passed && r.minus_passed for r in ratios)
    record=Dict("status"=>passed ? "complete_research_curve_accepted" : "complete_scan_with_failed_points",
        "point_count"=>length(es),"accepted_ratio_pairs"=>count(r->r.plus_passed && r.minus_passed,ratios),
        "config_identity"=>identity,"source_hashes"=>hashes,"production_default"=>false,
        "thermal_target"=>"infinity","meson_feedback"=>false,"q_tail_global_analytic_proof"=>false,
        "kernel_dependency"=>"single_source_reviewed_analysis_modules_opt_in_adapter")
    if make_plot
        Base.include(@__MODULE__,joinpath(ROOT,"scripts","relaxtime","workflow","charged_gbu_plot.jl"))
        Base.invokelatest(() -> render_ratio(ratios,output))
    end
    record["output_hashes"]=output_hashes(output)
    writejson(joinpath(output,"manifest.json"),record)
    return (output=output,status=record["status"],accepted=record["accepted_ratio_pairs"],points=length(es))
end

function run_job(job,c)
    path=joinpath(job.dir,"$(job.ch).json")
    progress=s->begin println("[charged-gbu] s=$(job.energy) $(job.ch) $(s)");flush(stdout);end
    result=try
        channel_density(job.bg,job.ch,c,progress;
            checkpoint=rows->CSV.write(joinpath(job.dir,"$(job.ch)_shells.csv"),rows))
    catch err
        err isa InterruptException && rethrow()
        (passed=false,status="evaluation_failed",density=NaN,reason=sprint(showerror,err))
    end
    writejson(path,result);progress("status=$(result.status)")
    return result.status
end
end
