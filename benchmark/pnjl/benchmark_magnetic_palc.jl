#!/usr/bin/env julia

"""
    benchmark_magnetic_palc.jl magnetic|nonmagnetic OUTPUT [WARM_REPETITIONS]

Same-host warm benchmarks for the validated magnetic PALC probe and the current
nonmagnetic production path. Run each mode in its own source checkout/process.
The root environment supplies all dependencies. BLAS and Julia use one thread.

Magnetic timings separate anchor, continuation, target correction, target output
and all-vertex diagnostic output. The latter is excluded from solve_seconds.
The first full pass is warmup; subsequent passes rotate cases. Verification,
profiling, JSON writing and pre-pass GC are outside the timed stages; GC caused
by the stages remains included. No source algorithm or solver default is changed.

MAGNETIC_REFERENCE_RESULT may identify the existing p128 validation JSON.
PALC_BENCH_SOURCE_SHA and PALC_BENCH_CPU_AFFINITY document the source and host.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BENCH_MODE = isempty(ARGS) ? "" : ARGS[1]
if BENCH_MODE == "magnetic"
    include(joinpath(BENCH_ROOT, "scripts", "analysis", "pnjl_magnetic", "compare_single_anchor_palc.jl"))
elseif BENCH_MODE == "nonmagnetic"
    include(joinpath(BENCH_ROOT, "scripts", "analysis", "pnjl_bifurcation_spike", "run_path_reuse_comparison.jl"))
else
    error("expected magnetic|nonmagnetic OUTPUT [WARM_REPETITIONS]")
end

module MagneticPALCBenchmark

using Dates, JSON, LinearAlgebra, Statistics, StaticArrays, ForwardDiff, Profile, SHA
import BifurcationKit
import BifurcationKit: @optic

const Models = Main.Models
const HBARC = 197.3269804
const TARGETS = [872.08, 900.0, 1000.0, 1050.0, 1100.0, 1150.0, 1200.0, 1500.0]

# BifurcationKit deep-copies the iterator/problem for both continuation directions.
# Share only the diagnostic counter, so counts cover both directions. This method
# exists solely in this benchmark process and changes no numerical state or policy.
Base.deepcopy_internal(work::Models.PALCWorkTelemetry, ::IdDict) = work

finite_json(x::AbstractFloat) = isfinite(x) ? x : nothing
finite_json(x::AbstractArray) = map(finite_json, x)
finite_json(x::NamedTuple) = Dict(string(k) => finite_json(v) for (k,v) in pairs(x))
finite_json(x::AbstractDict) = Dict(string(k) => finite_json(v) for (k,v) in pairs(x))
finite_json(x) = x

function write_result(path, result)
    mkpath(dirname(abspath(path)))
    open(path * ".pending", "w") do io
        JSON.print(io, finite_json(result), 2)
        println(io)
    end
    mv(path * ".pending", path; force=true)
end

function source_hashes()
    files = ["Project.toml", "Manifest.toml", "benchmark/pnjl/benchmark_magnetic_palc.jl",
        "src/models/solver/path/PALCBackend.jl", "src/models/solver/path/PathContinuation.jl",
        "src/models/pnjl_physics/PNJLModel.jl", "src/models/pnjl_physics/core/PNJLIntegrals.jl",
        "src/models/pnjl_physics/PNJLMagneticModel.jl", "src/models/pnjl_physics/MagneticGapSolver.jl",
        "src/models/pnjl_physics/core/MagneticIntegrals.jl",
        "src/models/pnjl_physics/core/MagneticThermodynamics.jl"]
    script = Main.BENCH_MODE == "magnetic" ?
        "scripts/analysis/pnjl_magnetic/compare_single_anchor_palc.jl" :
        "scripts/analysis/pnjl_bifurcation_spike/run_path_reuse_comparison.jl"
    push!(files, script)
    return Dict(file => bytes2hex(sha256(read(joinpath(Main.BENCH_ROOT,file))))
        for file in files if isfile(joinpath(Main.BENCH_ROOT,file)))
end

function metadata()
    return Dict{String,Any}(
        "schema" => "magnetic_palc_performance_v1", "mode" => Main.BENCH_MODE,
        "generated_at" => string(now()), "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "cpu" => Sys.CPU_NAME, "cpu_affinity" => get(ENV,"PALC_BENCH_CPU_AFFINITY","unspecified"),
        "source_sha" => get(ENV,"PALC_BENCH_SOURCE_SHA","working_copy"),
        "telemetry_contract" => "Diagnostic work counters are shared across BifurcationKit deepcopy; both continuation directions are counted without changing numerical state.",
        "source_hashes_before" => source_hashes(),
        "timing_contract" => "Warm same-process stage timing; first full pass excluded; cases rotated. Pre-pass GC, verification, profiling and JSON excluded; stage GC included. Magnetic all-vertex diagnostic output is reported separately. Nonmagnetic uses current production solve_path plus selected_path_results.",
        "runs" => Any[], "microbenchmarks" => Any[])
end

timing(run) = (;seconds=run.time, allocated_bytes=run.bytes, gc_seconds=run.gctime)

function measure(f; target_seconds=0.15, repetitions=3)
    f()
    trial = @elapsed f()
    samples = clamp(ceil(Int,target_seconds/max(trial,1e-7)),1,10000)
    rows = NamedTuple[]
    for _ in 1:repetitions
        GC.gc()
        run = @timed for _ in 1:samples
            f()
        end
        push!(rows,(;seconds=run.time/samples,allocated_bytes=run.bytes/samples,
            gc_seconds=run.gctime/samples))
    end
    return (;samples,repetitions,seconds=median(r.seconds for r in rows),
        allocated_bytes=median(r.allocated_bytes for r in rows),
        gc_seconds=median(r.gc_seconds for r in rows),
        min_seconds=minimum(r.seconds for r in rows),max_seconds=maximum(r.seconds for r in rows))
end

Base.@kwdef mutable struct ContinuationWork
    scalar_residual_calls::Int = 0
    dual_residual_calls::Int = 0
    jacobian_calls::Int = 0
end
Base.deepcopy_internal(work::ContinuationWork, ::IdDict) = work

# Mirrors the existing probe's controller exactly, adding counters only.
function counted_palc(model,cfg,anchor,work)
    raw = Main._residual(model,cfg)
    residual = function(x,p)
        if eltype(x) <: ForwardDiff.Dual
            work.dual_residual_calls += 1
        else
            work.scalar_residual_calls += 1
        end
        raw(x,p)
    end
    p0 = (;muq_fm=Main._muq_fm(cfg.muB_max_MeV))
    jac = function(x,p)
        work.jacobian_calls += 1
        cfg.jacobian_method === :forward ?
            ForwardDiff.jacobian(y->residual(y,p),Float64.(x)) : Main._jacobian(residual,x,p)
    end
    problem = BifurcationKit.BifurcationProblem(residual,collect(anchor),p0,(@optic _.muq_fm);
        J=jac,record_from_solution=(x,p;kwargs...)->(;solution=Float64.(x)))
    options = BifurcationKit.ContinuationPar(
        p_min=Main._muq_fm(cfg.muB_min_MeV),p_max=Main._muq_fm(cfg.muB_max_MeV),
        ds=cfg.ds,dsmax=cfg.dsmax,max_steps=cfg.max_steps,detect_bifurcation=0,
        newton_options=BifurcationKit.NewtonPar(tol=cfg.newton_tol,
            max_iterations=cfg.newton_max_iterations,verbose=false,linesearch=cfg.linesearch))
    return BifurcationKit.continuation(problem,BifurcationKit.PALC(),options;bothside=true,verbosity=0)
end

function probe_config(p_num; T_MeV=5.0,n_max=64)
    return Main.ProbeConfig(;T_MeV,p_num,n_max,pz_max=40.0,ds=.002,dsmax=.006,
        max_steps=1800,anchor_side=:high,jacobian_method=:forward,linesearch=true,
        newton_tol=1e-10,newton_max_iterations=28,compare_catalog=false,target_muB_MeV=TARGETS)
end

function correct_targets(model,cfg,branch)
    targets = NamedTuple[]
    failures = NamedTuple[]
    solves=0; attempts=0; failed_attempts=0; iterations=0
    for muB in cfg.target_muB_MeV
        roots = Models.MagneticGapCandidate[]
        for seed in Main._crossing_seeds(branch,Main._muq_fm(muB))
            solves += 1
            result = try
                Main._solve_one(model,Main._T_fm(cfg.T_MeV),muB,seed,cfg)
            catch err
                err isa InterruptException && rethrow()
                push!(failures,(;muB_MeV=muB,message=sprint(showerror,err)))
                nothing
            end
            result === nothing && continue
            attempts += result.attempt_count
            failed_attempts += result.failed_attempts
            for candidate in result.candidates
                iterations += candidate.iterations
                any(norm(candidate.x_state-r.x_state)<=cfg.root_merge_tol for r in roots) || push!(roots,candidate)
            end
        end
        push!(targets,(;muB_MeV=muB,roots))
    end
    return (;targets,solves,attempts,failed_attempts,iterations,failures)
end

function verify_magnetic(reference,cfg,summary,targets)
    expected = reference["targets"]
    distances = Float64[]
    counts = Int[]
    for (actual,expected_target) in zip(targets,expected)
        actual.muB_MeV == expected_target["muB_MeV"] || error("target ordering changed")
        baseline = expected_target["palc_roots"]
        length(actual.roots) == length(baseline) || error("magnetic root count changed at $(actual.muB_MeV)")
        push!(counts,length(actual.roots))
        for root in actual.roots
            push!(distances,minimum(norm(root.x_state-r["x_state"]) for r in baseline))
        end
        for root in baseline
            push!(distances,minimum(norm(root["x_state"]-r.x_state) for r in actual.roots))
        end
    end
    summary.folds == reference["branch"]["folds"] || error("magnetic fold count changed")
    summary.muB_min_MeV <= cfg.muB_min_MeV+1e-6 || error("lower endpoint not reached")
    summary.muB_max_MeV >= cfg.muB_max_MeV-1e-6 || error("upper endpoint not reached")
    maximum(distances) <= 1e-6 || error("magnetic reference roots changed")
    return (;root_counts=counts,total_roots=sum(counts),reference_max_state_distance=maximum(distances),
        residual_max=maximum(r.residual_norm for t in targets for r in t.roots),
        all_physical=all(r.physical for t in targets for r in t.roots))
end

function magnetic_run(eB,p_num,repetition,reference)
    cfg = probe_config(p_num)
    model = Main._model(eB,:current_mfir)
    work = ContinuationWork()
    GC.gc()
    anchor = @timed Main._solve_one(model,Main._T_fm(cfg.T_MeV),cfg.muB_max_MeV,Models.QUARK_SEED_5,cfg)
    root = anchor.value.candidates[anchor.value.selected_index]
    continuation = @timed counted_palc(model,cfg,root.x_state,work)
    correction = @timed correct_targets(model,cfg,continuation.value)
    target_reconstruction = @timed [(;muB_MeV=t.muB_MeV,
        roots=[Main._state_record(model,r.x_state,t.muB_MeV,cfg) for r in t.roots])
        for t in correction.value.targets]
    vertex_reconstruction = @timed [Main._state_record(model,p.solution,Float64(p.param)*3HBARC,cfg)
        for p in continuation.value.branch]
    summary = Main._branch_summary(continuation.value,cfg)
    validation = verify_magnetic(reference,cfg,summary,correction.value.targets)
    isempty(correction.value.failures) || error("failed magnetic target correction")
    stages = (;anchor=timing(anchor),continuation=timing(continuation),
        correction=timing(correction),target_reconstruction=timing(target_reconstruction),
        vertex_reconstruction=timing(vertex_reconstruction))
    solve_seconds = sum(s.seconds for s in (stages.anchor,stages.continuation,stages.correction,stages.target_reconstruction))
    record = (;name="magnetic_e$(eB)_p$(p_num)",eB_GeV2=eB,T_MeV=cfg.T_MeV,p_num,
        n_max=cfg.n_max,zeta_num=model.magnetic.zeta_num,repetition,warmup=repetition==0,
        solve_seconds,total_with_vertex_seconds=solve_seconds+stages.vertex_reconstruction.seconds,stages,
        branch=summary,validation,work=(;scalar_residual_calls=work.scalar_residual_calls,
            dual_residual_calls=work.dual_residual_calls,jacobian_calls=work.jacobian_calls,
            anchor_iterations=root.iterations,anchor_attempts=anchor.value.attempt_count,
            correction_solves=correction.value.solves,correction_attempts=correction.value.attempts,
            correction_failed_attempts=correction.value.failed_attempts,correction_iterations=correction.value.iterations))
    println("DONE ",record.name," repeat=",repetition," solve_s=",solve_seconds,
        " points=",summary.points," folds=",summary.folds," roots=",validation.total_roots);flush(stdout)
    return record,(;model,cfg,anchor=root.x_state,branch=continuation.value)
end

function potential_part(model,x,T,mu,cfg,part)
    thermo = Models.MagneticThermodynamics
    masses = Models.calculate_mass_vec(model,SVector(x[1],x[2],x[3]))
    params = model.base.params
    if part === :thermal
        return sum(thermo.omegat_flavor_landau(masses[i],mu[i],T,x[4],x[5],
            thermo.QUARK_CHARGE_ABS[i],model.magnetic.eB_fm2;
            n_max=cfg.n_max,p_num=cfg.p_num,pz_max=cfg.pz_max) for i in 1:3)
    elseif part === :magnetic_vacuum
        return sum(thermo.omega_magnetic_mfir(masses[i],thermo.QUARK_CHARGE_ABS[i],
            model.magnetic.eB_fm2;zeta_num=model.magnetic.zeta_num,N_c=params.N_color) for i in 1:3)
    elseif part === :zero_field_vacuum
        return sum(-2params.N_color*thermo.PNJLCore.vacuum_integral_with_cutoff(
            masses[i],convert(typeof(masses[i]),params.Λ_inv_fm)) for i in 1:3)
    elseif part === :local_potentials
        return Models.calculate_chiral(model,SVector(x[1],x[2],x[3]))+
            Models.polyakov_potential(model,x[4],x[5],T)
    end
    error("unknown potential part")
end

function active_nodes(model,x,T,mu,cfg)
    masses = Models.calculate_mass_vec(model,SVector(x[1],x[2],x[3]))
    charges = Models.MagneticThermodynamics.QUARK_CHARGE_ABS
    return [begin
        active=0;segments=0
        for n in 0:cfg.n_max
            mn = sqrt(masses[f]^2+2n*charges[f]*model.magnetic.eB_fm2+1e-24)
            abs(mu[f])+40T <= mn && continue
            active += 1
            lo=0.0;tmax=asinh(cfg.pz_max/mn)
            for width in (-12,-6,0,6,12,40)
                energy=abs(mu[f])+width*T
                energy<=mn && continue
                hi=min(asinh(sqrt(energy^2-mn^2)/mn),tmax)
                hi>lo || continue
                segments+=1;lo=hi
                hi>=tmax && break
            end
        end
        (;flavor=f,active_levels=active,segments,nodes=segments*cfg.p_num)
    end for f in 1:3]
end

function magnetic_micro(model,cfg,x_input,muB;label)
    x = collect(Float64.(x_input))
    T=Main._T_fm(cfg.T_MeV);mu=Models.normalize_mu_vec(Main._muq_fm(muB));p=(;muq_fm=mu[1])
    residual=Main._residual(model,cfg)
    gap=y->residual(y,p)
    jac=y->ForwardDiff.jacobian(gap,y)
    fd=()->Main._jacobian(residual,x,p)
    rows=Dict{String,Any}("residual"=>measure(()->gap(x)),
        "jacobian_ad"=>measure(()->jac(x)),"jacobian_fd"=>measure(fd),
        "output_state_record"=>measure(()->Main._state_record(model,x,muB,cfg)))
    parts=(:thermal,:magnetic_vacuum,:zero_field_vacuum,:local_potentials)
    summed=zeros(5)
    for part in parts
        f=y->potential_part(model,y,T,mu,cfg,part)
        g=y->ForwardDiff.gradient(f,y)
        summed .+= g(x)
        rows[string(part)*"_gradient"]=measure(()->g(x))
        rows[string(part)*"_hessian"]=measure(()->ForwardDiff.jacobian(g,x))
    end
    split_error=norm(summed-gap(x))
    split_error<=1e-9 || error("potential decomposition mismatch: $(split_error)")
    a=jac(x);d=fd()
    return (;label,T_MeV=cfg.T_MeV,muB_MeV=muB,eB_GeV2=model.magnetic.eB_fm2*(HBARC/1000)^2,
        p_num=cfg.p_num,n_max=cfg.n_max,zeta_num=model.magnetic.zeta_num,
        x_state=x,active_work=active_nodes(model,x,T,mu,cfg),measurements=rows,
        component_gradient_error=split_error,jacobian_fd_ad_max_abs=maximum(abs,a-d),
        jacobian_fd_ad_relative=norm(a-d)/max(norm(a),eps()))
end

function profile_summary(f)
    f()
    Profile.init(n=5_000_000,delay=0.001)
    Profile.clear()
    @profile f()
    data,lookup=Profile.retrieve()
    counts=Dict("thermal"=>0,"zeta"=>0,"polyakov_distribution"=>0,
        "atan_series"=>0,"forwarddiff"=>0,"linear_algebra"=>0)
    total=0
    frames=String[]
    for ip in data
        if ip==0
            if !isempty(frames)
                total+=1
                any(s->occursin("_thermal_level_pair",s),frames) && (counts["thermal"]+=1)
                any(s->occursin("zeta_prime_minus_one",s),frames) && (counts["zeta"]+=1)
                any(s->occursin("_scaled_polyakov_terms",s)||occursin("_polyakov_log_and_net_density",s),frames) && (counts["polyakov_distribution"]+=1)
                any(s->occursin("_atan_series_small",s),frames) && (counts["atan_series"]+=1)
                any(s->occursin("ForwardDiff",s),frames) && (counts["forwarddiff"]+=1)
                any(s->occursin("LinearAlgebra",s),frames) && (counts["linear_algebra"]+=1)
            end
            empty!(frames)
        elseif haskey(lookup,ip)
            for frame in lookup[ip]
                push!(frames,string(frame.func," ",frame.file))
            end
        end
    end
    return (;sampled_stacks=total,inclusive_stack_counts=counts,
        inclusive_fractions=Dict(k=>v/max(total,1) for (k,v) in counts),
        note="Inclusive stack presence; categories overlap. Windows sampling interval is approximate.")
end

function run_magnetic(output,repetitions)
    reference_path=get(ENV,"MAGNETIC_REFERENCE_RESULT",joinpath(Main.BENCH_ROOT,
        "data","outputs","results","pnjl","magnetic_palc_validation","20261112","corrected_palc_p128.json"))
    reference=JSON.parsefile(reference_path)
    data=metadata()
    data["reference_sha256"]=bytes2hex(sha256(read(reference_path)))
    data["controls"]=(;T_MeV=5,muB_range_MeV=[500,1700],targets_muB_MeV=TARGETS,
        n_max=64,pz_max=40,ds=.002,dsmax=.006,max_steps=1800,newton_tol=1e-10,
        newton_max_iterations=28,jacobian="forward",linesearch=true)
    cases=[(.05,128),(.155,128),(.4,128),(.155,64),(.4,64)]
    profile_context=nothing
    for repetition in 0:repetitions
        for (eB,p_num) in circshift(cases,repetition)
            expected=only(filter(c->c["eB_GeV2"]==eB,reference["cases"]))
            println("START magnetic eB=",eB," p=",p_num," repeat=",repetition);flush(stdout)
            record,context=magnetic_run(eB,p_num,repetition,expected)
            push!(data["runs"],record)
            eB==.155 && p_num==128 && (profile_context=context)
            write_result(output,data)
        end
    end
    for eB in (.05,.155,.4)
        expected=only(filter(c->c["eB_GeV2"]==eB,reference["cases"]))
        target=only(filter(t->t["muB_MeV"]==1500.,expected["targets"]))
        x=only(target["palc_roots"])["x_state"]
        model=Main._model(eB,:current_mfir)
        push!(data["microbenchmarks"],magnetic_micro(model,probe_config(128),x,1500.;label="cold_dense_e$(eB)"))
        write_result(output,data)
    end
    expected=only(filter(c->c["eB_GeV2"]==.155,reference["cases"]))
    roots=only(filter(t->t["muB_MeV"]==1000.,expected["targets"]))["palc_roots"]
    root=roots[argmax(r["masses_MeV"][1] for r in roots)]
    model=Main._model(.155,:current_mfir)
    push!(data["microbenchmarks"],magnetic_micro(model,probe_config(128),root["x_state"],1000.;label="cold_large_mass_e0.155"))
    write_result(output,data)
    for T_MeV in (120.,200.)
        muB=930.
        n_max=Models.MagneticThermodynamics.resolve_magnetic_nmax(T_MeV/HBARC,
            Models.normalize_mu_vec(muB/(3HBARC)),model.magnetic)
        cfg=probe_config(128;T_MeV,n_max)
        solved=Main._solve_one(model,T_MeV/HBARC,muB,Models.QUARK_SEED_5,cfg)
        x=solved.candidates[solved.selected_index].x_state
        push!(data["microbenchmarks"],magnetic_micro(model,cfg,x,muB;label="warm_T$(Int(T_MeV))_e0.155"))
        write_result(output,data)
    end
    c=profile_context
    data["continuation_profile"]=profile_summary(()->Main._palc(c.model,c.cfg,c.anchor,c.cfg.muB_max_MeV))
    data["source_hashes_after"]=source_hashes()
    data["source_unchanged"]=data["source_hashes_before"]==data["source_hashes_after"]
    data["complete"]=true
    write_result(output,data)
    return data
end

function compact_nonmagnetic(run,path,name,repetition)
    roots=Main.PALCEndpointComparison.all_root_records(Models.create_model(:PNJL),path,run.result)
    all(p!==nothing && p.converged for p in run.selected) || error("nonmagnetic selected solve failed")
    counts=[length(t.roots) for t in roots]
    name===:tmu_t120_fine && sum(counts)!=111 && error("T120 fine window did not recover the known 111 roots")
    return (;name,repetition,warmup=repetition==0,T_MeV=path.T_fm*HBARC,
        target_count=length(path.mu_values),root_counts=counts,total_roots=sum(counts),
        seconds=run.seconds,allocated_bytes=run.allocated_bytes,gc_seconds=run.gc_seconds,
        work=run.work,residual_max=maximum(r.residual for t in roots for r in t.roots),
        diagnostics=Models.to_namedtuple(run.result.diagnostics))
end

function nonmagnetic_micro(model,T_MeV)
    value=310/HBARC
    path=Models.FixedMuPath(T_MeV/HBARC,[value];p_num=24,t_num=8)
    options=Models._path_options(model,path;iterations=300)
    residual=Models._path_residual(model,path,options)
    seed=first(Models._endpoint_anchor_attempts(path)).seed
    root=Models._path_root(model,path,value,seed,residual,options,Models.PALCWorkTelemetry())
    root!==nothing && root.converged || error("nonmagnetic microbenchmark anchor failed")
    x=root.solution;out=similar(x)
    ad=Models._palc_jacobian(residual,x,value,Models.PALCWorkTelemetry(),:forwarddiff)
    fd=Models._palc_jacobian(residual,x,value,Models.PALCWorkTelemetry(),:finite_difference)
    a=copy(ad(x,(;value)));d=copy(fd(x,(;value)))
    return (;name="matched_T$(Int(T_MeV))_muq310",T_MeV,muB_MeV=930.,
        p_num=24,t_num=8,residual_norm=norm(residual(x,value)),measurements=(
            residual=measure(()->residual(x,value)),
            residual_inplace=measure(()->Models._path_residual!(out,residual,x,value)),
            jacobian_ad=measure(()->ad(x,(;value))),
            jacobian_fd=measure(()->fd(x,(;value)))),
        jacobian_fd_ad_relative=norm(a-d)/norm(a))
end

function run_nonmagnetic(output,repetitions)
    ENV["PALC_TIMING_SAMPLES"]="5"
    model=Models.create_model(:PNJL)
    cases=filter(c->c[1] in (:tmu_t120_xi0,:tmu_t120_fine,:tmu_t200_xi0),
        Main.PALCEndpointComparison.anchor_comparison_cases())
    data=metadata();data["timing_samples_per_pass"]=5
    data["controls"]=(;p_num=24,t_num=8,iterations=300,strategy="current default PALC",
        target_predictor="cubic",jacobian="finite_difference")
    for repetition in 0:repetitions
        for (name,path) in circshift(cases,repetition)
            println("START nonmagnetic ",name," repeat=",repetition);flush(stdout)
            run=Main.PALCPathReuseComparison.run_method(model,path,:default)
            record=compact_nonmagnetic(run,path,name,repetition)
            push!(data["runs"],record)
            println("DONE ",name," seconds=",run.seconds," roots=",record.total_roots);flush(stdout)
            write_result(output,data)
        end
    end
    for (name,path) in cases
        run=Main.PALCPathReuseComparison.run_method(model,path,:default)
        push!(data["microbenchmarks"],(;name,
            probe=Main.PALCPathReuseComparison.residual_microbenchmark(model,path,run.result;samples=1000)))
    end
    for T_MeV in (120.,200.)
        push!(data["microbenchmarks"],nonmagnetic_micro(model,T_MeV))
    end
    data["source_hashes_after"]=source_hashes()
    data["source_unchanged"]=data["source_hashes_before"]==data["source_hashes_after"]
    data["complete"]=true
    write_result(output,data)
    return data
end

function main(args)
    length(args)>=2 || error("expected MODE OUTPUT [WARM_REPETITIONS]")
    repetitions=length(args)>=3 ? parse(Int,args[3]) : 3
    1<=repetitions<=6 || throw(ArgumentError("warm repetitions must be 1..6"))
    BLAS.set_num_threads(1)
    return args[1]=="magnetic" ? run_magnetic(args[2],repetitions) : run_nonmagnetic(args[2],repetitions)
end

end # module

MagneticPALCBenchmark.main(ARGS)
