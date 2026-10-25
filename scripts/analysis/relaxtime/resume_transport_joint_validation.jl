#!/usr/bin/env julia
# Recover completed checkpoints into a new run; never rewrite the interrupted run.
include(joinpath(@__DIR__, "run_transport_joint_validation.jl"))

checkpoint_rows(directory, name) = isfile(joinpath(directory,name)) ?
    NamedTuple[NamedTuple(row) for row in CSV.File(joinpath(directory,name))] : NamedTuple[]

function verify_joint_checkpoint(directory)
    manifest=JSON.parsefile(joinpath(directory,"manifest.json"))
    isempty(checkpoint_rows(directory,"failures.csv")) ||
        error("checkpoint contains numerical failures; investigate them before recovery")
    for file in manifest["source_files"]
        path=replace(file["path"],'\\'=>'/')
        (startswith(path,"src/") || startswith(path,"config/models/") ||
            startswith(path,"config/physics/") || path in ("Project.toml","Manifest.toml") ||
            (startswith(path,"scripts/analysis/relaxtime/") &&
             path!="scripts/analysis/relaxtime/run_transport_joint_validation.jl")) || continue
        bytes2hex(sha256(read(joinpath(ROOT,file["path"]))))==file["sha256"] ||
            error("checkpoint numerical source changed: $path")
    end
    current=replace(read(joinpath(@__DIR__,"run_transport_joint_validation.jl"),String),"\r\n"=>"\n")
    original=replace(read(joinpath(directory,"source","scripts","analysis","relaxtime",
        "run_transport_joint_validation.jl"),String),"\r\n"=>"\n")
    for name in ("joint_background","joint_transport","independent_joint_rate","joint_observables")
        pattern=Regex("(?ms)^function "*name*"\\(.*?(?=^function |\\z)")
        a,b=match(pattern,current),match(pattern,original)
        a!==nothing && b!==nothing && strip(a.match)==strip(b.match) ||
            error("checkpoint implementation changed: $name")
    end
    for file in manifest["input_files"]
        bytes2hex(sha256(read(joinpath(directory,file["snapshot"]))))==file["sha256"] ||
            error("checkpoint input snapshot changed")
    end
    configs=manifest["parameters"]["configurations"]
    configs["reference"]==JSON.parse(JSON.json(joint_config("reference"))) ||
        error("checkpoint reference configuration changed")
    candidate=merge(joint_config("candidate"),(p=Int(configs["candidate"]["p"]),))
    configs["candidate"]==JSON.parse(JSON.json(candidate)) || error("checkpoint candidate configuration changed")
    names=(Tuple(Symbol(:tau_,s) for s in SPECIES)...,:eta_over_s,:zeta_over_s,:sigma_over_T)
    required=Set(String.(Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
    rates=checkpoint_rows(directory,"rates.csv")
    observables=checkpoint_rows(directory,"observables.csv")
    comparisons=checkpoint_rows(directory,"comparisons.csv")
    backgrounds=checkpoint_rows(directory,"background.csv")
    length(Set((r.point_id,r.configuration,r.channel) for r in rates))==length(rates) ||
        error("duplicate checkpoint rate")
    length(Set((r.point_id,r.configuration) for r in observables))==length(observables) ||
        error("duplicate checkpoint observables")
    length(Set((r.point_id,r.observable) for r in comparisons))==length(comparisons) ||
        error("duplicate checkpoint comparison")
    for r in rates
        r.configuration in ("candidate","reference") && r.channel in required || error("unknown checkpoint rate")
        tolerance=r.configuration=="candidate" ? 1e-3 : 2e-4
        isfinite(r.rate) && r.rate>0 && isfinite(r.estimated_or_reference_l1) &&
            0<=r.estimated_or_reference_l1<=tolerance || error("unqualified checkpoint rate")
    end
    for row in observables
        selected=filter(r->r.point_id==row.point_id && r.configuration==row.configuration,rates)
        length(selected)==length(required) && Set(r.channel for r in selected)==required ||
            error("checkpoint observable lacks all 17 rates")
        all(v->isfinite(v) && v>0,(getproperty(row,k) for k in names)) &&
            isfinite(row.density_relative_error) && 0<=row.density_relative_error<=1e-3 ||
            error("invalid checkpoint observable or density")
    end
    for row in comparisons
        name=Symbol(row.observable);name in names || error("unknown checkpoint observable")
        c=only(r for r in observables if r.point_id==row.point_id && r.configuration=="candidate")
        r=only(r for r in observables if r.point_id==row.point_id && r.configuration=="reference")
        row.candidate==getproperty(c,name) && row.reference==getproperty(r,name) ||
            error("checkpoint comparison does not match observables")
        change=row.candidate/row.reference-1
        isapprox(change,row.relative_change;rtol=1e-12,atol=1e-14) &&
            row.passed==(abs(change)<=0.01) || error("checkpoint qualification flag is inconsistent")
    end
    points=Set(String(r.point_id) for r in backgrounds)
    length(points)==manifest["parameters"]["points"] ||
        error("some checkpoint points have not started; use the original plan to recover those points")
    for point in points
        selected=filter(r->r.point_id==point,comparisons)
        isempty(selected) || (length(selected)==length(names) && Set(Symbol(r.observable) for r in selected)==Set(names)) ||
            error("partially written checkpoint comparison")
        # Completed candidates were already consumed through the native workflow.
        # Do not reconstruct their cache identity from a two-column sigma table.
        count(r->r.point_id==point && r.configuration=="candidate",observables)==1 ||
            error("native candidate checkpoint is incomplete; rerun that candidate through its full workflow")
    end
    evidence=(directory=abspath(directory),manifest_sha256=bytes2hex(sha256(read(joinpath(directory,"manifest.json")))),
        files=[(name=name,sha256=bytes2hex(sha256(read(joinpath(directory,name)))))
            for name in ("background.csv","rates.csv","observables.csv","comparisons.csv")])
    return (;manifest,candidate,rates,observables,comparisons,backgrounds,points,evidence)
end

function resume_joint_validation(directory,outdir)
    old=verify_joint_checkpoint(directory)
    inputs=joinpath(directory,"inputs")
    rows=filter(r->String(r.point_id) in old.points,collect(CSV.File(joinpath(inputs,"review_points.csv"))))
    length(rows)==length(old.points) || error("checkpoint points differ from the input snapshot")
    completed=Set(String(r.point_id) for r in old.comparisons)
    params=merge(old.manifest["parameters"],Dict("purpose"=>"resume interrupted joint qualification from verified checkpoints",
        "recovered_from"=>old.evidence,"completed_points_reused"=>length(completed),
        "rates_reused"=>length(old.rates),"timing"=>"mixed retained and recovery diagnostics; not production timing"))
    start_run(outdir,inputs,@__FILE__,params;extra_sources=[joinpath(@__DIR__,f) for f in
        ("run_transport_joint_validation.jl","run_transport_background_tail_validation.jl",
         "run_transport_background_validation.jl","TerminalTransportQuadrature.jl","ThermalWeightedSampling.jl",
         "ThermalAverageCandidates.jl","ThresholdCacheCandidates.jl","StreamedThermalValidation.jl")])
    rates=old.rates;observables=old.observables;comparisons=old.comparisons;failures=NamedTuple[]
    for (name,data) in (("rates.csv",rates),("observables.csv",observables),
        ("comparisons.csv",comparisons),("background.csv",old.backgrounds))
        write_rows(joinpath(outdir,name),data)
    end
    new_rates=0
    for row in rows
        point=String(row.point_id)
        if point in completed
            println("REUSED COMPLETE ",point);flush(stdout);continue
        end
        println("RESUME REFERENCE ",point);flush(stdout)
        try
            config=joint_config("reference")
            background=joint_background(row,config)
            qp,tp=background.state.quark_params,background.state.thermo_params
            stored=only(r for r in old.backgrounds if r.point_id==point && r.configuration=="reference")
            actual=(m_u=qp.m.u,m_d=qp.m.d,m_s=qp.m.s,Phi=tp.Φ,Phibar=tp.Φbar,A_u=qp.A.u,A_d=qp.A.d,A_s=qp.A.s)
            all(k->isapprox(getproperty(actual,k),getproperty(stored,k);rtol=1e-12,atol=1e-12),keys(actual)) ||
                error("reconstructed reference background differs from the checkpoint")
            values=Dict(Symbol(r.channel)=>Float64(r.rate) for r in rates if r.point_id==point && r.configuration=="reference")
            for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
                haskey(values,process) && continue
                Main.RelaxTime.PolarizationCache.reset_cache!();GC.gc()
                id=point*"_"*String(process)*"_reference"
                measured=@timed independent_joint_rate(process,background.state,config;
                    history_path=joinpath(outdir,"histories",id*".csv"),table_path=joinpath(outdir,"tables",id*".csv"))
                value=measured.value;values[process]=value.rate;new_rates+=1
                push!(rates,(point_id=point,configuration="reference",channel=String(process),rate=value.rate,
                    estimated_or_reference_l1=value.reference_l1,reference_order=value.order,
                    evaluations=value.evaluations,seconds=measured.time))
                write_rows(joinpath(outdir,"rates.csv"),rates)
                println("RATE reference ",process," ",value.rate," seconds=",measured.time);flush(stdout)
            end
            measured=@timed joint_transport(row,background,config;tau_kwargs=(existing_rates=NamedTuple(values),))
            reference=joint_observables(measured.value.result)
            push!(observables,merge((point_id=point,configuration="reference",
                density_relative_error=measured.value.density_error,transport_seconds=measured.time),reference))
            write_rows(joinpath(outdir,"observables.csv"),observables)
            candidate=only(r for r in observables if r.point_id==point && r.configuration=="candidate")
            for name in keys(reference)
                c,r=getproperty(candidate,name),getproperty(reference,name)
                push!(comparisons,(point_id=point,observable=String(name),candidate=c,reference=r,
                    relative_change=c/r-1,passed=abs(c/r-1)<=0.01))
            end
            write_rows(joinpath(outdir,"comparisons.csv"),comparisons)
            println("POINT ",point," max_relative=",maximum(abs(r.relative_change) for r in comparisons if r.point_id==point));flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=point,error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
            println("FAILED ",sprint(showerror,err));flush(stdout)
        end
        GC.gc()
    end
    passed=isempty(failures) && length(comparisons)==9length(rows) && all(r.passed for r in comparisons)
    finish_run(outdir;points=length(rows),comparisons=length(comparisons),failures=length(failures),
        completed_points_reused=length(completed),rates_reused=length(old.rates)-new_rates,
        new_reference_rates=new_rates,joint_target_passed=passed)
    passed || error("resumed joint qualification did not pass")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)>=2 && iseven(length(ARGS)) || error("Usage: resume_transport_joint_validation.jl SOURCE_RUN NEW_OUTPUT [SOURCE_RUN NEW_OUTPUT ...]")
    for i in 1:2:length(ARGS)
        resume_joint_validation(abspath(ARGS[i]),abspath(ARGS[i+1]))
    end
end
