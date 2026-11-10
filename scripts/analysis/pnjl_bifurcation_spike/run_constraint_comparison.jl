module PALCConstraintComparison

using JSON
using LinearAlgebra
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end
using Main.Models
const HBARC = Main.Constants_PNJL.ħc_MeV_fm

function row_record(model, path, value, result)
    result === nothing && return (; coordinate=value, converged=false)
    options = Models._path_options(model, path)
    residual = Models._path_residual(model, path, options)
    full_residual = residual(result.solution, value)
    return (; coordinate=value, converged=result.converged,
        residual=norm(full_residual), constraint_error=norm(full_residual[6:end]),
        solution=collect(result.solution), mu_vec=collect(result.mu_vec),
        pressure=result.pressure, rho_norm=result.rho_norm, entropy=result.entropy,
        energy=result.energy, masses=collect(result.masses))
end

function finite_json(value)
    value isa AbstractFloat && return isfinite(value) ? value : nothing
    value isa NamedTuple && return Dict(string(k)=>finite_json(v) for (k,v) in pairs(value))
    value isa AbstractDict && return Dict(string(k)=>finite_json(v) for (k,v) in pairs(value))
    value isa AbstractArray && return map(finite_json, value)
    value isa Symbol && return string(value)
    return value
end

function write_results(destination, records)
    mkpath(dirname(destination))
    open(destination, "w") do io
        JSON.print(io, finite_json((; source_sha=strip(read(`git rev-parse HEAD`, String)),
            actions_run_id=get(ENV, "GITHUB_RUN_ID", ""), p_num=24, t_num=8,
            residual_tolerance=1e-6, records)), 2)
        println(io)
    end
end

function run_comparison(destination; repetitions=2, scope="all")
    LinearAlgebra.BLAS.set_num_threads(1)
    model = Models.create_model(:PNJL)
    cases = [
        (:entropy, Models.ConstraintTPath(Models.FixedEntropy(0.5), collect(100.0:1.0:120.0)./HBARC; p_num=24, t_num=8)),
        (:sigma, Models.ConstraintTPath(Models.FixedSigma(12.0759), collect(100.0:2.0:140.0)./HBARC; p_num=24, t_num=8)),
        (:bqs, Models.ConstraintTPath(Models.FixedMuBConservedCharges(900.0/HBARC), collect(100.0:3.0:160.0)./HBARC; p_num=24, t_num=8)),
    ]
    records = NamedTuple[]
    for (name,path) in cases
        scope in ("all", string(name)) || continue
        for method in (:palc,:legacy), repetition in 1:repetitions
            println("START ", name, "/", method, " round ", repetition)
            flush(stdout)
            work = Models.PALCWorkTelemetry()
            started = time()
            result = Models.solve_path(model, path; continuation=method,
                iterations=300, (method===:palc ? (;path_work=work) : (;))...)
            selected = Models.selected_path_results(model, result; iterations=300)
            seconds = time()-started
            rows = [row_record(model,path,value,point) for (value,point) in zip(path.T_values,selected)]
            push!(records, (; name, method, repetition, seconds,
                config=Models.to_namedtuple(path), rows,
                work=Models.palc_work_snapshot(work),
                diagnostics=Models.to_namedtuple(result.diagnostics),
                candidate_counts=[s.candidate_branch_count for s in result.selections]))
            write_results(destination,records)
            println("DONE ", name, "/", method, " ", seconds, " s")
            flush(stdout)
        end
    end
    if scope in ("all","single_seed")
        for (name,grid) in ((:tmu_full,collect(0.0:15.0:450.0)),
                            (:tmu_coexistence,collect(280.0:2.0:320.0)))
            path = Models.FixedMuPath(120.0/HBARC, grid./HBARC; p_num=24,t_num=8)
            for method in (:palc,:single_seed), repetition in 1:repetitions
                println("START ", name, "/", method, " round ", repetition)
                flush(stdout)
                value = first(path.mu_values)
                seed = Models.get_seed(Models.DefaultSeed(), [path.T_fm,value], Models.FixedMu())
                anchor = method===:palc ? Models.SeedPoolAnchor() : Models.SingleSeedAnchor(value,seed)
                work = Models.PALCWorkTelemetry()
                started = time()
                result = Models.solve_path(model,path;
                    continuation_strategy=Models.PALCContinuation(),anchor_strategy=anchor,
                    path_work=work,iterations=300)
                if method === :single_seed
                    work.anchor_solves == 1 || error("single-seed experiment must attempt exactly one anchor")
                    work.fallback_solves == 0 || error("single-seed experiment must not use catalog fallback")
                end
                selected = Models.selected_path_results(model,result; iterations=300)
                seconds = time()-started
                rows = [row_record(model,path,value,point) for (value,point) in zip(path.mu_values,selected)]
                push!(records, (; name,method,repetition,seconds,config=Models.to_namedtuple(path),rows,
                    work=Models.palc_work_snapshot(work),
                    diagnostics=Models.to_namedtuple(result.diagnostics),
                    candidate_counts=[s.candidate_branch_count for s in result.selections]))
                write_results(destination,records)
                println("DONE ", name, "/", method, " ", seconds, " s")
                flush(stdout)
            end
        end
    end
    return records
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PALCConstraintComparison.run_comparison(ARGS[1];
        repetitions=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 2,
        scope=length(ARGS)>=3 ? ARGS[3] : "all")
end
