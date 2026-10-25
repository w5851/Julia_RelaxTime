#!/usr/bin/env julia
# Summarize converged outer-grid references through the actual tau assembler.
using CSV, JSON, SHA
include(joinpath(@__DIR__,"ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon

function summarize_outer(inputdir,outdir;repair_dir=nothing)
    basedir=dirname(inputdir)
    sources=NamedTuple[];latest=Dict{Tuple{String,String,String},Any}()
    directories=[joinpath(basedir,"outer_thermal_qualification_s$(shard)_v1") for shard in 0:1]
    repair_dir===nothing || push!(directories,repair_dir)
    unresolved=Set{Tuple{String,String}}()
    for dir in directories
        status=JSON.parsefile(joinpath(dir,"status.json"))
        status["complete"] && (repair_dir!==nothing || status["failures"]==0) ||
            error("outer reference shard is incomplete or failed: $dir")
        sha256(read(joinpath(dir,"inputs","review_points.csv")))==
            sha256(read(joinpath(inputdir,"review_points.csv"))) || error("outer source input differs")
        status["failures"]==0 || union!(unresolved,
            ((String(r.point_id),String(r.channel)) for r in CSV.File(joinpath(dir,"failures.csv"))))
        path=joinpath(dir,"rates.csv")
        push!(sources,(path=abspath(path),sha256=bytes2hex(sha256(read(path)))))
        for row in CSV.File(path)
            key=(String(row.point_id),String(row.channel),String(row.grid))
            if dir==repair_dir || !haskey(latest,key) || row.reference_order>latest[key].reference_order
                latest[key]=row
            end
        end
    end
    for (point,channel) in unresolved
        all(haskey(latest,(point,channel,g)) && latest[(point,channel,g)].reference_converged
            for g in ("base","cost_candidate","joint")) || error("unresolved outer failure: $point $channel")
    end
    points=Dict(String(r.point_id)=>r for r in CSV.File(joinpath(inputdir,"review_points.csv")))
    selected=sort!(unique(first(k) for k in keys(latest)))
    start_run(outdir,inputdir,@__FILE__,(purpose="full six-species tau response to outer quadrature and cache error",
        source_rates=sources,resolved_failures=collect(unresolved),points=selected,relative_target=0.01,
        reference="independently validated sigma table on p192/c28/phi54, density p256",
        tau_densities="same retained background; semi-infinite p128 scale5 and polar32",
        scope="fixed masses, A, couplings and t24; background and terminal transport are separate gates"))
    species=(:u,:d,:s,:ubar,:dbar,:sbar)
    channels=Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
    p,pw=AS._cached_semi_infinite_momentum_grid(128,5.0)
    rate_rows=NamedTuple[];tau_rows=NamedTuple[];density_rows=NamedTuple[]
    for point in selected
        state=retained_state(points[point]);qp,tp=state.quark_params,state.thermo_params
        densities=NamedTuple{species}(map(species) do s
            AS.number_density(s,AS.get_mass(s,qp),AS.get_mu(s,qp),tp.T,tp.Φ,tp.Φbar,tp.ξ;
                p_grid=p,p_w=pw,angle_nodes=32)
        end)
        append!(density_rows,[(point_id=point,species=String(s),density=getproperty(densities,s)) for s in species])
        for channel in channels
            ref=latest[(point,String(channel),"joint")]
            ref.reference_converged || error("reference not converged: $point $channel")
            for grid in ("base","cost_candidate","joint")
                row=latest[(point,String(channel),grid)]
                row.reference_converged || error("unconverged final reference")
                push!(rate_rows,(point_id=point,channel=String(channel),grid=grid,
                    reference_order=row.reference_order,candidate_rate=row.candidate_rate,
                    reference_rate=row.reference_rate,joint_reference=ref.reference_rate,
                    outer_relative_change=row.reference_rate/ref.reference_rate-1,
                    cache_and_outer_relative_change=row.candidate_rate/ref.reference_rate-1,
                    cache_validated=row.cache_validated,cache_l1_relative=row.cache_l1_relative,
                    reference_l1_relative=row.reference_l1_relative))
            end
        end
        reference_rates=NamedTuple{channels}(Tuple(Float64(latest[(point,String(c),"joint")].reference_rate) for c in channels))
        reference_inverse=Main.RelaxTime.RelaxationTime.relaxation_rates(densities,reference_rates)
        for grid in ("base","cost_candidate","joint"),method in ("reference_table","adaptive_cache")
            field=method=="reference_table" ? :reference_rate : :candidate_rate
            rates=NamedTuple{channels}(Tuple(Float64(getproperty(latest[(point,String(c),grid)],field)) for c in channels))
            inverse=Main.RelaxTime.RelaxationTime.relaxation_rates(densities,rates)
            for s in species
                tau=1/getproperty(inverse,s);reference_tau=1/getproperty(reference_inverse,s)
                push!(tau_rows,(point_id=point,grid=grid,method=method,species=String(s),
                    tau=tau,reference_tau=reference_tau,relative_change=tau/reference_tau-1))
            end
        end
    end
    write_rows(joinpath(outdir,"rates.csv"),rate_rows)
    write_rows(joinpath(outdir,"tau.csv"),tau_rows)
    write_rows(joinpath(outdir,"densities.csv"),density_rows)
    summary=[(grid=g,method=m,max_relative_change=maximum(abs(r.relative_change) for r in tau_rows if r.grid==g && r.method==m))
        for g in ("base","cost_candidate","joint") for m in ("reference_table","adaptive_cache")]
    write_json(joinpath(outdir,"summary.json"),(tau_summary=summary,
        cache_qualification=all(r.cache_validated for r in rate_rows),
        max_cache_l1_relative=maximum(r.cache_l1_relative for r in rate_rows),
        max_reference_l1_relative=maximum(r.reference_l1_relative for r in rate_rows)))
    finish_run(outdir;points=length(selected),channels=length(selected)*length(channels),
        tau_comparisons=length(tau_rows),cache_qualification=all(r.cache_validated for r in rate_rows),
        tau_target_passed=all(r.max_relative_change<=0.01 for r in summary))
    foreach(println,summary)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=3 || error("Usage: summarize_thermal_outer_validation.jl INPUT_DIR NEW_OUTPUT_DIR [REPAIR_DIR]")
    summarize_outer(abspath(ARGS[1]),abspath(ARGS[2]);repair_dir=length(ARGS)==3 ? abspath(ARGS[3]) : nothing)
end
