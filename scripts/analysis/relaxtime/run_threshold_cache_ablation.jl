#!/usr/bin/env julia
include(joinpath(@__DIR__, "ThresholdUpgradeCommon.jl"))
using .ThresholdUpgradeCommon, CSV

function run_ablation(inputdir, outdir)
    start_run(outdir, inputdir, @__FILE__, (
        purpose="two-factor grid/subtraction ablation at retained backgrounds",
        xis=[0.37,0.36,0.38], channels=["uubar_to_uubar","uubar_to_ddbar"],
        rate_config=rate_config(), added_grid_levels=[96,192],
        subtraction_coefficient="fixed to stock A within each case; no refitting on dense grid",
        cache_validation="analysis-injected caches; full inputs and raw knots retained",
        reference="direct q-coordinate J at 128/256 nodes and t24/t48"))
    # Use the snapshot inputs throughout this run.
    inputdir = joinpath(outdir, "inputs")
    plan = collect(CSV.File(joinpath(inputdir, "channel_plan.csv")))
    rows = NamedTuple[]
    reference_rows = NamedTuple[]
    started = time()
    for xi in (0.37, 0.36, 0.38), process in (:uubar_to_uubar, :uubar_to_ddbar)
        println("CASE xi=", xi, " channel=", process); flush(stdout)
        row = retained_row(inputdir, xi)
        state = retained_state(row)
        qp,tp,kc = state.quark_params,state.thermo_params,state.K_coeffs
        id = replace("$(row.review_mode)_$(row.plot_panel)_$(row.plot_series)_xi$(row.xi)",
                     "."=>"p", "-"=>"m")
        job = only(j for j in plan if String(j.point_id)==id && Symbol(j.channel)==process)
        s0 = 4qp.m.u^2
        stock = AS.build_w0cdf_pchip_cache(process, qp, tp, kc; N=560,
            p_cutoff=Main.Constants_PNJL.Λ_inv_fm, n_sigma_points=24,
            threshold_subtraction=true, asym_window=0.6, asym_fit_min_points=8,
            asym_extra_points=10, sigma_cache_policy=:validated_anchored)
        original = Set(stock.s_vals)
        raw = Dict(s=>AS._exact_cached_sigma_value(stock,s,stock.sigma_vals[i])
                   for (i,s) in enumerate(stock.s_vals))
        sparse_off = raw_cache(process, stock.s_vals, [raw[s] for s in stock.s_vals], stock)
        # Same cumulative levels as the retained dense192 countercheck.
        for level in (96,192)
            xs = vcat(10.0 .^ range(log10(s0*1e-10),log10(0.6);length=level),
                      [0.6*(i/level)^2 for i in 1:level])
            for x in xs
                s = s0+x
                haskey(raw,s) || (raw[s] = CS.total_cross_section(process,s,qp,tp,kc;n_points=24))
            end
        end
        grid = sort!(collect(keys(raw)))
        values = [raw[s] for s in grid]
        dense_off = raw_cache(process, grid, values, stock)
        dense_on = raw_cache(process, grid, values, stock; subtract=true)
        write_rows(joinpath(outdir,"knots","xi$(xi)_$(process).csv"),
            [(s=s,x=s-s0,sigma=raw[s],stock_node=s in original) for s in grid])
        println("DIRECT J ", xi, " ", process); flush(stdout)
        j128 = direct_integral(process,state;order=128)
        j256 = direct_integral(process,state;order=256)
        j128t48 = direct_integral(process,state;order=128,t_nodes=48)
        push!(reference_rows,(xi=xi,channel=String(process),J128=j128,J256=j256,
            J128t48=j128t48,relative_order_change=j256/j128-1,
            relative_t_change=j128t48/j128-1))
        write_rows(joinpath(outdir,"direct_reference.csv"),reference_rows)
        for (method,cache,subtract) in (
            ("stock_on",stock,true),("stock_off",sparse_off,false),
            ("dense_on",dense_on,true),("dense_off",dense_off,false))
            println("RATE ",xi," ",process," ",method); flush(stdout)
            cfg = merge(rate_config(),(threshold_subtraction=subtract,require_cache_fingerprint=false))
            rate_timed = @timed AS.average_scattering_rate(process,qp,tp,kc;cfg...,cs_cache=cache)
            rate = rate_timed.value
            if method=="stock_on"
                isapprox(rate,Float64(job.source_rate);rtol=1e-9,atol=0.0) ||
                    error("Stock rate does not reproduce retained source at xi=$xi $process")
            end
            j12 = local_integral(cache,state;order=12)
            j24 = local_integral(cache,state;order=24)
            push!(rows,(xi=xi,channel=String(process),method=method,
                nodes=length(cache.s_vals),subtraction=subtract,template_A=stock.asym_A,
                rate=rate,source_rate=Float64(job.source_rate),
                source_contribution=Float64(job.source_contribution),
                source_total=Float64(job.source_total),
                J=j24,J_relative_error=j24/j256-1,
                cache_quadrature_change=j24/j12-1,
                rate_seconds=rate_timed.time,rate_bytes=rate_timed.bytes))
            write_rows(joinpath(outdir,"ablation.csv"),rows)
            write_rows(joinpath(outdir,"partial_tau.csv"),reconstruct_tau(rows))
            println("RESULT rate=",rate," Jerr=",j24/j256-1);flush(stdout)
        end
    end
    finish_run(outdir; cases=6, comparisons=length(rows), elapsed_seconds=time()-started)
    println("COMPLETE ",outdir)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("Usage: run_threshold_cache_ablation.jl INPUT_DIR NEW_OUTPUT_DIR")
    run_ablation(abspath(ARGS[1]),abspath(ARGS[2]))
end
