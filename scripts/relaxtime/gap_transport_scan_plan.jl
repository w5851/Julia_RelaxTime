module GapTransportScanPlan

function build_scan_plan(opts)
    T_values = collect(range(opts.tmin_mev; stop=opts.tmax_mev, step=opts.tstep_mev))
    muB_values = unique(sort(Float64.(collect(range(opts.mubmin_mev; stop=opts.mubmax_mev, step=opts.mubstep_mev)))))
    xi_continuity_mode = length(T_values) == 1 && length(muB_values) == 1 && length(opts.xi_values) > 1
    total = length(opts.xi_values) * length(T_values) * length(muB_values)
    return (
        T_values=T_values,
        muB_values=muB_values,
        xi_continuity_mode=xi_continuity_mode,
        total=total,
    )
end

@inline function _maybe_collect_gc(opts, done::Int)
    if opts.gc_every_n > 0 && done % opts.gc_every_n == 0
        GC.gc()
    end
    return nothing
end

@inline function _maybe_print_progress(done::Int, total::Int, T_mev::Float64, muB_mev::Float64, xi::Float64)
    if done % 10 == 0
        println("progress: $(done)/$(total) (last muB=$(muB_mev) MeV, T=$(T_mev) MeV, xi=$(xi))")
    end
    return nothing
end

function execute_scan_plan!(point_runner!::Function, opts, plan, existing)
    stats_success = 0
    stats_error = 0
    stats_skipped = 0
    done = 0

    if plan.xi_continuity_mode
        T_mev = plan.T_values[1]
        muB_mev = plan.muB_values[1]
        previous_solution = nothing
        previous_phase = :unknown

        for xi in opts.xi_values
            done += 1
            key = (T_mev, muB_mev, xi)
            if opts.resume && (key in existing)
                stats_skipped += 1
                continue
            end

            point_result = point_runner!(T_mev, muB_mev, xi, previous_solution, previous_phase)
            previous_solution = point_result.next_solution
            previous_phase = point_result.next_phase
            if point_result.success
                stats_success += 1
            else
                stats_error += 1
            end

            _maybe_collect_gc(opts, done)
            _maybe_print_progress(done, plan.total, T_mev, muB_mev, xi)
        end
    else
        for xi in opts.xi_values
            for muB_mev in plan.muB_values
                previous_solution = nothing
                previous_phase = :unknown
                for T_mev in plan.T_values
                    done += 1
                    key = (T_mev, muB_mev, xi)
                    if opts.resume && (key in existing)
                        stats_skipped += 1
                        continue
                    end

                    point_result = point_runner!(T_mev, muB_mev, xi, previous_solution, previous_phase)
                    previous_solution = point_result.next_solution
                    previous_phase = point_result.next_phase
                    if point_result.success
                        stats_success += 1
                    else
                        stats_error += 1
                    end

                    _maybe_collect_gc(opts, done)
                    _maybe_print_progress(done, plan.total, T_mev, muB_mev, xi)
                end
            end
        end
    end

    return (
        success=stats_success,
        error=stats_error,
        skipped=stats_skipped,
        done=done,
        total=plan.total,
    )
end

export build_scan_plan
export execute_scan_plan!

end # module GapTransportScanPlan
