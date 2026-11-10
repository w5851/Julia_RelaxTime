const PRODUCTION_COMPARISON_LOAD_START = time()
include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))

module PALCProductionComparison
using Main.Models
using LinearAlgebra
using Statistics
import JSON

const HBARC = Main.Constants_PNJL.ħc_MeV_fm

function scan_contract_check(route, strategy, temperatures, anisotropies, grid, scan_file)
    bytes_before = read(scan_file)
    callbacks = Ref(0)
    common = (; T_values=temperatures, xi_values=anisotropies, continuation=strategy,
        output_path=scan_file, overwrite=false, resume=true, p_num=24, t_num=8,
        iterations=300, residual_norm_max=1e-6,
        progress_cb=(args...) -> (callbacks[] += 1))
    stats = route === :tmu ? Models.run_tmu_scan(; common..., mu_values=grid) :
        Models.run_trho_scan(; common..., rho_values=grid, reverse_rho=true,
            constraint_mode=route === :asym ? :fixed_asymmetric_rho : :fixed_rho,
            seed_policy=:candidates, evaluate_all_attempts=true)
    stats.skipped == length(temperatures)*length(anisotropies)*length(grid) || error("resume coverage changed")
    callbacks[] == 0 || error("resume recomputed complete rows")
    read(scan_file) == bytes_before || error("resume changed existing CSV")
    return (; skipped=stats.skipped, byte_preserving=true)
end

function point_record(coordinate, result, route, T, xi)
    result === nothing && return (; T, xi, coordinate, converged=false)
    # TmuScan retains its historical physical baryon-density field.
    rho = route === :tmu ? result.rho_norm / Main.Constants_PNJL.ρ0_inv_fm3 : result.rho_norm
    nB = rho * Main.Constants_PNJL.ρ0_inv_fm3
    chemical_term = route === :asym ?
        3*nB*(0.876*result.mu_vec[1]+result.mu_vec[2])/1.876 : 3*nB*mean(result.mu_vec)
    flavor_density = Models.model_rho(Models.create_model(:PNJL), result.x_state, result.mu_vec,
        T / HBARC; xi, p_num=24, t_num=8)
    density_error = abs(sum(flavor_density)/3 - nB)
    constraint_error = if route === :tmu
        maximum(abs.(result.mu_vec .- coordinate/HBARC))
    elseif route === :rho
        max(abs(rho-coordinate), maximum(abs.(result.mu_vec .- mean(result.mu_vec))))
    else
        max(abs(rho-coordinate), abs(flavor_density[1]/flavor_density[2]-0.876), abs(flavor_density[3]))
    end
    return (; T, xi, coordinate, converged=result.converged,
        solution=collect(result.solution), mu_vec=collect(result.mu_vec),
        pressure=result.pressure, rho_norm=rho, entropy=result.entropy, energy=result.energy,
        masses=collect(result.masses), residual=result.residual_norm, density_error, constraint_error,
        free_energy=-result.pressure + chemical_term)
end

function phase_pipeline_check(strategy, destination, grid)
    started = time()
    directory = joinpath(dirname(destination), "phase_$(strategy)")
    result = Models.run_phase_pipeline(:PNJL; mode=:production,
        T_grid=[120.0], rho_grid=grid, xi=0.0, output_dir=directory,
        continuation=strategy, seed_policy=:candidates, memoize_uniform=true,
        p_num=24, t_num=8, iterations=300, cep_max_refine_level=0,
        rho_geometry_convergence=false, promote_reference=false)
    result.config_snapshot["continuation"] == String(strategy) || error("phase strategy was not propagated")
    csv = joinpath(directory, "trho_scan.csv")
    lines = readlines(csv)
    length(lines)-1 == length(grid) || error("phase pipeline changed the bounded grid")
    if strategy === :palc
        all(line -> occursin("continuation=palc", line), lines[2:end]) || error("phase session did not use PALC")
    end
    curves = Models.load_curves_from_trho_csv(csv; xi=0.0, min_points=3)
    haskey(curves, 120.0) || error("phase pipeline did not produce the T120 curve")
    mu, rho = curves[120.0]
    shape = Models.detect_s_shape(mu, rho)
    maxwell = Models.maxwell_construction(mu, rho; spinodal_hint=shape)
    return (; seconds=time()-started, point_count=length(rho), strategy,
        has_s_shape=shape.has_s_shape, maxwell_converged=maxwell.converged,
        mu_transition=maxwell.mu_transition, rho_hadron=maxwell.rho_hadron,
        rho_quark=maxwell.rho_quark)
end

function curve_record(rows, route)
    isempty(rows) && return (;)
    ordered = sort(rows; by=row -> row.coordinate)
    good = filter(row -> row.converged, ordered)
    length(good) == length(rows) || return (; valid=false)
    values = [row.coordinate for row in ordered]
    mu = route === :asym ?
        [3 * (0.876*row.mu_vec[1] + row.mu_vec[2]) / 1.876 * HBARC for row in ordered] :
        [mean(row.mu_vec) * HBARC for row in ordered]
    negative = route === :tmu ? Int[] : findall(<(0), diff(mu) ./ diff(values))
    phase = if route === :rho
        shape = Models.detect_s_shape(mu, values)
        boundary = Models.maxwell_construction(mu, values; spinodal_hint=shape)
        (; has_s_shape=shape.has_s_shape, maxwell_converged=boundary.converged,
            mu_transition=boundary.mu_transition, rho_hadron=boundary.rho_hadron,
            rho_quark=boundary.rho_quark, area_residual=boundary.area_residual)
    else
        (;)
    end
    return (; valid=true, negative_slope_intervals=[(values[i],values[i+1]) for i in negative],
        chemical_potential_min=minimum(mu), chemical_potential_max=maximum(mu), phase)
end

function run_case(route::Symbol, strategy::Symbol, destination::String;
    quick::Bool=false, repetitions::Int=2)
    LinearAlgebra.BLAS.set_num_threads(1)
    # Every strategy sees identical complete isotherms. T120 crosses the first-order interval.
    temperatures = quick ? [120.0,170.0] : [100.0,120.0,140.0,170.0,200.0]
    anisotropies = quick ? [0.0] : [0.0,0.4]
    grid = route === :tmu ? (quick ? collect(260.0:20.0:380.0) : collect(0.0:15.0:450.0)) :
        (quick ? collect(0.6:0.4:3.0) : collect(0.2:0.15:3.5))
    repetitions > 0 || error("repetitions must be positive")
    rounds = NamedTuple[]
    load_seconds = time() - Main.PRODUCTION_COMPARISON_LOAD_START
    for repetition in 1:repetitions
        rows = NamedTuple[]
        path_details = NamedTuple[]
        work = Models.SolverWorkTelemetry()
        path_work = Models.PALCWorkTelemetry()
        scan_file = joinpath(dirname(destination), "$(route)_$(strategy)_$(repetition).csv")
        common = (; T_values=temperatures, xi_values=anisotropies, continuation=strategy,
            output_path=scan_file, overwrite=true, resume=false, p_num=24, t_num=8,
            iterations=300, residual_norm_max=1e-6, work_telemetry=work,
            (strategy === :legacy ? (;) : (; path_work))...,
            path_diagnostics_cb=(coordinate, result) -> push!(path_details,
                (; coordinate, diagnostics=Models.to_namedtuple(result.diagnostics),
                    work=isempty(result.branches) ? (;) : get(first(result.branches).diagnostics, :work, (;)),
                    branch_count=length(result.branches),
                    candidate_counts=[selection.candidate_branch_count for selection in result.selections])))
        started = time()
        stats = if route === :tmu
            Models.run_tmu_scan(; common..., mu_values=grid,
                progress_cb=(point,result) -> push!(rows, point_record(point.mu,result,route,point.T,point.xi)))
        else
            Models.run_trho_scan(; common..., rho_values=grid, reverse_rho=true,
                constraint_mode=route === :asym ? :fixed_asymmetric_rho : :fixed_rho,
                # Match the prior production equilibrium's full-candidate selection.
                seed_policy=:candidates, evaluate_all_attempts=true,
                progress_cb=(point,result) -> push!(rows, point_record(point.rho,result,route,point.T,point.xi)))
        end
        elapsed = time() - started
        expected_order = route === :tmu ? grid : reverse(grid)
        for xi in anisotropies, T in temperatures
            [row.coordinate for row in rows if row.T==T && row.xi==xi] == expected_order ||
                error("callback order differs from the requested grid")
        end
        contract = scan_contract_check(route, strategy, temperatures, anisotropies, grid, scan_file)
        sort!(rows; by=row -> (row.xi,row.T,row.coordinate))
        length(rows) == length(temperatures)*length(anisotropies)*length(grid) ||
            error("production callback omitted a result")
        curves = [(; T,xi,curve_record(filter(row -> row.T==T && row.xi==xi, rows),route)...)
            for xi in anisotropies, T in temperatures]
        push!(rounds, (; repetition, seconds=elapsed, stats, rows, contract, curves=vec(curves),
            paths=path_details, work=Models.solver_work_snapshot(work),
            single_pass_work=Models.palc_work_snapshot(path_work)))
        payload = (; schema="palc_production_comparison_v1", route, strategy,
            source_sha=readchomp(`git rev-parse HEAD`), julia_version=string(VERSION),
            bifurcationkit_version=string(Base.pkgversion(Models.BifurcationKit)),
            actions_run_id=get(ENV,"GITHUB_RUN_ID",""),
            load_seconds,
            config=(; temperatures,anisotropies,grid,p_num=24,t_num=8,residual_norm_max=1e-6),rounds)
        open(destination,"w") do io
            JSON.print(io,payload)
        end
        println("$(route) $(strategy) round $(repetition): $(round(elapsed;digits=3))s, $(stats.success)/$(stats.total)")
        flush(stdout)
    end
    if route === :rho
        phase_check = phase_pipeline_check(strategy, destination, grid)
        payload = JSON.parsefile(destination)
        payload["phase_pipeline_check"] = phase_check
        open(io -> JSON.print(io, payload), destination, "w")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 3 || error("usage: route strategy output [quick|medium] [repetitions]")
    mkpath(dirname(abspath(ARGS[3])))
    PALCProductionComparison.run_case(Symbol(ARGS[1]),Symbol(ARGS[2]),abspath(ARGS[3]);
        quick=length(ARGS)>=4 && ARGS[4]=="quick",
        repetitions=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 2)
end
