#!/usr/bin/env julia

module TransportBaselineExport

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const POINTS = (
    (T=0.50, mu=0.00, xi=0.0), (T=0.50, mu=1.50, xi=0.0),
    (T=0.66, mu=0.50, xi=0.0), (T=0.66, mu=1.00, xi=0.0),
    (T=0.75, mu=0.00, xi=0.0), (T=0.75, mu=0.75, xi=0.0),
    (T=0.90, mu=0.00, xi=0.0), (T=0.90, mu=0.15, xi=0.0),
    (T=0.90, mu=0.30, xi=0.0), (T=1.05, mu=0.00, xi=0.0),
    (T=0.90, mu=0.00, xi=0.2), (T=0.90, mu=0.15, xi=0.2),
    (T=0.90, mu=0.00, xi=-0.2),
)

function parse_args(args; io=stdout)
    output = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--output", "--backend")
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            if arg == "--output"
                output === nothing || throw(ArgumentError("duplicate --output"))
                output = args[i]
            else
                lowercase(args[i]) == "models" || throw(ArgumentError("backend must be models"))
            end
        elseif arg in ("-h", "--help")
            println(io, "Usage: julia --project=. scripts/dev/export_transport_fixedpoint_baseline.jl --output <new-candidate.csv> [--backend models]")
            println(io, "Existing files are never overwritten; no numerical model is loaded for help or rejected paths.")
            return nothing
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end
    output === nothing && throw(ArgumentError("--output is required; specify a new candidate CSV"))
    isempty(strip(output)) && throw(ArgumentError("--output must not be empty"))
    return abspath(output)
end

function check_new_output(output)
    endswith(lowercase(output), ".csv") || throw(ArgumentError("--output must name a CSV file"))
    (ispath(output) || islink(output)) && throw(ArgumentError("output already exists: $output"))
end

function build_compute_point()
    models = Main.Models
    workflow = models.transport_workflow_module()
    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    solver = models.NLsolveGapSolver(method=:trust_region, jacobian=:finite, xtol=1e-10, ftol=1e-10)
    return pt -> workflow.solve_gap_and_transport(
        pt.T, pt.mu; xi=pt.xi, tau=tau, compute_tau=false, compute_bulk=true,
        p_num=8, t_num=4,
        transport_config=workflow.TransportIntegrationConfig(p_nodes=8, p_max=3.5),
        solver_backend=:models, models_solver=solver, models_residual_norm_max=1e-4,
        seed_state=workflow.HADRON_SEED_5,
    )
end

function load_compute_point()
    if !isdefined(Main, :Models)
        Base.include(Main, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
    end
    return Base.invokelatest(build_compute_point)
end

function validate_candidate(path)
    lines = readlines(path)
    length(lines) == length(POINTS) + 1 || error("candidate row count mismatch")
    first(lines) == "T,mu,xi,eta,sigma,zeta" || error("candidate header mismatch")
    for (line, point) in zip(lines[2:end], POINTS)
        columns = split(line, ',')
        length(columns) == 6 || error("candidate column count mismatch")
        values = parse.(Float64, columns)
        all(isfinite, values) || error("candidate contains nonfinite values")
        values[1:3] == [point.T, point.mu, point.xi] || error("candidate point order mismatch")
    end
    return nothing
end

function export_candidate(output, compute_point)
    check_new_output(output)
    mkpath(dirname(output))
    temporary, stream = mktemp(dirname(output))
    try
        println(stream, "T,mu,xi,eta,sigma,zeta")
        for pt in POINTS
            result = try
                Base.invokelatest(compute_point, pt)
            catch err
                error("transport computation failed at $pt: $(sprint(showerror, err))")
            end
            result.equilibrium.converged === true || error("equilibrium did not converge at $pt")
            values = (result.transport.eta, result.transport.sigma, result.transport.zeta)
            all(isfinite, values) || error("nonfinite transport result at $pt")
            @printf(stream, "%.6f,%.6f,%.6f,%.16e,%.16e,%.16e\n", pt.T, pt.mu, pt.xi, values...)
        end
        close(stream)
        validate_candidate(temporary)
        hardlink(temporary, output) # Atomic, no-clobber publication on the same filesystem.
    finally
        isopen(stream) && close(stream)
        isfile(temporary) && rm(temporary)
    end
    return output
end

function main(args=collect(String.(ARGS)); io=stdout, compute_point=nothing)
    try
        output = parse_args(args; io)
        output === nothing && return 0
        check_new_output(output)
        compute_point === nothing && (compute_point = load_compute_point())
        export_candidate(output, compute_point)
        println(io, "candidate exported to: $output\nbackend = models\npoints = $(length(POINTS))")
        return 0
    catch err
        println(io, "[transport-baseline] FAILED: ", sprint(showerror, err))
        return 1
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(TransportBaselineExport.main())
end
