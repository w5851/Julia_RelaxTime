"""Regenerate the twelve-point phase-guided transport baseline after a physics change.

`tests/baselines/relaxtime/baseline_thermal_weighted_transport_v3.csv` records
three numbers per (point, observable): `candidate`/`qualified_candidate` from the
default sampling profile and `reference`/`qualified_reference` from the tighter
companion profile.  Both are driven here through the same production entrypoint
the regression exercises, so a numeric commitment change (for example `#365`,
which corrected the strange-quark scattering partner) can be re-established
without hand-editing the fixture.

The reference configuration lives in
`config/models/pnjl/transport_sampling/thermal_weighted_reference_v1.toml`, which
reproduces `configurations.reference` of the recorded v3 metadata.

The script never writes into `tests/baselines/`: the target must not exist yet,
so every run produces a reviewable candidate.
"""

using CSV
using JSON3
using SHA

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _BASELINE_DIR = joinpath(_ROOT, "tests", "baselines", "relaxtime")
const _SOURCE_CSV = joinpath(_BASELINE_DIR, "baseline_thermal_weighted_transport_v3.csv")
const _REFERENCE_PROFILE = joinpath(_ROOT, "config", "models", "pnjl", "transport_sampling",
                                    "thermal_weighted_reference_v1.toml")

Base.include(Main, joinpath(_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl"))

const OBSERVABLES = ("tau_u", "tau_s", "tau_ubar", "tau_sbar",
                     "eta_over_s", "zeta_over_s", "sigma_over_T", "eta", "sigma")

function scan_point(; outdir, T_MeV, muB_MeV, xi, profile=nothing)
    isdir(outdir) && error("refusing to reuse an existing scan directory: $outdir")
    mkpath(outdir)
    args = ["--outdir", outdir, "--mode", "fixed-T-sparse-muB", "--channel-diagnostics",
            "--T-list", string(T_MeV), "--muB-list", string(muB_MeV), "--xi-list", string(xi)]
    profile === nothing || append!(args, ["--sampling-profile", profile])
    opts = Main.PhaseGuidedCLI.parse_args(args)
    context = Main.ProvenanceMetadata.new_run_context("transport_baseline_refresh", args)
    Main.run_phase_guided_scan(opts, context; phase_reference=Main._load_runtime_phase_reference(opts))
    path = joinpath(outdir, "phase_guided_transport_scan.csv")
    isfile(path) || error("scan produced no CSV at $path")
    scanned = collect(CSV.File(path; comment="#"))
    # A point the orchestrator skipped leaves a header-only CSV and logs
    # `SCAN POINT FAILED ... skipped`, so the failure surfaces here rather than
    # in the scan.  The usual causes are an exhausted thermal-measure event
    # budget or a multi-threaded run changing how many events are consumed.
    if isempty(scanned)
        error("scan skipped the point (empty CSV) at T=$T_MeV muB=$muB_MeV xi=$xi; " *
              "check sampling.max_events and run with JULIA_NUM_THREADS=1")
    end
    row = only(scanned)
    Bool(row.converged) || error("point did not converge: T=$T_MeV muB=$muB_MeV xi=$xi")
    return row
end

function main(args)
    length(args) == 1 || error("usage: export_phase_guided_twelve_point_baseline.jl <new-candidate.csv>")
    target = abspath(args[1])
    ispath(target) && error("refusing to overwrite existing target: $target")
    isfile(_SOURCE_CSV) || error("missing recorded baseline: $_SOURCE_CSV")
    isfile(_REFERENCE_PROFILE) || error("missing reference profile: $_REFERENCE_PROFILE")

    # The thermal measure draws from one seeded event stream: extra threads change
    # how many events are consumed and can exhaust `sampling.max_events` on points
    # that converge single-threaded (observed on quote-remote at 8 threads).  The
    # regression and CI both run single-threaded, so the baseline must too.
    Threads.nthreads() == 1 || error("run the baseline export with JULIA_NUM_THREADS=1 " *
                                     "(current nthreads = $(Threads.nthreads()))")
    recorded = collect(CSV.File(_SOURCE_CSV))
    points = unique((p.point_id, Float64(p.T_MeV), Float64(p.muB_MeV), Float64(p.xi)) for p in recorded)
    sort!(points, by=p -> p[1])
    println("  points = ", length(points), "  observables = ", length(OBSERVABLES))
    flush(stdout)

    workroot = mktempdir()
    rows = NamedTuple[]
    for (id, T, muB, xi) in points
        println("  computing $id (T=$T muB=$muB xi=$xi)"); flush(stdout)
        cand = scan_point(; outdir=joinpath(workroot, id, "candidate"), T_MeV=T, muB_MeV=muB, xi=xi)
        ref = scan_point(; outdir=joinpath(workroot, id, "reference"), T_MeV=T, muB_MeV=muB, xi=xi,
                          profile=_REFERENCE_PROFILE)
        for observable in OBSERVABLES
            candidate = Float64(getproperty(cand, Symbol(observable)))
            reference = Float64(getproperty(ref, Symbol(observable)))
            (isfinite(candidate) && candidate > 0) || error("non-finite candidate $observable at $id")
            (isfinite(reference) && reference > 0) || error("non-finite reference $observable at $id")
            relative = abs(candidate / reference - 1)
            relative <= 0.01 || println("    warning: $id $observable candidate/reference differ by ",
                                        round(100relative, digits=3), "%")
            push!(rows, (point_id=id, T_MeV=T, muB_MeV=muB, xi=xi, observable=observable,
                         candidate=candidate, reference=reference, qualified_candidate=candidate))
        end
    end

    sort!(rows, by=row -> (row.point_id, row.observable))
    CSV.write(target, rows; header=[:point_id, :T_MeV, :muB_MeV, :xi, :observable,
                                   :candidate, :reference, :qualified_candidate])
    println("  wrote $(length(rows)) rows to $target")
    println("  candidate sha256 = ", bytes2hex(sha256(read(target))))
end

main(ARGS)
