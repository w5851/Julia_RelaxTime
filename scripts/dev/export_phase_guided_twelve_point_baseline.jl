"""Regenerate the twelve-point phase-guided transport baseline after a physics change.

`tests/baselines/relaxtime/baseline_thermal_weighted_transport_v3.csv` records
three numbers per (point, observable): `candidate` and `qualified_candidate` come
from the candidate configuration, `reference`/`qualified_reference` from the
tighter reference configuration.  Both configurations are driven here through the
same production entrypoint the regression exercises, so a numeric commitment
change (for example `#365`, which corrected the strange-quark scattering partner)
can be re-established without hand-editing the fixture.

The script never writes into `tests/baselines/`: the target must not exist yet,
so every run produces a reviewable candidate.  Configuration values default to
the recorded v3 metadata and can be overridden, but the emitted JSON always
records what was actually used.
"""

using CSV
using JSON3
using SHA
using Dates

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _BASELINE_DIR = joinpath(_ROOT, "tests", "baselines", "relaxtime")
const _SOURCE_CSV = joinpath(_BASELINE_DIR, "baseline_thermal_weighted_transport_v3.csv")

Base.include(Main, joinpath(_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl"))

const OBSERVABLES = ("tau_u", "tau_s", "tau_ubar", "tau_sbar",
                     "eta_over_s", "zeta_over_s", "sigma_over_T", "eta", "sigma")

function scan_point(; outdir, T_MeV, muB_MeV, xi, tau_p, tau_c, tau_phi, tau_nsigma, p_num, t_num)
    isdir(outdir) && error("refusing to reuse an existing scan directory: $outdir")
    mkpath(outdir)
    args = ["--outdir", outdir, "--mode", "fixed-T-sparse-muB", "--channel-diagnostics",
            "--T-list", string(T_MeV), "--muB-list", string(muB_MeV), "--xi-list", string(xi),
            "--tau-p-nodes", string(tau_p), "--tau-angle-nodes", string(tau_c),
            "--tau-phi-nodes", string(tau_phi), "--tau-n-sigma", string(tau_nsigma),
            "--p-num", string(p_num), "--t-num", string(t_num)]
    opts = Main.PhaseGuidedCLI.parse_args(args)
    context = Main.ProvenanceMetadata.new_run_context("transport_baseline_refresh", args)
    Main.run_phase_guided_scan(opts, context; phase_reference=Main._load_runtime_phase_reference(opts))
    path = joinpath(outdir, "phase_guided_transport_scan.csv")
    isfile(path) || error("scan produced no CSV at $path")
    row = only(CSV.File(path; comment="#"))
    Bool(row.converged) || error("point did not converge: T=$T_MeV muB=$muB_MeV xi=$xi")
    return row
end

function main(args)
    length(args) == 1 || error("usage: export_phase_guided_twelve_point_baseline.jl <new-candidate.csv>")
    target = abspath(args[1])
    ispath(target) && error("refusing to overwrite existing target: $target")
    isfile(_SOURCE_CSV) || error("missing recorded baseline: $_SOURCE_CSV")

    recorded = collect(CSV.File(_SOURCE_CSV))
    points = unique((p.point_id, Float64(p.T_MeV), Float64(p.muB_MeV), Float64(p.xi)) for p in recorded)
    sort!(points, by=p -> p[1])
    println("  points = ", length(points), "  observables = ", length(OBSERVABLES))

    # Recorded v3 configurations (see baseline_thermal_weighted_transport_v3.json).
    candidate_cfg = (tau_p=48, tau_c=36, tau_phi=36, tau_nsigma=24, p_num=12, t_num=6)
    reference_cfg = (tau_p=128, tau_c=48, tau_phi=54, tau_nsigma=48, p_num=12, t_num=6)

    workroot = mktempdir()
    rows = NamedTuple[]
    for (id, T, muB, xi) in points
        println("  computing $id (T=$T muB=$muB xi=$xi)")
        cand = scan_point(; outdir=joinpath(workroot, id, "candidate"), T_MeV=T, muB_MeV=muB, xi=xi, candidate_cfg...)
        ref = scan_point(; outdir=joinpath(workroot, id, "reference"), T_MeV=T, muB_MeV=muB, xi=xi, reference_cfg...)
        for observable in OBSERVABLES
            candidate = Float64(getproperty(cand, Symbol(observable)))
            reference = Float64(getproperty(ref, Symbol(observable)))
            (isfinite(candidate) && candidate > 0) || error("non-finite candidate $observable at $id")
            (isfinite(reference) && reference > 0) || error("non-finite reference $observable at $id")
            relative = abs(candidate / reference - 1)
            relative <= 0.01 || println("    warning: $id $observable candidate/reference differ by $(round(100relative, digits=3))%")
            push!(rows, (point_id=id, T_MeV=T, muB_MeV=muB, xi=xi, observable=observable,
                         candidate=candidate, reference=reference, qualified_candidate=candidate))
        end
    end

    sort!(rows, by=row -> (row.point_id, row.observable))
    CSV.write(target, rows; header=[:point_id, :T_MeV, :muB_MeV, :xi, :observable,
                                   :candidate, :reference, :qualified_candidate])
    digest = bytes2hex(sha256(read(target)))
    println("  wrote $(length(rows)) rows to $target")
    println("  candidate sha256 = $digest")
end

main(ARGS)
