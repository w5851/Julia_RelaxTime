#!/usr/bin/env julia

isdefined(Main, :ChargedGBUContourScan) || include("run_charged_gbu_contour_scan.jl")

"""Finite-window sensitivity on four saved grid backgrounds; run on Actions only."""
module ChargedGBUQ0WindowAudit

using CSV, JSON3, Dates
using Main.Models

const S = Main.ChargedGBUContourScan
const B = S.SavedBackgrounds
const CASES = ((label="freezeout_140_425", T=140., muB=425.),
    (label="freezeout_80_725", T=80., muB=725.),
    (label="safe_160_0", T=160., muB=0.),
    (label="onset_145_375", T=145., muB=375.))

function variants()
    base = S.screening_settings(route="q0_lambda_reference")
    return [(label="window_1e-5", settings=base),
        (label="window_1e-4", settings=merge(base, (omega_lower_inv_fm=1e-4,))),
        (label="window_1e-6", settings=merge(base, (omega_lower_inv_fm=1e-6,))),
        (label="omega128", settings=merge(base, (omega_nodes=128,))),
        (label="profile128", settings=merge(base, (mesh=128, cut_nodes=64, tail_nodes=64,))),
        (label="q16", settings=merge(base, (q_nodes=16,))),
        (label="strict_zero_limit", settings=merge(base, (endpoint_policy="strict_zero_limit",)))]
end

function compare(base, other; atol=1e-12, rtol=1e-3)
    both = isfinite(base) && isfinite(other)
    absolute = both ? abs(other-base) : NaN
    relative = both ? absolute/max(abs(base),abs(other),1e-30) : NaN
    return (absolute_difference=absolute, relative_difference=relative,
        passed=both && absolute <= atol+rtol*max(abs(base),abs(other)))
end

function run(input_root, output)
    ispath(output) && throw(ArgumentError("refusing to overwrite finite-window audit"))
    manifests = sort([joinpath(d,f) for (d,_,fs) in walkdir(input_root) for f in fs if f == "manifest.json"])
    isempty(manifests) && throw(ArgumentError("no saved background manifests"))
    m = JSON3.read(read(first(manifests), String))
    Ts, mus = Float64.(m.T_grid), Float64.(m.muB_grid)
    saved = B.read_snapshot(input_root, Ts, mus)
    model = Models.create_model(:PNJL)
    rows = NamedTuple[]; ratios = NamedTuple[]; comparisons = NamedTuple[]
    records = Dict{Tuple{String,String,Symbol},Any}(); details = Dict{String,Any}(); inputs = NamedTuple[]
    mkpath(output)
    for case in CASES
        i, j = findfirst(==(case.T), Ts), findfirst(==(case.muB), mus)
        i !== nothing && j !== nothing || throw(ArgumentError("audit case missing from saved grid"))
        old = saved.points[(i,j)]
        old.background !== nothing || throw(ArgumentError("audit case lacks a saved background"))
        bg = B.restore_seed(model, Float64.(old.background.seed), case.T, case.muB, Float64(old.background.residual))
        push!(inputs, (label=case.label, bg=bg, seed=old.background.seed, source_scan_identity=old.scan_identity))
        for v in variants()
            channels = Dict{Symbol,Any}()
            for ch in S.CHANNELS
                result = S._channel_record(bg, ch, v.settings; route="q0_lambda_reference")
                records[(case.label,v.label,ch)] = result; channels[ch] = result
                push!(rows, (case=case.label, T_MeV=case.T, muB_MeV=case.muB,
                    variant=v.label, channel=String(ch), status=result.status,
                    density_inv_fm3=result.density_inv_fm3, warning_shells=result.endpoint_warning_shells,
                    static_imaginary_max=result.static_imaginary_max, reason=result.reason, elapsed_s=result.elapsed_s))
            end
            ratio = S._ratios(S.CHANNELS, channels)
            push!(ratios, merge((case=case.label, variant=v.label), ratio))
            details[case.label*"/"*v.label] = (settings=v.settings, channels=channels, ratios=ratio)
            println(case.label, " ", v.label, " plus=", ratio.Kplus_over_pi_plus,
                " statuses=", join([channels[ch].status for ch in S.CHANNELS], ",")); flush(stdout)
        end
    end
    # The window/quadrature gate tests the newly chosen prescription. Profile
    # and q refinement are reported separately as screening-resolution evidence.
    window_passed = true
    for case in CASES, v in variants()[2:end], ch in S.CHANNELS
        base = records[(case.label,"window_1e-5",ch)]; other = records[(case.label,v.label,ch)]
        c = compare(base.density_inv_fm3, other.density_inv_fm3)
        required = v.label in ("window_1e-4", "window_1e-6", "omega128")
        retained_failure = case.label == "onset_145_375" && ch in (:K_plus,:K_minus) &&
            !base.passed && !other.passed && occursin("onset",base.reason) && occursin("onset",other.reason)
        passed = c.passed || retained_failure
        required && (window_passed &= passed)
        push!(comparisons, merge((case=case.label, variant=v.label, channel=String(ch),
            required_window_check=required, retained_failure=retained_failure), c, (passed=passed,)))
    end
    for case in CASES, label in ("window_1e-4", "window_1e-6", "omega128")
        base = only(filter(r->r.case==case.label && r.variant=="window_1e-5",ratios))
        other = only(filter(r->r.case==case.label && r.variant==label,ratios))
        for field in (:Kplus_over_pi_plus, :Kminus_over_pi_minus)
            c = compare(getproperty(base,field),getproperty(other,field))
            retained_failure = case.label=="onset_145_375" && !isfinite(getproperty(base,field)) && !isfinite(getproperty(other,field))
            passed = c.passed || retained_failure; window_passed &= passed
            push!(comparisons, merge((case=case.label, variant=label, channel=String(field),
                required_window_check=true, retained_failure=retained_failure), c, (passed=passed,)))
        end
    end
    CSV.write(joinpath(output,"densities.csv"),rows)
    CSV.write(joinpath(output,"ratios.csv"),ratios)
    CSV.write(joinpath(output,"comparisons.csv"),comparisons)
    S._write_atomic(joinpath(output,"inputs.json"),(cases=inputs,solver_called=false))
    S._write_atomic(joinpath(output,"details.json"),details)
    output_hashes=Dict(name=>B.hashfile(joinpath(output,name)) for name in readdir(output))
    S._write_atomic(joinpath(output,"manifest.json"),(schema="charged_gbu_q0_window_audit_v1",
        generated_at_utc=string(Dates.now(Dates.UTC)), git_head=readchomp(`git -C $(S.ROOT) rev-parse HEAD`),
        source_hashes=merge(S.source_hashes(S.DEFAULT_CONFIG),
            Dict("scripts/analysis/relaxtime/audit_charged_gbu_q0_window.jl"=>B.hashfile(@__FILE__))),
        background_source=saved.provenance, cases=CASES, variants=variants(),
        window_checks_passed=window_passed, window_relative_tolerance=1e-3, window_absolute_tolerance=1e-12,
        output_hashes=output_hashes, density_prescription="finite_window_derivative",
        solver_called=false, production_authorized=false, diagnostic_only=true))
    window_passed || error("finite-window sensitivity check failed; inspect retained audit before scanning")
    return nothing
end

function main(args=ARGS)
    length(args)==4 && args[1]=="--input-root" && args[3]=="--output" ||
        throw(ArgumentError("use --input-root DIR --output DIR"))
    run(abspath(args[2]),abspath(args[4]))
end

end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && ChargedGBUQ0WindowAudit.main()
