"""Dense diagnostic freezeout scan for the charged GBU route.

This script extends the retained three-point FixedMuBConservedCharges scan with
intermediate energies.  It writes a new directory only and keeps both direct
finite-q and q=0 lambda-reference routes explicit.  The output is diagnostic
and never a production baseline.
"""
module CausalGBUDenseFreezeout

include("audit_causal_gbu_validation.jl")
using .CausalGBUValidation
const V = CausalGBUValidation
const R = V.CausalGBUResearch
using CSV, JSON3, SHA
using Main.Models

const ENERGIES = [3.0, 5.0, 7.7, 11.5, 19.6, 27.0, 39.0, 62.4, 130.0, 200.0]

function energy_grid(text)
    e = parse.(Float64, split(text, ','))
    all(x -> isfinite(x) && x > 0, e) && length(unique(e)) == length(e) ||
        throw(ArgumentError("energies must be finite, positive and unique"))
    return sort(e)
end

"""Keep failed or missing energies as NaN gaps; never silently omit a point."""
function ratio_table(direct, reference, energies)
    rows = NamedTuple[]
    for e in energies
        groups = [filter(r -> r.sqrt_s_NN_GeV == e, source) for source in (direct, reference)]
        for group in groups
            length(unique(r.channel for r in group)) == length(group) || error("duplicate channel at $(e)")
        end
        maps = [Dict(r.channel => r for r in group) for group in groups]
        function ratio(m, numerator, denominator)
            haskey(m, numerator) && haskey(m, denominator) || return (NaN, false)
            a, b = m[numerator], m[denominator]
            v = a.density_inv_fm3 / b.density_inv_fm3
            ok = a.passed && b.passed && isfinite(v) && a.density_inv_fm3 >= 0 && b.density_inv_fm3 > 0
            return (v, ok)
        end
        dp, dpg = ratio(maps[1], "K_plus", "pi_plus")
        dm, dmg = ratio(maps[1], "K_minus", "pi_minus")
        rp, rpg = ratio(maps[2], "K_plus", "pi_plus")
        rm, rmg = ratio(maps[2], "K_minus", "pi_minus")
        push!(rows, (sqrt_s_NN_GeV=e, direct_Kplus_over_pi_plus=dp,
            reference_Kplus_over_pi_plus=rp, direct_Kminus_over_pi_minus=dm,
            reference_Kminus_over_pi_minus=rm, direct_plus_passed=dpg, direct_minus_passed=dmg,
            reference_plus_passed=rpg, reference_minus_passed=rmg, production_authorized=false))
    end
    return rows
end

function main()
    base = joinpath(R.ROOT, "data", "outputs", "results", "relaxtime", "analysis", "charged_rpa_phase_backend")
    output = get(ENV, "GBU_DENSE_OUTPUT", joinpath(base, "fig4_like_freezeout_ratio_dense_20260905_v2"))
    energies = energy_grid(get(ENV, "GBU_DENSE_ENERGIES", join(ENERGIES, ',')))
    R.method_contract()
    hashes = R.start_output(output)
    model = Models.create_model(:PNJL)
    profile = Models.load_freezeout_profile(profile="default")
    points = Models.build_freezeout_scan_points(energies; profile=profile, traversal=:sqrts_descending)
    settings = R.Settings(mesh=512, ne=128, nw=4800, nq=24, lower=1e-5, upper=56.0, thermal=24.0)
    cache = Dict{Tuple{Float64,Float64},Any}()
    backgrounds = NamedTuple[]
    direct = NamedTuple[]
    reference = NamedTuple[]
    shells_direct = NamedTuple[]
    shells_reference = NamedTuple[]
    failures = NamedTuple{(:sqrt_s_NN_GeV,:channel,:route,:reason,:production_authorized),Tuple{Float64,String,String,String,Bool}}[]
    manifest = Dict{String,Any}("schema" => "charged_gbu_dense_freezeout_v2", "status" => "running",
        "energies_GeV" => energies, "background" => "FixedMuBConservedCharges quark-only BQS; rhoQ/rhoB=0.4; rhoS=0",
        "settings" => settings, "routes" => ["direct_finite_q", "q0_lambda_reference"],
        "solver_called" => true, "production_authorized" => false, "complete_curve_authorized" => false,
        "source_hashes" => hashes, "background_solves" => 0)
    function checkpoint()
        for (name, data) in (("backgrounds", backgrounds), ("direct_densities", direct),
                ("reference_densities", reference), ("direct_shells", shells_direct), ("reference_shells", shells_reference))
            isempty(data) || CSV.write(joinpath(output, name*".csv"), data)
        end
        CSV.write(joinpath(output, "failures.csv"), failures)
        CSV.write(joinpath(output, "ratio_comparison.csv"), ratio_table(direct, reference, energies))
        manifest["background_solves"] = length(cache)
        manifest["failed_evaluations"] = length(failures)
        manifest["completed_density_evaluations"] = length(direct) + length(reference)
        open(joinpath(output, "progress.json"), "w") do io
            JSON3.write(io, manifest)
        end
    end
    checkpoint()
    for pt in points
        key = (Float64(pt.T_MeV), Float64(pt.muB_MeV))
        bg = if haskey(cache, key)
            cache[key].bg
        else
            candidate = isempty(cache) ? nothing : argmin(k -> abs(k[1]-pt.T_MeV) + abs(k[2]-pt.muB_MeV), collect(keys(cache)))
            seed = candidate === nothing ? nothing : cache[candidate].seed
            try
                r = V.background(model, pt.T_MeV, pt.muB_MeV, 48; seed=seed)
                cache[key] = r
                push!(backgrounds, r.row)
                checkpoint()
                println("[gbu-dense-background] s=$(pt.sqrt_s_NN_GeV) residual=$(r.bg.residual)")
                flush(stdout)
                r.bg
            catch err
                err isa InterruptException && rethrow()
                push!(failures, (sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV, channel="all", route="background", reason=sprint(showerror, err), production_authorized=false))
                checkpoint()
                continue
            end
        end
        for channel in R.CHANNELS
            for (route, isref, target, shells) in (("direct_finite_q", false, direct, shells_direct), ("q0_lambda_reference", true, reference, shells_reference))
                try
                    d = R.density(bg, channel, settings; reference=isref)
                    append!(shells, [merge((sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV, background_p=48), r) for r in d.rows])
                    push!(target, (sqrt_s_NN_GeV=Float64(pt.sqrt_s_NN_GeV), T_MeV=Float64(pt.T_MeV), muB_MeV=Float64(pt.muB_MeV), background_p=48,
                        channel=String(channel), route=route, density_inv_fm3=d.density, passed=d.passed,
                        failed_shells=count(r -> !r.passed, d.rows), production_authorized=false))
                    println("[gbu-dense] s=$(pt.sqrt_s_NN_GeV) $(route) $(channel) n=$(d.density) pass=$(d.passed)")
                catch err
                    err isa InterruptException && rethrow()
                    push!(failures, (sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV, channel=String(channel), route=route, reason=sprint(showerror, err), production_authorized=false))
                    println(stderr, "[gbu-dense-failed] s=$(pt.sqrt_s_NN_GeV) $(route) $(channel): $(sprint(showerror, err))")
                end
                checkpoint()
                flush(stdout)
            end
        end
    end
    manifest["status"] = "diagnostic_completed"
    manifest["all_conditional_gates"] = length(direct)==length(reference)==4length(energies) &&
        isempty(failures) && all(r.passed for r in vcat(direct, reference))
    checkpoint()
    R.finish_output(output, hashes, manifest)
    println("[gbu-dense] output: $(output)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
end
