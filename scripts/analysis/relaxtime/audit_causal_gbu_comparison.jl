"""Same-input temp7/q0 and q0-extrapolation/finite-q comparisons. Diagnostic only."""
module CausalGBUComparison
include("causal_gbu_research_utils.jl")
using .CausalGBUResearch
using CSV, JSON3
const R = CausalGBUResearch

function settings(;kwargs...)
    d = R.method_contract()["numerics"]
    base = (mesh=d["segment_nodes"],np=d["momentum_nodes"],nx=d["angle_nodes"],ne=d["energy_nodes"],
        nw=d["omega_nodes"],nr=d["root_nodes"],nq=d["q_nodes"],qmax=d["qmax_inv_fm"],
        thermal=d["thermal_max_inv_fm"],lower=d["omega_min_inv_fm"],upper=d["omega_max_inv_fm"])
    return R.Settings(;merge(base,(;kwargs...))...)
end

function main()
    base = joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input = get(ENV,"GBU_RESEARCH_INPUT",joinpath(base,"negative_density_phase_fig2_like"))
    output = get(ENV,"GBU_RESEARCH_OUTPUT",joinpath(base,"method_v1_comparison"))
    old = get(ENV,"GBU_TEMP7_SOURCE","")
    isfile(old) || error("GBU_TEMP7_SOURCE must explicitly name audit_fig2_phase_shift.jl")
    old_hash = R.hashfile(old)
    bg = R.frozen_background(input)
    oracle = Module(:Temp7MatchedOracle)
    Core.eval(oracle, :(const RelaxTime = Main.RelaxTime))
    Base.include(oracle,old)
    hashes = R.start_output(output)
    rows = NamedTuple[]
    s = settings()
    ps,pws = R.gauleg(0.0,bg.vacuum,s.np)
    for channel in R.CHANNELS
        a,b = R.charged_rpa_spec(channel).pair
        shift = bg.mu[a]-bg.mu[b]
        matched = R.bubble_at(bg,channel,0.0,settings(thermal=bg.vacuum))
        extended = R.bubble_at(bg,channel,0.0,s)
        # Pi_temp7=pol_temp7/2 by the charged-ladder trace normalization.
        i1 = sum(Base.invokelatest(oracle.i1_lit,bg.m[f],bg.mu[f],bg.T,bg.Phi,bg.PhiBar,ps,pws) for f in (a,b))
        lo,hi = abs(bg.m[a]-bg.m[b]),bg.m[a]+bg.m[b]
        lambdas = sort!(unique(vcat(collect(range(lo+0.001,hi-0.001;length=32)),
            hi .+ [0.01,0.03,0.1,0.3,0.7,1.0], [0.07,0.17,0.37])))
        for lambda in lambdas
            w = lambda-shift
            i2 = Base.invokelatest(oracle.i2_lit_pv_b0,w,bg.m[a],bg.m[b],bg.mu[a],bg.mu[b],bg.T,bg.Phi,bg.PhiBar)
            oldpi = 2*(i1-(lambda^2-(bg.m[a]-bg.m[b])^2)*i2)
            newpi = R.cauchy_transform(matched.profile,lambda)
            extpi = R.cauchy_transform(extended.profile,lambda)
            dold,dnew = 1-4bg.coupling[channel]*oldpi,1-4bg.coupling[channel]*newpi
            push!(rows,(channel=String(channel),pair="$(a),$(b)",T_MeV=bg.T_MeV,
                mu1_inv_fm=bg.mu[a],mu2_inv_fm=bg.mu[b],q_inv_fm=0.0,k0_inv_fm=w,lambda_inv_fm=lambda,
                region=lambda<lo ? "landau" : lambda<hi ? "gap" : "unitary",
                temp7_pi_real=real(oldpi),temp7_pi_imag=imag(oldpi),new_matched_pi_real=real(newpi),
                new_matched_pi_imag=imag(newpi),real_difference=real(newpi-oldpi),imag_difference=imag(newpi-oldpi),
                new_extended_pi_real=real(extpi),new_extended_pi_imag=imag(extpi),
                thermal_prescription_difference=abs(extpi-newpi),
                temp7_minus_arg_inverse=-angle(dold),new_minus_arg_inverse=-angle(dnew),
                production_authorized=false))
        end
        println("[gbu-comparison] matched q0 $(channel)")
    end
    CSV.write(joinpath(output,"matched_q0.csv"),rows)
    summaries,shellrows = NamedTuple[],NamedTuple[]
    for channel in R.CHANNELS, ref in (false,true)
        result = R.density(bg,channel,s;reference=ref)
        append!(shellrows,result.rows)
        push!(summaries,(channel=String(channel),route=ref ? "q0_lambda_reference" : "direct_finite_q",
            density_inv_fm3=result.density,passed=result.passed,failed_shells=count(r->!r.passed,result.rows),
            production_authorized=false))
        CSV.write(joinpath(output,"density_comparison.csv"),summaries)
        CSV.write(joinpath(output,"shell_comparison.csv"),shellrows)
        println("[gbu-comparison] $(channel) reference=$(ref) n=$(result.density) pass=$(result.passed)")
    end
    R.hashfile(old)==old_hash || error("temp7 source changed")
    R.finish_output(output,hashes,Dict("status"=>"research_comparison","background"=>bg,
        "settings"=>s,"temp7_source"=>old,"temp7_sha256"=>old_hash,
        "temp7_dependency"=>"current retained OneLoopIntegrals.B0, not a historical binary",
        "q0_matching"=>"all flavors, masses, T, chemical potentials, Phi/PhiBar, K and vacuum/thermal cutoff equal",
        "normalization"=>"Pi_project=pol_temp7/2; D_project=1-4K*Pi; no fitted sign or fold",
        "reference"=>"lambda-invariant timelike q0 extrapolation; not literal external-k0 temp7 mapping",
        "input_hashes"=>Dict(f=>R.hashfile(joinpath(input,f)) for f in ("plot_manifest.json","charged_phase_profile_detail.csv")),
        "solver_called"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
