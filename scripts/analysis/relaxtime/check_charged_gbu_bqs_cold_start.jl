"""Diagnose low-density cold-start stopping error using the SAME joint BQS residual.

Only analysis tolerances are tightened. No production solver or gate is changed,
and the original failed cross-check remains immutable evidence.
"""
module ChargedGBUBQSColdStart
using JSON3, LinearAlgebra, NLsolve
include("compute_charged_gbu_bqs_point.jl")
const B=ChargedGBUBQSPoint
const M=B.M
const W=B.W

function main(input,output)
    ispath(output) && error("refusing to overwrite $output")
    nominal=W.readjson(joinpath(input,"background.json"))
    cold=W.readjson(joinpath(input,"background_cold.json"))
    model=M.create_model(:PNJL)
    h=Main.Constants_PNJL.ħc_MeV_fm
    pt=(T_MeV=Float64(nominal.bg.T_MeV),muB_MeV=Float64(nominal.bg.muB_MeV))
    T=pt.T_MeV/h
    mode=M.FixedMuBConservedCharges(pt.muB_MeV/h,.4,0.)
    params=M.GapParams(T,M.cached_nodes(48,8),0.;p_num=48,t_num=8,model_kind=:PNJL)
    residual! = M.build_residual!(mode,params)
    config=W.validate_config(B.TOML.parsefile(W.DEFAULT_CONFIG))
    rows=NamedTuple[]
    for (name,b) in (("cold",cold),("continuation",nominal))
        x0=Float64.(b.seed);f0=zeros(8);residual!(f0,x0)
        # Preserve unified eight-dimensional solve; tighter ftol, never relaxed.
        sol=NLsolve.nlsolve(residual!,x0;autodiff=:forward,method=:trust_region,
            xtol=1e-14,ftol=1e-14,iterations=1000)
        x=Float64.(sol.zero);f=zeros(8);residual!(f,x)
        verified=W.background(model,pt,config;seed=x)
        d=B.diagnostics(model,verified,192)
        initial_rho=M.model_rho(model,x0[1:5],x0[6:8],T;p_num=192,t_num=8,xi=0.)
        initial_charges=M.conserved_densities_from_flavor(initial_rho)
        push!(rows,(start=name,initial_residual=f0,final_residual=f,
            initial_rho_B=initial_charges.rho_B,initial_rho_S=initial_charges.rho_S,
            initial_Q_over_B=initial_charges.rho_Q/initial_charges.rho_B,
            f_converged=sol.f_converged,iterations=sol.iterations,
            initial_solution=x0,polished_solution=x,verified_solution=verified.seed,
            maximum_difference_from_nominal=maximum(abs,verified.seed-Float64.(nominal.seed)),
            diagnostics=d))
    end
    delta=maximum(abs,rows[1].verified_solution-rows[2].verified_solution)
    passed=delta<1e-6 && all(r->r.f_converged && norm(r.final_residual)<1e-12 &&
        r.maximum_difference_from_nominal<1e-6 && abs(r.diagnostics.Q_over_B-.4)<1e-4 &&
        abs(r.diagnostics.rho_S_fm3)<1e-7,rows)
    result=(passed=passed,maximum_polished_seed_difference=delta,rows=rows,
        ftol=1e-14,xtol=1e-14,original_solver_unchanged=true,original_failed_check_preserved=true,
        input_hashes=Dict(f=>W.hashfile(joinpath(input,f)) for f in ("background.json","background_cold.json")),
        scope="same_8D_residual_cold_start_polish_not_a_new_density_integral_or_global_branch_proof")
    W.writejson(output,result)
    println(JSON3.write(result))
    return passed
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && (main(ARGS...) || exit(1))
end
