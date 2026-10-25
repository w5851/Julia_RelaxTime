#!/usr/bin/env julia
using CSV, SHA
if !isdefined(Main,:TerminalTransportQuadrature)
    include(joinpath(@__DIR__,"TerminalTransportQuadrature.jl"))
end
using .TerminalTransportQuadrature, .ThresholdUpgradeCommon

function terminal_configs(profile)
    profile in ("cutoff", "cost", "historical") || throw(ArgumentError("unknown terminal profile: $profile"))
    reference=[(name="semi128_c48",n=128,pmax=10.0,c=48,semi=true),
        (name="semi256_c64",n=256,pmax=10.0,c=64,semi=true)]
    profile=="historical" && return vcat([
        (name="finite8_n24_c8",n=24,pmax=8.0,c=8,semi=false),
        (name="finite8_n96_c48",n=96,pmax=8.0,c=48,semi=false)],reference)
    profile=="cost" && return vcat([
        (name="semi16_c8",n=16,pmax=10.0,c=8,semi=true),
        (name="semi24_c12",n=24,pmax=10.0,c=12,semi=true),
        (name="semi32_c16",n=32,pmax=10.0,c=16,semi=true),
        (name="semi48_c24",n=48,pmax=10.0,c=24,semi=true),
        (name="semi64_c32",n=64,pmax=10.0,c=32,semi=true),
        (name="semi32_c32",n=32,pmax=10.0,c=32,semi=true),
        (name="semi64_c16",n=64,pmax=10.0,c=16,semi=true),
        (name="semi32_c16_scale5",n=32,pmax=5.0,c=16,semi=true),
        (name="semi32_c16_scale20",n=32,pmax=20.0,c=16,semi=true)],reference)
    return vcat([(name="finite15_n16_c4",n=16,pmax=15.0,c=4,semi=false),
        (name="finite15_n32_c16",n=32,pmax=15.0,c=16,semi=false),
        (name="finite15_n64_c32",n=64,pmax=15.0,c=32,semi=false),
        (name="finite15_n128_c64",n=128,pmax=15.0,c=64,semi=false),
        (name="finite20_n96_c48",n=96,pmax=20.0,c=48,semi=false),
        (name="finite25_n128_c64",n=128,pmax=25.0,c=64,semi=false)],reference)
end

function terminal_validation(inputdir,outdir;profile="cutoff")
    configs=terminal_configs(profile)
    start_run(outdir,inputdir,@__FILE__,(purpose="terminal RTA moment and bulk Gram convergence",
        profile=profile,configurations=configs,relative_target=0.01,background_and_tau="fixed retained production state",
        bulk_test="largest relative quadratic-form difference; every bracket at the fixed masses",
        reference="semi-infinite 256 momentum nodes and 64 polar nodes; checked against 128/48"),
        extra_sources=[joinpath(@__DIR__,"TerminalTransportQuadrature.jl")])
    rows=collect(CSV.File(joinpath(inputdir,"review_points.csv")))
    results=NamedTuple[];gram_rows=NamedTuple[];failures=NamedTuple[]
    for row in rows
        try
            qp=(m=(u=Float64(row.m_u),d=Float64(row.m_d),s=Float64(row.m_s)),
                μ=(u=Float64(row.muq_fm),d=Float64(row.muq_fm),s=Float64(row.muq_fm)))
            tp=(T=Float64(row.T_fm),Φ=Float64(row.Phi),Φbar=Float64(row.Phibar),ξ=Float64(row.xi))
            tau=NamedTuple{SPECIES}(Tuple(Float64(getproperty(row,Symbol(:tau_,s))) for s in SPECIES))
            calculated=map(configs) do config
                p,pw=config.semi ? AS._cached_semi_infinite_momentum_grid(config.n,config.pmax) :
                    Main.RelaxTime.GaussLegendre.gauleg(0.0,config.pmax,config.n)
                ct,cw=Main.RelaxTime.GaussLegendre.gauleg(-1.0,1.0,config.c)
                terminal_grams(qp,tp,p,pw,ct,cw)
            end
            reference=last(calculated)
            ref=transport_from_grams(reference,qp,tp,tau)
            for (config,grams) in zip(configs,calculated)
                value=transport_from_grams(grams,qp,tp,tau)
                eta_bound=maximum(abs(getproperty(grams,s)[1,1]/getproperty(reference,s)[1,1]-1) for s in SPECIES)
                sigma_bound=maximum(abs(getproperty(grams,s)[1,2]/getproperty(reference,s)[1,2]-1) for s in SPECIES)
                bulk_bound=maximum(gram_relative_change(getproperty(grams,s),getproperty(reference,s)) for s in SPECIES)
                push!(results,(point_id=String(row.point_id),muB_MeV=Float64(row.muB_MeV),
                    alpha_T=Float64(row.alpha_T),xi=Float64(row.xi),configuration=config.name,
                    eta=value.eta,sigma=value.sigma,eta_reference=ref.eta,sigma_reference=ref.sigma,
                    retained_eta=Float64(row.eta),retained_sigma=Float64(row.sigma),
                    eta_relative_difference=abs(value.eta/ref.eta-1),
                    sigma_relative_difference=abs(value.sigma/ref.sigma-1),
                    eta_any_positive_tau_bound=eta_bound,sigma_any_positive_tau_bound=sigma_bound,
                    bulk_any_bracket_bound=bulk_bound,
                    all_transport_moments_within_target=max(eta_bound,sigma_bound,bulk_bound)<=0.01))
                for s in SPECIES,i in 1:3,j in i:3
                    push!(gram_rows,(point_id=String(row.point_id),configuration=config.name,
                        species=String(s),i=i,j=j,moment=getproperty(grams,s)[i,j]))
                end
            end
            write_rows(joinpath(outdir,"moments.csv"),results)
            write_rows(joinpath(outdir,"grams.csv"),gram_rows)
            checked=results[end-length(configs)+1:end]
            candidate_name=profile=="cost" ? "semi32_c16" :
                profile=="historical" ? "finite8_n24_c8" : "finite15_n128_c64"
            candidate=only(filter(r->r.configuration==candidate_name,checked))
            reference_check=only(filter(r->r.configuration=="semi128_c48",checked))
            println("POINT ",row.point_id," ",candidate_name,"_bulk_bound=",candidate.bulk_any_bracket_bound,
                " semi128_bound=",reference_check.bulk_any_bracket_bound);flush(stdout)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(point_id=String(row.point_id),error=sprint(showerror,err,catch_backtrace())))
            write_rows(joinpath(outdir,"failures.csv"),failures)
        end
    end
    reference_checks=filter(r->r.configuration=="semi128_c48",results)
    finish_run(outdir;points=length(rows),failures=length(failures),
        semi_infinite_reference_checked=isempty(failures) && !isempty(rows) &&
            length(reference_checks)==length(rows) &&
            all(r->r.all_transport_moments_within_target,reference_checks))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    2<=length(ARGS)<=3 || error("Usage: run_terminal_transport_validation.jl INPUT_DIR NEW_OUTPUT_DIR [cutoff|cost|historical]")
    terminal_validation(abspath(ARGS[1]),abspath(ARGS[2]);profile=length(ARGS)==3 ? ARGS[3] : "cutoff")
end
