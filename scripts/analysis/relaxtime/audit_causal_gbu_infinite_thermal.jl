"""Retained-background infinite-thermal candidate audit. Never updates production."""
module CausalGBUInfiniteThermalAudit
include("causal_gbu_infinite_yield.jl")
const Y=CausalGBUInfiniteYield
const P=Y.P
const I=Y.I
const R=Y.R
using CSV,JSON3

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_INFINITE_OUTPUT",joinpath(base,"direction_b_infinite_20260907"))
    qs=parse.(Float64,split(get(ENV,"GBU_INFINITE_Q","0,0.1,1.4,3.2,6,8"),','))
    channels=Symbol.(split(get(ENV,"GBU_INFINITE_CHANNELS","pi_plus,pi_minus,K_plus,K_minus"),','))
    mesh=parse(Int,get(ENV,"GBU_INFINITE_MESH","512"))
    hashes=R.start_output(output)
    rows,probes,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    for ch in channels,q in qs
        try
            k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=36.)
            p=P.profile(k;mesh=mesh)
            fine=P.profile(k;mesh=2mesh)
            zlist=[k.shift+0.8im,2+0.8im,complex(k.threshold),complex(k.threshold+0.01),
                complex((k.threshold+k.landau)/2)]
            for z in zlist
                raw=I.polarization(k,z;nodes=128)
                value=P.polarization(fine,z)
                radial=imag(z)>0 ? I.radial_polarization(bg,ch,q,z;nodes=192) : value
                push!(probes,(channel=String(ch),q_inv_fm=q,lambda_real=real(z),eta=imag(z),
                    raw_real=real(raw),raw_imag=imag(raw),profile_error=abs(value-raw),
                    mesh_change=abs(value-P.polarization(p,z)),radial_error=abs(radial-raw),
                    passed=max(abs(value-raw),abs(radial-raw))<1e-6))
            end
            a=Y.shell(p;nodes=64);b=Y.shell(fine;nodes=128)
            push!(rows,(channel=String(ch),q_inv_fm=q,density=b.density,bound=b.bound,
                landau=b.landau,pair=b.pair,roots=b.root_count,negative_roots=b.negative_root_count,
                density_change=abs(a.density-b.density),relative_change=abs(a.density-b.density)/max(abs(b.density),1e-12),
                threshold_phase=b.threshold_phase,landau_phase=b.landau_phase,high_phase=b.high_phase,
                static_inverse=b.static_inverse,omega_tail_conditional_bound=b.omega_tail_conditional_bound))
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(channel=String(ch),q_inv_fm=q,reason=sprint(showerror,err)))
        end
        isempty(rows) || CSV.write(joinpath(output,"shells.csv"),rows)
        isempty(probes) || CSV.write(joinpath(output,"probes.csv"),probes)
        isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
        println("[infinite] $(ch) q=$(q) shells=$(length(rows)) failures=$(length(failures))");flush(stdout)
    end
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    R.finish_output(output,hashes,Dict("background"=>bg,"shell_rows"=>length(rows),
        "failures"=>length(failures),"probe_rows"=>length(probes),"mesh"=>[mesh,2mesh],
        "thermal_target"=>"infinite_internal_momentum","solver_called"=>false,
        "probe_tolerance_inv_fm2"=>1e-6,"all_probes_passed"=>all(r.passed for r in probes),
        "full_UHP_certified"=>false,"q_integral_computed"=>false,
        "status"=>isempty(failures) && all(r.passed for r in probes) ? "local_candidate_passed" : "candidate_has_failures"))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
