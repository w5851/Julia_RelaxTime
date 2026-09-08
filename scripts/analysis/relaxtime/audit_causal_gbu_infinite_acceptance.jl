"""Automatic representative acceptance; no new equilibrium solves or baseline writes."""
module CausalGBUInfiniteAcceptanceAudit
include("causal_gbu_infinite_acceptance.jl")
const A=CausalGBUInfiniteAcceptance
const Y=A.Y
const P=A.P
const I=A.I
const R=A.R
using CSV,JSON3

function local_case(bg,ch,q,mesh)
    k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=36.)
    p=P.profile(k;mesh=mesh)
    fine=P.profile(k;mesh=2mesh)
    b=Y.shell(fine;nodes=96)
    coarse=Y.shell(p;nodes=64)
    contour=A.contour_audit(fine;nodes=32)
    contour2=A.contour_audit(fine;nodes=64,indent=3e-6,radius=30.)
    i,j=R.charged_rpa_spec(ch).pair
    zs=[complex((k.landau+k.threshold)/2),complex(k.threshold),2+0.8im]
    probe_error=maximum(abs(P.polarization(fine,z)-I.polarization(k,z;nodes=128)) for z in zs)
    finite=P.profile(Y.finite_kernel(bg,ch,q,10.;cut_nodes=96);mesh=2mesh)
    finiteval=q>0 ? Y.shell(finite;nodes=96).density : 0.
    etas=NamedTuple[]
    if q>0
        for eta in (0.003,0.001,0.0003),lower in (1e-3,3e-4)
            value=A.eta_shell(fine,eta;lower=lower,nodes=96)
            push!(etas,(channel=String(ch),q_inv_fm=q,eta_inv_fm=eta,lower_inv_fm=lower,
                density=value,pv_density=b.density,relative_difference=abs(value-b.density)/max(abs(b.density),1e-12)))
        end
    end
    row=(channel=String(ch),q_inv_fm=q,density=b.density,bound=b.bound,landau=b.landau,pair=b.pair,
        positive_roots=b.root_count,negative_roots=b.negative_root_count,
        profile_error_inv_fm2=probe_error,mesh_density_difference=abs(b.density-coarse.density),
        finite_L10_density=finiteval,finite_L10_difference=abs(finiteval-b.density),
        relative_L10_difference=abs(finiteval-b.density)/max(abs(b.density),1e-12),
        upper_half_plane_count=contour.count,contour_passed=contour.passed && contour2.passed,
        root_disks_passed=contour.root_disks_passed && contour2.root_disks_passed,
        contour_max_step=max(contour.max_step,contour2.max_step),
        tail_conditional_bound=b.omega_tail_conditional_bound,
        numerics_passed=probe_error<1e-6 && abs(b.density-coarse.density)<0.01max(abs(b.density),1e-12),
        global_UV_certified=false,production_authorized=false)
    return row,etas
end

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    tag=get(ENV,"GBU_INFINITY_BACKGROUND","T170")
    bg=tag=="T170" ? R.frozen_background(joinpath(base,"negative_density_phase_fig2_like")) :
        tag=="T205" ? R.saved_background(joinpath(base,"method_v1_mott"),205.625,24) :
        R.saved_background(joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3"),parse(Float64,tag),48)
    output=get(ENV,"GBU_INFINITY_OUTPUT",joinpath(base,"direction_b_infinite_acceptance_20260907"))
    qs=parse.(Float64,split(get(ENV,"GBU_INFINITY_Q","1.4,3.2"),','))
    channels=Symbol.(split(get(ENV,"GBU_INFINITY_CHANNELS","pi_plus,K_plus"),','))
    mesh=parse(Int,get(ENV,"GBU_INFINITY_MESH","512"))
    hashes=R.start_output(output)
    rows,etas,failures=NamedTuple[],NamedTuple[],NamedTuple[]
    for ch in channels,q in qs
        try
            row,e=local_case(bg,ch,q,mesh)
            push!(rows,row);append!(etas,e)
        catch err
            err isa InterruptException && rethrow()
            push!(failures,(channel=String(ch),q_inv_fm=q,reason=sprint(showerror,err)))
        end
        isempty(rows) || CSV.write(joinpath(output,"acceptance.csv"),rows)
        isempty(etas) || CSV.write(joinpath(output,"eta_limits.csv"),etas)
        isempty(failures) || CSV.write(joinpath(output,"failures.csv"),failures)
        println("[infinite-acceptance] $(ch) q=$(q) completed=$(length(rows)) failed=$(length(failures))");flush(stdout)
    end
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    passed=isempty(failures) && all(r.numerics_passed && r.contour_passed for r in rows)
    R.finish_output(output,hashes,Dict("status"=>passed ? "representative_numerical_acceptance" : "acceptance_failed",
        "background"=>bg,"rows"=>length(rows),"failures"=>length(failures),"eta_rows"=>length(etas),
        "passed"=>passed,"mesh"=>[mesh,2mesh],"solver_called"=>false,"thermal_target"=>"infinity",
        "global_UV_certified"=>false,"complete_freezeout_curve_authorized"=>false))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
