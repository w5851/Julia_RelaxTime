"""Direction-B static prerequisites on retained BQS input, never a solver."""
module CausalGBUDirectionBAudit
include("causal_gbu_direction_b.jl")
const B=CausalGBUDirectionB
const R=B.R
isdefined(Main,:Models) || Base.include(Main,joinpath(R.ROOT,"src","models","Models.jl"))
using Main.Models
using StaticArrays,ForwardDiff,CSV,JSON3

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    output=get(ENV,"GBU_DIRECTION_B_OUTPUT",joinpath(base,"direction_b_static_20260906"))
    hashes=R.start_output(output)
    rows=NamedTuple[]
    mvec=Float64[bg.m...]
    muvec=Float64[bg.mu...]
    for occupation in ("fermi","pnjl"),tc in (bg.vacuum,10.,24.),n in (64,128)
        phi,bar=occupation=="fermi" ? (1.,1.) : (bg.Phi,bg.PhiBar)
        params=merge(Main.Constants_PNJL.pnjl_constants(),(thermal_p_max_inv_fm=tc,))
        model=Models.PNJLModel(params)
        omega(v)=Models.vacuum_contribution(model,SVector{3}(v[1:3]))+
            Models.thermal_contribution(model,SVector{3}(v[1:3]),phi,bar,SVector{3}(v[4:6]),bg.T;
                p_num=n,t_num=8,xi=0.)
        x=vcat(mvec,muvec)
        d=ForwardDiff.gradient(omega,x)
        h=ForwardDiff.hessian(omega,x)
        states=[B.static_flavor(mvec[i],muvec[i],bg.T,phi,bar,bg.vacuum,tc;nodes=n) for i in 1:3]
        potential_error=abs(omega(x)-sum(s.omega_inv_fm4 for s in states))
        ps,ws=R.gauleg(0.,tc,n)
        for i in 1:3
            single=Main.RelaxTime.OneLoopIntegrals.A(mvec[i],muvec[i],bg.T,phi,bar,ps,ws)
            cond_from_A=3mvec[i]*single/(4pi^2)
            mass_error=abs(d[i]-cond_from_A)
            density_error=abs(-d[i+3]-states[i].density_inv_fm3)
            mixed_error=abs(h[i,i+3]-h[i+3,i])
            derivative_error=abs(d[i]-states[i].condensate_inv_fm3)
            push!(rows,(occupation=occupation,thermal_cutoff_inv_fm=tc,nodes=n,
                flavor=String((:u,:d,:s)[i]),potential_error_inv_fm4=potential_error,
                mass_derivative_error_inv_fm3=mass_error,density_derivative_error_inv_fm3=density_error,
                independent_derivative_error_inv_fm3=derivative_error,mixed_derivative_error_inv_fm2=mixed_error,
                condensate_inv_fm3=cond_from_A,density_inv_fm3=states[i].density_inv_fm3,
                passed=potential_error<1e-8 && max(mass_error,density_error,derivative_error)<1e-8 && mixed_error<1e-8,
                is_stationarity_test=false,is_dynamic_hessian_test=false,solver_called=false,production_authorized=false))
        end
        println("[direction-b] $(occupation) thermal=$(tc) n=$(n) static identities checked")
    end
    CSV.write(joinpath(output,"static_identities.csv"),rows)
    toy=B.toy_uhp_pole(3/(8pi),49.,0.24)
    open(joinpath(output,"synthetic_subtraction.json"),"w") do io
        JSON3.write(io,merge(toy,(c=3/(8pi),S_inv_fm2=49.,K_fm2=0.24,
            model="linear_synthetic_tail_only_not_project_bubble")))
    end
    paperdir=joinpath(R.ROOT,"tmp","direction_b_review_20260906")
    papers=Dict(joinpath(paperdir,f)=>R.hashfile(joinpath(paperdir,f)) for f in
        ("2105.14323.pdf","2102.02844.pdf","hep-ph_9509363.pdf"))
    all(R.hashfile(joinpath(bg.input_directory,p))==h for (p,h) in bg.input_hashes) || error("input drift")
    R.finish_output(output,hashes,Dict("status"=>"direction_b_prerequisite_diagnostic","background"=>bg,
        "static_identity_rows"=>length(rows),"all_static_identities_passed"=>all(r.passed for r in rows),
        "paper_hashes"=>papers,"full_stationarity_certified"=>false,"direction_b_response_selected"=>false,
        "upstream_thermal_extension_ruled_out"=>false,"synthetic_example_is_project_pole"=>false,
        "solver_called"=>false,"meson_density_computed"=>false,"production_authorized"=>false,
        "limitations"=>["Static differentiability does not uniquely determine a dynamic response",
            "Positive restored cut and matched low coefficients do not guarantee UHP stability",
            "No actual response regulator, matching coefficients or upstream equilibrium was changed"]))
    println("[direction-b] $(count(r->r.passed,rows))/$(length(rows)) identities; no response selected")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
