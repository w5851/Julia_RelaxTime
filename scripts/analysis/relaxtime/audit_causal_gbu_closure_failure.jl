"""Retain independent high-precision evidence for a failed support-gap audit."""
module CausalGBUClosureFailure
include("causal_gbu_regulator_checks.jl")
const C=CausalGBURegulatorChecks
const R=C.R
using CSV, JSON3

# Retained pre-fix real-axis expression, only for reproducing its cancellation.
function legacy_pv(p,w)
    x,y=p.energy,p.imaginary
    j=searchsortedlast(x,w)
    rho=1<=j<length(x) ? y[j]+(y[j+1]-y[j])*(w-x[j])/(x[j+1]-x[j]) : 0.0
    v=0.0
    for i in 1:length(x)-1
        a,b=x[i],x[i+1]
        slope=(y[i+1]-y[i])/(b-a)
        coefficient=y[i]+slope*(w-a)-rho
        v+=slope*(b-a)
        w!=a && w!=b && (v+=coefficient*log(abs((b-w)/(a-w))))
    end
    rho!=0 && (v+=rho*log(abs((last(x)-w)/(first(x)-w))))
    return complex(v/pi,rho)
end

# Independent high precision oracle for the SAME piecewise-linear spectrum.
# This is only evaluated in an analytic gap, where rho(w)=0 exactly.
function big_gap_value(profile,w)
    setprecision(256) do
        x,y=BigFloat.(profile.energy),BigFloat.(profile.imaginary)
        z=BigFloat(w)
        val=BigFloat(0)
        for i in 1:length(x)-1
            a,b=x[i],x[i+1]
            (y[i]==0 && y[i+1]==0) && continue
            (z<a || z>b) || error("oracle requires analytic gap")
            slope=(y[i+1]-y[i])/(b-a)
            val+=slope*(b-a)+(y[i]+slope*(z-a))*log(abs((b-z)/(a-z)))
        end
        return Float64(val/big(pi))
    end
end

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    output=get(ENV,"GBU_FAILURE_OUTPUT",joinpath(base,"regulator_closure_failure_probe_20260905"))
    bg=R.frozen_background(joinpath(base,"negative_density_phase_fig2_like"))
    hashes=R.start_output(output)
    b=C.centered_profile(bg,:K_plus,3.;mesh=256,ne=128)
    cuts=C.support_intervals(3.,bg.m.u,bg.m.s,bg.vacuum,24.,:centered)
    gap=first(C.analytic_gaps(cuts,b.shift))
    records,scan_status=NamedTuple[],NamedTuple[]
    for variant in (:legacy,:stable),nr in (64,128,256,512)
        inverse=variant===:stable ? b.inverse : w->1-4bg.coupling[:K_plus]*legacy_pv(b.profile,w+b.shift)
        result=R.certify_gap_roots((w,_)->inverse(w),3.,[(gap[1]+1e-6,gap[2]-1e-6)];
            physical_sheet=true,real_axis=true,omega_nodes=nr)
        # Reproduce the original domain, including its doubled margin; do not
        # conflate the arithmetic comparison with the separately fixed margin.
        push!(scan_status,(variant=String(variant),nodes=nr,count=result.count,passed=result.passed,
            status=String(result.status),rejected=join((String(r.status) for r in result.rejected),';'),
            effective_margin_inv_fm=2e-6))
        for root in result.roots
            w=root.omega_inv_fm
            oracle=1-4bg.coupling[:K_plus]*big_gap_value(b.profile,w+b.shift)
            push!(records,(variant=String(variant),nodes=nr,root_inv_fm=w,float_residual=root.residual,
                oracle_inverse=oracle,slope=root.slope,passed=result.passed))
        end
    end
    probes=NamedTuple[]
    for w in 3.7810782 .+ (-1e-6,-1e-7,-1e-8,0.,1e-8,1e-7,1e-6)
        v=real(b.inverse(w))
        oracle=1-4bg.coupling[:K_plus]*big_gap_value(b.profile,w+b.shift)
        push!(probes,(k0_inv_fm=w,float_inverse=v,oracle_inverse=oracle,error=v-oracle))
    end
    x,y=b.profile.energy,b.profile.imaginary
    order=sortperm(diff(x))
    cells=[(left=x[i],right=x[i+1],width=x[i+1]-x[i],yl=y[i],yr=y[i+1],slope=(y[i+1]-y[i])/(x[i+1]-x[i])) for i in order[1:20]]
    CSV.write(joinpath(output,"roots.csv"),records)
    CSV.write(joinpath(output,"scan_status.csv"),scan_status)
    CSV.write(joinpath(output,"gap_values.csv"),probes)
    CSV.write(joinpath(output,"small_cells.csv"),cells)
    R.finish_output(output,hashes,Dict("status"=>"diagnostic_failure_probe","background"=>bg,
        "density_computed"=>false,"oracle"=>"BigFloat256 same piecewise spectrum analytic gap"))
    println(scan_status);println(probes);println(first(cells))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
