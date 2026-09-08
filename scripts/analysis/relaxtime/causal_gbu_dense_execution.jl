"""Analysis-only scan execution helpers; preserve the retained integration mesh."""
module CausalGBUDenseExecution
include("causal_gbu_research_utils.jl")
const R=CausalGBUResearch

"""Reference shell without building a finite-q loop that is subsequently discarded.

This duplicates the reference branch of R.shell deliberately as an auditable
execution adapter. Numeric parity against that slow oracle is tested.
"""
function reference_shell(bg,channel,q,s,ref)
    b=ref.bubble
    inverse=w->begin
        v=R.q0_reference_coordinate(w,q,b.shift)
        v===nothing ? complex(1.0,0.0) : b.inverse(v)
    end
    a,c=R.charged_rpa_spec(channel).pair
    threshold=hypot(q,bg.m[a]+bg.m[c])-b.shift
    ld=hypot(q,bg.m[a]-bg.m[c])-b.shift
    roots=[hypot(r.omega_inv_fm+b.shift,q)-b.shift for r in ref.gap.roots]
    boundary_ok=abs(angle(b.inverse(-b.shift)))<s.phase_tol
    s.upper>threshold+0.01 || error("upper endpoint is below threshold")
    limit=R.threshold_phase_limit(inverse,threshold;phase_tol=s.phase_tol)
    uw=sort!(unique(vcat(threshold .+ [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3],
        collect(range(threshold+0.002,s.upper;length=s.nw)))))
    lw=ld>s.lower+1e-8 ? exp.(range(log(s.lower),log(ld-1e-8);length=s.nw)) : Float64[]
    onset=q-b.shift
    !isempty(lw) && s.lower<onset<ld && (lw=sort!(unique(vcat(lw,[onset,onset+1e-8]))))
    up,lp=[-angle(inverse(w)) for w in uw],[-angle(inverse(w)) for w in lw]
    branches=maximum(abs,diff(up))<pi && (length(lp)<2 || maximum(abs,diff(lp))<pi)
    u=R.bu_phase_integral_parts(uw,up,bg.T;weight=:gbu)
    l=isempty(lw) ? nothing : R.bu_phase_integral_parts(lw,lp,bg.T;weight=:gbu)
    bose_ok=all(>(0),roots)
    bound=bose_ok ? sum((1/expm1(w/bg.T) for w in roots);init=0.0) : NaN
    landau=l===nothing ? 0.0 : l.derivative
    lr=limit.phase/pi-length(roots)
    gates=(gap=ref.gap.passed,bose=bose_ok,branch=branches,levinson=abs(lr)<s.phase_tol,
        threshold_limit=limit.passed,upper=abs(last(up))<s.phase_tol,reference_boundary=boundary_ok)
    failed=join([String(k) for (k,v) in pairs(gates) if !v],';')
    return (channel=String(channel),q_inv_fm=Float64(q),route="q0_lambda_reference",bound=bound,
        unitary=u.derivative,landau=landau,shell_inv_fm2=q^2/(2pi^2)*(bound+u.derivative+landau),
        gap_count=length(roots),contour_count=ref.gap.contour_count,contour_passed=ref.gap.contour_passed,
        static_passed=ref.gap.static_passed,root_nodes_stable=ref.gap.root_nodes_stable,levinson_residual=lr,
        threshold_phase_over_pi=limit.phase/pi,threshold_offset=limit.offset,
        threshold_phase_at_1e7_over_pi=first(limit.phases)/pi,threshold_limit_change=limit.phase_change_over_pi,
        threshold_limit_passed=limit.passed,upper_phase=last(up),branch_passed=branches,
        reference_boundary_passed=boundary_ok,failed_gates=failed,passed=isempty(failed),
        status=isempty(failed) ? "conditional_window" : "gate_failed",production_authorized=false)
end

function reference_density(bg,channel,s)
    q,ws=R.gauleg(0.0,s.qmax,s.nq)
    b=R.bubble_at(bg,channel,0.0,s)
    ref=(bubble=b,gap=R.gap_audit(b,s))
    rows=[reference_shell(bg,channel,v,s,ref) for v in q]
    value=sum(ws[i]*rows[i].shell_inv_fm2 for i in eachindex(q))
    return (density=value,passed=all(r.passed for r in rows) && isfinite(value) && value>=0,rows=rows)
end
end
