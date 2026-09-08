"""Automatic Mott-aware outer integration checks, independent of phase unwrap."""
module CausalGBUInfiniteQGate
include("causal_gbu_infinite_acceptance.jl")
const A=CausalGBUInfiniteAcceptance
const Y=A.Y
const P=A.P
const I=A.I
const R=A.R

function mott_momenta(bg,ch;maximum_q=8.,nodes=64)
    function threshold_inverse(q)
        k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=max(36.,q+28bg.T+2.))
        return real(1-4bg.coupling[ch]*I.polarization(k,k.threshold;nodes=nodes))
    end
    roots=Y.bracket_roots(threshold_inverse,0.,maximum_q;nodes=16)
    checks=map(roots) do q
        k=I.kernel(bg,ch,q;cut_nodes=96,split_inv_fm=max(36.,q+28bg.T+2.))
        fine=real(1-4bg.coupling[ch]*I.polarization(k,k.threshold;nodes=2nodes))
        (q_inv_fm=q,threshold_inverse=fine,passed=abs(fine)<1e-6)
    end
    return (roots=roots,checks=checks,passed=all(r.passed for r in checks))
end

"""Check cut real-axis crossings of F using independent raw rho sign changes.

Between consecutive zeros of rho, F remains in one open half-plane. Positive
real F at non-gap crossings excludes an additional winding on these segments.
Root rectangles separately check both analytic gaps of this numerical kernel.
"""
function topology(p;scan_nodes=32)
    k=p.kernel;d,U=k.landau,k.threshold
    roots=Y.gap_roots(p)
    cs=Float64[]
    xs=sort!(unique!(vcat([a+(b-a)*sinpi(j/(2scan_nodes))^2
        for (a,b) in zip(k.edges[1:end-1],k.edges[2:end]) for j in 0:scan_nodes])))
    rs=k.rho.(xs)
    for j in 1:length(xs)-1
        rs[j]*rs[j+1]<0 || continue
        a,b=xs[j],xs[j+1]
        # Discontinuous q0 vacuum edge is checked by the explicit contour audit.
        k.q==0 && continue
        fl=rs[j]
        for _ in 1:48
            mid=(a+b)/2; fm=k.rho(mid)
            if fl*fm<=0; b=mid; else; a=mid;fl=fm;end
        end
        push!(cs,(a+b)/2)
    end
    values=[real(P.inverse(p,x)) for x in cs]
    reflect(z)=imag(z)>=0 ? P.inverse(p,z) : conj(P.inverse(p,conj(z)))
    windows=[Y.Limits.gap_window(a,b,roots[j]) for (j,(a,b)) in enumerate(((-U,-d),(d,U)))]
    disks=map(windows) do window
        a,b=window.left,window.right
        h=min(0.02,(b-a)/20)
        R.contour_count(reflect,a,b,-h,h;nodes=32,max_nodes=2048)
    end
    return (cut_crossings=length(cs),minimum_crossing_inverse=minimum(values;init=1.),
        positive_count=disks[2].count,negative_count=disks[1].count,
        passed=all(>(0),values) && all(disks[j].passed && disks[j].count==length(roots[j]) for j in 1:2),
        cutoff_scan_nodes=scan_nodes,uniform_interval_certificate=false)
end
end
