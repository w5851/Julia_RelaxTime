"""Signed real-axis gap checks for a compact spectral interpolant.

Completeness here is for ALL real gaps of the specified interpolant at one q,
not for cut-embedded zeros, UHP instabilities, or an exact continuum kernel.
The root count does not inspect or unwrap the density phase.
"""
module CausalGBUGapCompleteness
include("causal_gbu_regulator_checks.jl")
const C=CausalGBURegulatorChecks
const R=C.R
const CS=Main.RelaxTime.CausalSpectralBubble

"""Decompose sign-changing cells into endpoint triangles, without new knots.

u*(b-x)/(b-a)+v*(x-a)/(b-a) is exact even when the sign crossing cannot be
represented between adjacent Float64 energies. Neither endpoint is clipped.
"""
function signed_cells(p)
    cells=NTuple{4,Float64}[]
    for i in 1:length(p.energy)-1
        a,b,u,v=p.energy[i],p.energy[i+1],p.imaginary[i],p.imaginary[i+1]
        u==v==0 && continue
        if u!=0 && v!=0 && signbit(u)!=signbit(v)
            push!(cells,(a,b,u,0.),(a,b,0.,v))
        else
            push!(cells,(a,b,u,v))
        end
    end
    return cells
end

function cell_real(cell,w)
    a,b,u,v=cell
    if w==a
        u==0 || return copysign(Inf,u+v)
        return v-u
    elseif w==b
        v==0 || return -copysign(Inf,u+v)
        return v-u
    end
    a<w<b && throw(ArgumentError("interval touches a nonzero spectral cell"))
    t=(b-a)/(a-w)
    lg=log1p(t)
    return u*lg+(v-u)*CS._cauchy_cell_moment(t,lg)
end

"""Cauchy enclosure on a real gap: each sign-definite cell is monotone there.

This is a conservative Float64 numerical bound, with an explicit rounding
allowance, not an interval-arithmetic proof of the un-interpolated kernel.
"""
function inverse_enclosure(cells,K,left,right)
    all(isfinite,(K,left,right)) && left<=right || throw(ArgumentError("invalid gap enclosure"))
    low,high,size=0.,0.,0.
    for cell in cells
        a,b,_,_=cell
        max(a,left)<min(b,right) && throw(ArgumentError("enclosure is not an analytic gap"))
        u,v=cell_real(cell,left),cell_real(cell,right)
        low+=min(u,v)/pi
        high+=max(u,v)/pi
        size+=max(abs(u),abs(v))/pi
    end
    if !all(isfinite,(low,high,size))
        return (lower=-Inf,upper=Inf,roundoff=Inf,excludes_zero=false)
    end
    error=256eps(Float64)*max(size,1.)*max(length(cells),1)
    lower,upper=extrema((1-4K*(low-error),1-4K*(high+error)))
    return (lower=lower,upper=upper,roundoff=4abs(K)*error,
        excludes_zero=lower>0 || upper<0)
end

"""Outside |lambda|>=R, |4K Pi|<=4|K| integral|rho|/pi/(R-U)<1."""
function tail_exclusion(cells,K)
    isfinite(K) || throw(ArgumentError("finite coupling required"))
    U=maximum((max(abs(c[1]),abs(c[2])) for c in cells);init=0.)
    mass=sum(c -> (c[2]-c[1])*(abs(c[3])+abs(c[4]))/(2pi),cells;init=0.)
    mass*=1+256eps(Float64)*max(length(cells),1)
    radius=U+max(1.,8abs(K)*mass)
    bound=4abs(K)*mass/(radius-U)
    return (radius_inv_fm=radius,support_radius_inv_fm=U,l1_mass_inv_fm3=mass,
        inverse_deviation_bound=bound,passed=bound<1)
end

function profile_support(cells)
    return C.merge_intervals([(c[1],c[2]) for c in cells])
end

function uncovered_width(intervals,cuts)
    excess=0.
    for (left,right) in intervals
        pieces=[(left,right)]
        for (a,b) in cuts
            next=Tuple{Float64,Float64}[]
            for (l,r) in pieces
                l<min(a,r) && push!(next,(l,min(a,r)))
                max(b,l)<r && push!(next,(max(b,l),r))
            end
            pieces=next
        end
        excess=max(excess,maximum((r-l for (l,r) in pieces);init=0.))
    end
    return excess
end

"""Require geometry and interpolant support to agree, in BOTH directions.

Missing numerical spectral weight must not create a supposedly physical gap.
Only coordinate-roundoff-sized endpoint differences are admitted and reported.
"""
function support_geometry_check(cells,cuts)
    represented=profile_support(cells)
    excess=uncovered_width(represented,cuts)
    missing=uncovered_width(cuts,represented)
    scale=maximum((max(abs(c[1]),abs(c[2])) for c in cells);init=1.)
    allowance=128eps(Float64)*max(scale,1.)
    return (max_excess_width_inv_fm=excess,max_missing_width_inv_fm=missing,
        roundoff_allowance_inv_fm=allowance,passed=max(excess,missing)<=allowance)
end

function count_all_real_gaps(p,K,q;root_nodes=128,margin=1e-6)
    root_nodes>=8 && isfinite(margin) && margin>0 || throw(ArgumentError("invalid root settings"))
    cells=signed_cells(p)
    tail=tail_exclusion(cells,K)
    cuts=profile_support(cells)
    gaps=C.analytic_gaps(cuts,0.;lower=-tail.radius_inv_fm,upper=tail.radius_inv_fm)
    inverse(z)=1-4K*R.cauchy_transform(p,z)
    reflection(z)=imag(z)>=0 ? inverse(z) : conj(inverse(conj(z)))
    rows=NamedTuple[]
    for (id,(left,right)) in enumerate(gaps)
        enclosure=inverse_enclosure(cells,K,left,right)
        if enclosure.excludes_zero
            push!(rows,(gap_id=id,left_lambda_inv_fm=left,right_lambda_inv_fm=right,
                root_count=0,contour_count=0,count_method="sign_enclosure",
                margin_inv_fm=0.,endpoints_excluded=true,root_nodes_stable=true,
                roots_lambda_inv_fm="",max_root_shift_inv_fm=0.,passed=true,
                enclosure_lower=enclosure.lower,enclosure_upper=enclosure.upper,status="no_real_zero"))
            continue
        end
        # Account for the excluded strips; decrease margin only to RESOLVE a root,
        # never to omit a difficult boundary. Same residual/angle criteria retained.
        used_margin=min(margin,(right-left)/8)
        edge_ok=false
        for _ in 1:5
            left+used_margin<right-used_margin || break
            ea=inverse_enclosure(cells,K,left,left+used_margin)
            eb=inverse_enclosure(cells,K,right-used_margin,right)
            edge_ok=ea.excludes_zero && eb.excludes_zero
            edge_ok && break
            used_margin/=10
        end
        lo=R.certify_gap_roots((w,_)->inverse(w),q,[(left,right)];physical_sheet=true,
            real_axis=true,omega_nodes=root_nodes,endpoint_margin=used_margin)
        hi=R.certify_gap_roots((w,_)->inverse(w),q,[(left,right)];physical_sheet=true,
            real_axis=true,omega_nodes=2root_nodes,endpoint_margin=used_margin)
        stable=lo.passed && hi.passed && lo.count==hi.count
        shift=lo.count==hi.count ? maximum((abs(lo.roots[i].omega_inv_fm-hi.roots[i].omega_inv_fm)
            for i in eachindex(lo.roots));init=0.) : Inf
        stable &= shift<1e-7
        h=min(0.02,(right-left-2used_margin)/20)
        c=R.contour_count(reflection,left+used_margin,right-used_margin,-h,h;nodes=16)
        passed=edge_ok && stable && c.passed && c.count==hi.count
        push!(rows,(gap_id=id,left_lambda_inv_fm=left,right_lambda_inv_fm=right,
            root_count=hi.count,contour_count=c.count,count_method="bisection_and_argument_principle",
            margin_inv_fm=used_margin,endpoints_excluded=edge_ok,root_nodes_stable=stable,
            roots_lambda_inv_fm=join((r.omega_inv_fm for r in hi.roots),';'),
            max_root_shift_inv_fm=shift,passed=passed,enclosure_lower=enclosure.lower,
            enclosure_upper=enclosure.upper,status=passed ? "all_open_gap_roots_checked" : "unresolved_gap_or_endpoint"))
    end
    return (rows=rows,tail=tail,passed=tail.passed && !isempty(rows) && all(r.passed for r in rows),
        count=sum(r.root_count for r in rows),
        scope="all_real_gaps_of_compact_interpolant_not_cut_or_UHP_zeros")
end
end
