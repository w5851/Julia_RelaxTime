# Analysis-only exact grouping of identical angular geometry.
if !isdefined(Main, :ThermalAverageCandidates)
    include(joinpath(@__DIR__, "ThermalAverageCandidates.jl"))
end

module GroupedThermalProjection

using Main.ThermalAverageCandidates: ThermalPlan, relative_velocity
using Main.ThresholdUpgradeCommon: AS

export AngularGroup, angular_groups, grouped_project_linear_rates, compensated_primary_rate

const AngularGroup = AS.ThermalAngularQuadrature.AngularGroup

"""Group only bitwise-equal relative-angle vectors, then equal phi cosines.

No symmetry is assumed for occupations or custom grids. Every polar pair and
every azimuthal weight remains in its group. This changes floating-point
addition order, but neither quadrature nodes nor function arguments.
"""
function angular_groups(plan::ThermalPlan)
    return AS.ThermalAngularQuadrature.angular_groups(plan.cos_relative,plan.phiw)
end

"""Independent full-node primary integral with Kahan summation.

Use this to distinguish accumulated roundoff in a long positive sum from
changes in quadrature or angular grouping.
"""
function compensated_primary_rate(plan::ThermalPlan, query::F, bounds) where {F}
    mi,mj=plan.mi,plan.mj
    p,pw,ei,ej,fi,fj=plan.p,plan.pw,plan.ei,plan.ej,plan.fi,plan.fj
    cw,aw,cr=plan.cw,plan.phiw,plan.cos_relative
    lower,upper=bounds; total=0.0; correction=0.0
    @inbounds for i in eachindex(p),j in eachindex(p)
        pi,pj=p[i],p[j]
        for a in eachindex(cw),b in eachindex(cw)
            f1,f2=fi[i,a],fj[j,b]
            (f1==0 || f2==0) && continue
            for k in eachindex(aw)
                s=mi^2+mj^2+2*(ei[i]*ej[j]-pi*pj*cr[k,a,b])
                (s<=lower || s>=upper) && continue
                v=relative_velocity(s,mi,mj,Val(:equal_mass))
                (v==0 || v>2) && continue
                base=pw[i]*pw[j]*cw[a]*cw[b]*aw[k]*pi^2*pj^2*f1*f2*v
                term=base*query(s)-correction
                next=total+term
                correction=(next-total)-term
                total=next
            end
        end
    end
    return plan.prefactor*total
end

"""Same union-grid projection with exactly repeated angular arguments grouped."""
function grouped_project_linear_rates(plan::ThermalPlan, knots::Vector{Float64}, s0::Float64,
    primary::F, queries::Tuple, weighted_until::Tuple, bounds::Tuple{Float64,Float64};
    groups=angular_groups(plan)) where {F}
    length(queries)==length(weighted_until) || throw(ArgumentError("representation lengths differ"))
    issorted(knots) && all(diff(knots).>0) || throw(ArgumentError("strictly increasing union grid required"))
    left=zeros(length(knots)-1); right=copy(left)
    inv_left=copy(left); inv_right=copy(left)
    fallback=map(_->0.0,queries); original=0.0; original_correction=0.0
    mi,mj=plan.mi,plan.mj
    p,pw,ei,ej,fi,fj=plan.p,plan.pw,plan.ei,plan.ej,plan.fi,plan.fj
    cw=plan.cw; lower,upper=bounds
    @inbounds for i in eachindex(p),j in eachindex(p)
        pi,pj=p[i],p[j]
        momentum_weight=pw[i]*pw[j]*pi^2*pj^2
        for group in groups
            angular_weight=0.0
            for (a,b) in group.pairs
                angular_weight+=cw[a]*cw[b]*fi[i,a]*fj[j,b]
            end
            angular_weight==0.0 && continue
            for k in eachindex(group.cosines)
                s=mi^2+mj^2+2*(ei[i]*ej[j]-pi*pj*group.cosines[k])
                (s<=lower || s>=upper) && continue
                v=relative_velocity(s,mi,mj,Val(:equal_mass))
                (v==0 || v>2) && continue
                base=momentum_weight*angular_weight*group.weights[k]*v
                primary_term=base*primary(s)-original_correction
                primary_next=original+primary_term
                original_correction=(primary_next-original)-primary_term
                original=primary_next
                if s<knots[1] || s>knots[end]
                    fallback=map((sum,q)->sum+base*q(s),fallback,queries)
                    continue
                end
                cell=max(1,searchsortedfirst(knots,s)-1)
                aa,bb=knots[cell],knots[cell+1]
                t=(s-aa)/(bb-aa)
                wl=base*(1-t); wr=base*t
                left[cell]+=wl; right[cell]+=wr
                inv_left[cell]+=wl/(s-s0); inv_right[cell]+=wr/(s-s0)
            end
        end
    end
    values=map(queries,weighted_until,fallback) do query,weighted_end,remainder
        y=query.(knots)
        total=remainder
        for k in eachindex(left)
            if (knots[k]+knots[k+1])/2<=weighted_end
                total+=inv_left[k]*(knots[k]-s0)*y[k]+inv_right[k]*(knots[k+1]-s0)*y[k+1]
            else
                total+=left[k]*y[k]+right[k]*y[k+1]
            end
        end
        plan.prefactor*total
    end
    return (plan.prefactor*original,values...)
end

end
