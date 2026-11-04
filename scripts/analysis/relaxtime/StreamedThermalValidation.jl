# Memory-bounded checks on outer grids finer than a retained discrete measure.
if !isdefined(Main, :ThermalAverageCandidates)
    include(joinpath(@__DIR__, "ThermalAverageCandidates.jl"))
end

module StreamedThermalValidation

using Main.ThermalAverageCandidates: ThermalPlan, relative_velocity
using Main.ThresholdUpgradeCommon: AS

export streamed_integrals, thermal_support, TableComparison
export compile_comparison

struct TableComparison{C,R,P}
    candidate::C
    reference::R
    previous::P
end

@inline function (comparison::TableComparison)(s)
    candidate_value=comparison.candidate(s)
    reference_value=comparison.reference(s)
    previous_value=comparison.previous(s)
    return (candidate_value,reference_value,abs(candidate_value-reference_value),
        abs(previous_value-reference_value))
end

struct CompiledTableComparison
    nodes::Vector{Float64}
    values::Vector{NTuple{3,Float64}}
    left::Vector{NTuple{3,Float64}}
    right::Vector{NTuple{3,Float64}}
    weighted::Vector{NTuple{3,Bool}}
    s0::NTuple{3,Float64}
end

"""Use one interval lookup for three piecewise linear/x-sigma tables.

The union contains every knot of every table. It changes only floating-point
evaluation order, without binning thermal nodes or approximating their weights.
"""
function compile_comparison(comparison::TableComparison)
    tables=(comparison.candidate,comparison.reference,comparison.previous)
    all(t->t isa AS.ThermalWeightedSampling.WeightedSigmaTable,tables) ||
        throw(ArgumentError("compiled comparison requires WeightedSigmaTable inputs"))
    lo=maximum(t->first(t.nodes),tables);hi=minimum(t->last(t.nodes),tables)
    lo<hi || throw(ArgumentError("comparison tables have no common interval"))
    nodes=sort!(unique(vcat([lo,hi],tables[1].nodes,tables[2].nodes,tables[3].nodes)))
    filter!(s->lo<=s<=hi,nodes)
    values=[map(t->t(s),tables) for s in nodes]
    left=NTuple{3,Float64}[];right=NTuple{3,Float64}[];weighted=NTuple{3,Bool}[]
    for i in 1:length(nodes)-1
        flags=map(t->t.nodes[searchsortedfirst(t.nodes,nodes[i+1])]<=t.weighted_end,tables)
        push!(weighted,flags)
        push!(left,ntuple(k->values[i][k]*(flags[k] ? nodes[i]-tables[k].s0 : 1.0),3))
        push!(right,ntuple(k->values[i+1][k]*(flags[k] ? nodes[i+1]-tables[k].s0 : 1.0),3))
    end
    return CompiledTableComparison(nodes,values,left,right,weighted,map(t->t.s0,tables))
end

@inline function (query::CompiledTableComparison)(s)
    first(query.nodes)<=s<=last(query.nodes) || throw(DomainError(s,"compiled table support exceeded"))
    i=searchsortedfirst(query.nodes,s)
    values=if query.nodes[i]==s
        @inbounds query.values[i]
    else
        @inbounds t=(s-query.nodes[i-1])/(query.nodes[i]-query.nodes[i-1])
        @inbounds ntuple(k->((1-t)*query.left[i-1][k]+t*query.right[i-1][k])/
            (query.weighted[i-1][k] ? s-query.s0[k] : 1.0),Val(3))
    end
    c,r,p=values
    return (c,r,abs(c-r),abs(p-r))
end

mutable struct IntegralAccumulator{F,N}
    query::F
    total::NTuple{N,Float64}
    correction::NTuple{N,Float64}
end

@inline function visit!(acc::IntegralAccumulator{F,N}, s, w) where {F,N}
    values = acc.query(s)
    terms = ntuple(i -> w * values[i] - acc.correction[i], Val(N))
    next = ntuple(i -> acc.total[i] + terms[i], Val(N))
    acc.correction = ntuple(i -> (next[i] - acc.total[i]) - terms[i], Val(N))
    acc.total = next
    return nothing
end

mutable struct SupportAccumulator
    lower::Float64
    upper::Float64
end

@inline function visit!(acc::SupportAccumulator, s, w)
    acc.lower = min(acc.lower, s)
    acc.upper = max(acc.upper, s)
    return nothing
end

function visit_thermal!(acc, plan::ThermalPlan, bounds; fold_equal_mass=true)
    plan.prefactor > 0 || throw(ArgumentError("positive thermal normalization required"))
    groups = AS.angular_groups(plan.cos_relative, plan.phiw)
    lo, hi = bounds
    p, pw, ei, ej = plan.p, plan.pw, plan.ei, plan.ej
    folded = fold_equal_mass && plan.mi == plan.mj && ei == ej
    events = 0
    @inbounds for i in eachindex(p), j in (folded ? (i:length(p)) : eachindex(p))
        momentum = pw[i] * pw[j] * p[i]^2 * p[j]^2 * plan.prefactor
        reverse_momentum = folded && i != j ?
            pw[j] * pw[i] * p[j]^2 * p[i]^2 * plan.prefactor : 0.0
        for group in groups
            angular = 0.0
            reverse_angular = 0.0
            for (a, b) in group.pairs
                angular += plan.cw[a] * plan.cw[b] * plan.fi[i,a] * plan.fj[j,b]
                if folded && i != j
                    reverse_angular += plan.cw[a] * plan.cw[b] * plan.fi[j,a] * plan.fj[i,b]
                end
            end
            combined = momentum * angular + reverse_momentum * reverse_angular
            combined == 0 && continue
            for k in eachindex(group.cosines)
                s = plan.mi^2 + plan.mj^2 +
                    2.0 * (ei[i] * ej[j] - p[i] * p[j] * group.cosines[k])
                lo < s < hi || continue
                v = relative_velocity(s, plan.mi, plan.mj, Val(:equal_mass))
                0 < v <= 2 || continue
                w = combined * group.weights[k] * v
                w == 0 && continue
                isfinite(w) && w > 0 || throw(DomainError(w, "invalid thermal weight"))
                visit!(acc, s, w)
                events += 1
            end
        end
    end
    events > 0 || throw(ArgumentError("thermal quadrature has no supported positive weight"))
    return events
end

"""Integrate a tuple-valued query at every supported outer quadrature node.

No histogram, interpolation of the thermal measure, or stored event array is
used. The caller can integrate two tables and their absolute difference in
one pass, thereby checking cache reuse on the *new* outer grid. This does not
grant a table validity outside the support that was actually checked.
"""
function streamed_integrals(plan::ThermalPlan, query::F, bounds, ::Val{N};
    fold_equal_mass=true) where {F,N}
    N > 0 || throw(ArgumentError("at least one integral required"))
    zero = ntuple(_ -> 0.0, Val(N))
    acc = IntegralAccumulator{F,N}(query, zero, zero)
    events = visit_thermal!(acc, plan, bounds; fold_equal_mass=fold_equal_mass)
    all(isfinite, acc.total) || throw(DomainError(acc.total, "nonfinite thermal integral"))
    return (values=acc.total, events=events)
end

"""Find the actual positive-weight support without allocating an event array."""
function thermal_support(plan::ThermalPlan, bounds; fold_equal_mass=true)
    acc = SupportAccumulator(Inf, -Inf)
    events = visit_thermal!(acc, plan, bounds; fold_equal_mass=fold_equal_mass)
    return (lower=acc.lower, upper=acc.upper, events=events)
end

end
