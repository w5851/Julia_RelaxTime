# Positive thermal-weighted interpolation and physical feature discovery.
module ThermalWeightedSampling

const AS = parentmodule(@__MODULE__)
using ...PairLoopQuadrature: PairLoopPlan, pair_b0_pv

export ThermalMeasure, thermal_measure, integrate_measure, WeightedSigmaTable,
       thermal_seed_grid, build_thermal_weighted_table, validated_table,
       real_axis_feature_seeds, scattering_feature_seeds

struct ThermalMeasure
    nodes::Vector{Float64}
    weights::Vector{Float64}
    events::Int
    function ThermalMeasure(nodes,weights,events=length(nodes))
        length(nodes)==length(weights) && !isempty(nodes) ||
            throw(ArgumentError("nonempty matched thermal nodes and weights required"))
        issorted(nodes) && all(isfinite,nodes) && all(diff(nodes).>0) ||
            throw(ArgumentError("thermal nodes must be finite and strictly increasing"))
        all(w->isfinite(w) && w>0,weights) ||
            throw(ArgumentError("thermal weights must be finite and positive"))
        new(nodes,weights,events)
    end
end

"""Exact discrete thermal measure, after combining identical angular arguments.

Weights contain the density normalization. Sorting/combining changes only
floating-point addition order. This is not a histogram or a bin-centre rule.
The measure is tied to the supplied background, outer grid and support bounds.
"""
function thermal_measure(plan,bounds;max_events=30_000_000,fold_equal_mass=true,
    relative_velocity=:cm)
    plan.prefactor>0 || throw(ArgumentError("positive thermal normalization required"))
    velocity_mode=AS._normalize_relative_velocity_policy(relative_velocity)
    groups=AS.angular_groups(plan.cos_relative,plan.phiw)
    nodes=Float64[]; weights=Float64[]
    sizehint!(nodes,500_000);sizehint!(weights,500_000)
    lo,hi=bounds; p,pw,ei,ej=plan.p,plan.pw,plan.ei,plan.ej
    folded=fold_equal_mass && plan.mi==plan.mj && ei==ej
    @inbounds for i in eachindex(p),j in (folded ? (i:length(p)) : eachindex(p))
        momentum=pw[i]*pw[j]*p[i]^2*p[j]^2*plan.prefactor
        reverse_momentum=folded && i!=j ? pw[j]*pw[i]*p[j]^2*p[i]^2*plan.prefactor : 0.0
        for group in groups
            angular=0.0;reverse_angular=0.0
            for (a,b) in group.pairs
                angular+=plan.cw[a]*plan.cw[b]*plan.fi[i,a]*plan.fj[j,b]
                if folded && i!=j
                    reverse_angular+=plan.cw[a]*plan.cw[b]*plan.fi[j,a]*plan.fj[i,b]
                end
            end
            combined=momentum*angular+reverse_momentum*reverse_angular
            combined==0 && continue
            for k in eachindex(group.cosines)
                s=plan.mi^2+plan.mj^2+2.0*(ei[i]*ej[j]-p[i]*p[j]*group.cosines[k])
                lo<s<hi || continue
                v=AS._relative_velocity_for_mode(velocity_mode,s,plan.mi,plan.mj,ei[i],ej[j])
                0<v<=2 || continue
                w=combined*group.weights[k]*v
                w==0 && continue
                isfinite(w) && w>0 || throw(DomainError(w,"nonpositive thermal measure"))
                length(nodes)<max_events || throw(ArgumentError("thermal measure event budget exhausted"))
                push!(nodes,s);push!(weights,w)
            end
        end
    end
    isempty(nodes) && throw(ArgumentError("thermal quadrature has no supported nonzero nodes"))
    events=length(nodes)
    order=sortperm(nodes;alg=QuickSort)
    sorted_nodes=nodes[order];sorted_weights=weights[order]
    count=0;correction=0.0
    @inbounds for i in eachindex(sorted_nodes)
        if count==0 || sorted_nodes[i]!=sorted_nodes[count]
            count+=1
            sorted_nodes[count]=sorted_nodes[i];sorted_weights[count]=sorted_weights[i]
            correction=0.0
        else
            term=sorted_weights[i]-correction
            next=sorted_weights[count]+term
            correction=(next-sorted_weights[count])-term
            sorted_weights[count]=next
        end
    end
    resize!(sorted_nodes,count);resize!(sorted_weights,count)
    return ThermalMeasure(sorted_nodes,sorted_weights,events)
end

function integrate_measure(measure::ThermalMeasure,query::F) where {F}
    total=0.0;correction=0.0
    @inbounds for i in eachindex(measure.nodes)
        term=measure.weights[i]*query(measure.nodes[i])-correction
        next=total+term
        correction=(next-total)-term;total=next
    end
    return total
end

struct WeightedSigmaTable
    nodes::Vector{Float64}
    sigma::Vector{Float64}
    s0::Float64
    weighted_end::Float64
end

@inline function interpolate_segment(a,b,fa,fb,s,s0,weighted)
    t=(s-a)/(b-a)
    return weighted ? ((1-t)*(a-s0)*fa+t*(b-s0)*fb)/(s-s0) : (1-t)*fa+t*fb
end

@inline function (table::WeightedSigmaTable)(s::Float64)
    first(table.nodes)<=s<=last(table.nodes) ||
        throw(DomainError(s,"table is bound to its validated thermal quadrature support"))
    i=searchsortedfirst(table.nodes,s)
    i==1 && return table.sigma[1]
    table.nodes[i]==s && return table.sigma[i]
    @inbounds return interpolate_segment(table.nodes[i-1],table.nodes[i],table.sigma[i-1],
        table.sigma[i],s,table.s0,table.nodes[i]<=table.weighted_end)
end

"""Generic seeds: threshold logarithm, square-root spacing, tail and branch points.

Feature locations are physical thresholds, never xi-specific exception lists.
The strategy does not use historical reference nodes or their failures.
"""
function thermal_seed_grid(measure::ThermalMeasure,s0::Float64;
    weighted_window=0.6,thresholds=Float64[],order=48)
    order>=8 || throw(ArgumentError("seed order must be at least eight"))
    a,b=first(measure.nodes),last(measure.nodes)
    s0<a<b || throw(ArgumentError("positive threshold offsets required"))
    xa,xb=a-s0,b-s0;split=clamp(s0+weighted_window,a,b)
    nodes=Float64[a,b,split]
    append!(nodes,s0 .+ 10.0.^range(log10(xa),log10(xb);length=order))
    append!(nodes,s0 .+ range(sqrt(xa),sqrt(xb);length=2order+1).^2)
    append!(nodes,range(split,b;length=2order+1))
    for threshold in thresholds
        a<threshold<b || continue
        push!(nodes,threshold)
        radius=min(threshold-a,b-threshold,1.0)
        tiny=max(128eps(threshold),radius*1e-8)
        if tiny<radius
            for offset in 10.0.^range(log10(tiny),log10(radius);length=16)
                a<threshold-offset<b && push!(nodes,threshold-offset)
                a<threshold+offset<b && push!(nodes,threshold+offset)
            end
        end
    end
    return sort!(unique(filter(s->a<=s<=b,nodes)))
end

struct SamplingBudgetExceeded <: Exception end

struct SamplingCell
    a::Float64
    b::Float64
    fa::Float64
    fb::Float64
    rate::Float64
    error::Float64
    split::Float64
end

function assess_cell(sample,measure,a,b,s0,weighted_end,guard_phase)
    fa,fb=sample(a),sample(b)
    weighted=b<=weighted_end
    first_i=searchsortedfirst(measure.nodes,a)
    a>first(measure.nodes) && (first_i=searchsortedlast(measure.nodes,a)+1)
    last_i=searchsortedlast(measure.nodes,b)
    first_i>last_i && return SamplingCell(a,b,fa,fb,0.0,0.0,(a+b)/2)
    # Quarter probes are reused after subdivision; irrationally shifted guard
    # probes and a supported thermal node provide independent locations.
    exact_leaf=last_i-first_i+1<=16
    probes=if exact_leaf
        # Once few supported nodes remain, evaluate all of them. The discrete
        # interpolation error is then known; off-support curvature need not
        # keep a conservative envelope refining forever.
        copy(measure.nodes[first_i:last_i])
    else
        fractions=Float64[0.25,0.5,0.75]
        guard_phase>=1 && push!(fractions,sqrt(2.0)-1)
        guard_phase>=2 && push!(fractions,pi-3)
        v=[a+(b-a)*f for f in fractions]
        push!(v,measure.nodes[(first_i+last_i)÷2])
        v
    end
    probes=sort!(unique(filter(s->a<s<b,probes)))
    values=sample.(probes)
    max_defect=0.0;split=(a+b)/2
    for (s,actual) in zip(probes,values)
        approximate=interpolate_segment(a,b,fa,fb,s,s0,weighted)
        defect=abs(actual-approximate)*(weighted ? s-s0 : 1.0)
        max_defect=max(max_defect,defect)
    end
    xs=vcat(a,probes,b);ys=vcat(fa,values,fb)
    coarse=0.0;defect_integral=0.0;effective_mass=0.0;k=1;largest_contribution=0.0
    @inbounds for i in first_i:last_i
        s,w=measure.nodes[i],measure.weights[i]
        while k<length(xs)-1 && s>xs[k+1]
            k+=1
        end
        old=interpolate_segment(a,b,fa,fb,s,s0,weighted)
        fine=interpolate_segment(xs[k],xs[k+1],ys[k],ys[k+1],s,s0,weighted)
        coarse+=w*old
        defect_integral+=w*abs(fine-old)
        if exact_leaf && a<s<b && w*abs(fine-old)>largest_contribution
            largest_contribution=w*abs(fine-old);split=s
        end
        effective_mass+=weighted ? w/(s-s0) : w
    end
    # This is an estimator, not a proof for an arbitrary unsampled function.
    # Both terms use the full discrete thermal weights and avoid cancellation.
    estimated=exact_leaf ? 2defect_integral : max(4defect_integral,2max_defect*effective_mass)
    all(isfinite,(coarse,defect_integral,effective_mass,estimated)) ||
        throw(DomainError((a,b,coarse,estimated),"nonfinite thermal interpolation integral"))
    return SamplingCell(a,b,fa,fb,coarse,estimated,split)
end

"""Refine the thermal-weighted interpolation defect, then run two shifted audits.

Every background uses the same rule. Failed audits enter the same refinement
loop, with all previous direct samples retained. Exhausted budgets return
`converged=false`, never a silently accepted cache. Convergence pertains to
cache error at this fixed outer quadrature; separate outer/t/B0 validation is
still necessary. No finite black-box sample set can prove absence of every
arbitrarily narrow peak; independent validation remains part of qualification.
"""
function build_thermal_weighted_table(f,measure::ThermalMeasure,s0::Float64;
    rtol=1e-3,atol=1e-12,weighted_window=0.6,thresholds=Float64[],seed_order=48,
    extra_seeds=Float64[],max_evaluations=40_000,max_nodes=12_000,max_rounds=40)
    isfinite(rtol) && isfinite(atol) && rtol>=0 && atol>=0 && rtol+atol>0 ||
        throw(ArgumentError("finite nonnegative tolerances with one positive required"))
    max_evaluations>0 && max_nodes>=2 && max_rounds>0 ||
        throw(ArgumentError("positive sampling budgets required"))
    isfinite(weighted_window) && weighted_window>=0 ||
        throw(ArgumentError("nonnegative finite weighted window required"))
    values=Dict{Float64,Float64}()
    function sample(s)
        get!(values,s) do
            length(values)<max_evaluations || throw(SamplingBudgetExceeded())
            v=Float64(f(s))
            isfinite(v) && v>=0 || throw(DomainError((s,v),"invalid direct cross section"))
            v
        end
    end
    seeds=thermal_seed_grid(measure,s0;weighted_window=weighted_window,
        thresholds=thresholds,order=seed_order)
    append!(seeds,filter(s->first(measure.nodes)<s<last(measure.nodes),extra_seeds))
    sort!(unique!(seeds))
    split=s0+weighted_window
    cells=SamplingCell[];history=NamedTuple[];phase=0;converged=false;reason="round_budget"
    try
        length(seeds)<=max_nodes || throw(SamplingBudgetExceeded())
        cells=[assess_cell(sample,measure,a,b,s0,split,phase) for (a,b) in zip(seeds[1:end-1],seeds[2:end])]
        for round in 1:max_rounds
            rate=sum(c.rate for c in cells);error=sum(c.error for c in cells)
            isfinite(rate) && isfinite(error) ||
                throw(DomainError((rate,error),"nonfinite thermal error total"))
            tolerance=atol+rtol*abs(rate)
            push!(history,(round=round,guard_phase=phase,nodes=length(cells)+1,
                evaluations=length(values),rate=rate,estimated_error=error,tolerance=tolerance))
            if error<=tolerance
                if phase==2
                    converged=true;reason="thermal_defect_and_two_shifted_audits";break
                end
                phase+=1
                cells=[assess_cell(sample,measure,c.a,c.b,s0,split,phase) for c in cells]
                continue
            end
            length(cells)+1<max_nodes || (reason="node_budget";break)
            # Refine the cells responsible for at least half of the current
            # estimated error. Priority is the contribution to the thermal rate.
            ranked=sortperm([c.error for c in cells];rev=true)
            chosen=Set{Int}();accounted=0.0
            for i in ranked
                length(cells)+length(chosen)+1>=max_nodes && break
                push!(chosen,i);accounted+=cells[i].error
                accounted>=0.5error && break
            end
            isempty(chosen) && (reason="node_budget";break)
            refined=SamplingCell[]
            for (i,c) in enumerate(cells)
                if i in chosen && c.a<c.split<c.b
                    push!(refined,assess_cell(sample,measure,c.a,c.split,s0,split,phase))
                    push!(refined,assess_cell(sample,measure,c.split,c.b,s0,split,phase))
                else
                    push!(refined,c)
                end
            end
            cells=refined
        end
    catch err
        err isa SamplingBudgetExceeded || rethrow()
        reason="evaluation_or_seed_budget"
    end
    nodes=isempty(cells) ? Float64[] : vcat(first(cells).a,[c.b for c in cells])
    sigma=isempty(cells) ? Float64[] : vcat(first(cells).fa,[c.fb for c in cells])
    table=WeightedSigmaTable(nodes,sigma,s0,split)
    return (table=table,converged=converged,reason=reason,guard_phase=phase,
        rate=isempty(cells) ? NaN : sum(c.rate for c in cells),
        estimated_error=isempty(cells) ? Inf : sum(c.error for c in cells),
        evaluations=length(values),samples=values,history=history)
end

function validated_table(result)
    result.converged || throw(ArgumentError("thermal cache did not converge: $(result.reason)"))
    return result.table
end

"""Locate real-axis denominator crossings and seed their estimated width.

These are numerical feature anchors, not complex poles or Mott diagnoses.
Sign-change scanning is combined with the independent thermal interpolation
audits; it is not an existence theorem for arbitrary unresolved roots.
"""
function _regular_feature_seeds(denominator,grid)
    samples=denominator.(grid)
    invalid=findall(z->!isfinite(real(z)) || !isfinite(imag(z)),samples)
    isempty(invalid) || throw(DomainError([(s=grid[i],value=samples[i]) for i in invalid],
        "nonfinite propagator feature scan"))
    roots=NamedTuple[];seeds=Float64[]
    for i in 1:length(grid)-1
        a,b=grid[i],grid[i+1];fa,fb=real(samples[i]),real(samples[i+1])
        (fa==0 || fb==0 || signbit(fa)!=signbit(fb)) || continue
        if fa==0
            b=a
        elseif fb==0
            a=b
        else
            for _ in 1:48
                mid=(a+b)/2
                a<mid<b || break
                fm=real(denominator(mid))
                if signbit(fm)==signbit(fa)
                    a=mid;fa=fm
                else
                    b=mid
                end
            end
        end
        root=(a+b)/2
        h=max(256eps(root),1e-6*max(1.0,abs(root)))
        left=max(first(grid),root-h);right=min(last(grid),root+h)
        slope=(real(denominator(right))-real(denominator(left)))/(right-left)
        value=denominator(root)
        width=abs(slope)>0 ? max(abs(imag(value)),1e-12)/abs(slope) : grid[i+1]-grid[i]
        width=max(width,128eps(root))
        push!(roots,(location=root,width=width,real_denominator=real(value),imag_denominator=imag(value)))
        push!(seeds,root)
        for factor in (-16.0,-8.0,-4.0,-2.0,-1.0,-0.5,-0.25,0.25,0.5,1.0,2.0,4.0,8.0,16.0)
            s=root+factor*width
            first(grid)<s<last(grid) && push!(seeds,s)
        end
    end
    return (seeds=sort!(unique(seeds)),features=roots)
end

function real_axis_feature_seeds(denominator,grid;breakpoints=Float64[])
    length(grid)>=2 && issorted(grid) && all(isfinite,grid) ||
        throw(ArgumentError("feature grid must contain ordered finite points"))
    lower,upper=first(grid),last(grid)
    lower<upper || throw(ArgumentError("feature grid must span a nonempty interval"))
    all(isfinite,breakpoints) || throw(ArgumentError("feature breakpoints must be finite"))
    cuts=sort!(unique(filter(s->lower<=s<=upper,Float64.(breakpoints))))
    isempty(cuts) && return _regular_feature_seeds(denominator,grid)
    # A sharp regulator endpoint is a known singularity, not a denominator
    # zero. Keep both one-sided limits and never bisect across the singularity.
    # The margin protects the subsequent s -> q^2 calculation from rounding
    # onto the endpoint; it does not replace any cross-section evaluation.
    margin(c)=512eps(max(abs(c),1.0))
    probes=Float64[]
    for cut in cuts, exponent in -12:-2, side in (-1.0,1.0)
        s=cut+side*max(margin(cut),max(1.0,abs(cut))*10.0^exponent)
        lower<=s<=upper && all(c->abs(s-c)>=margin(c),cuts) && push!(probes,s)
    end
    augmented=sort!(unique(vcat(Float64.(grid),probes)))
    filter!(s->all(c->abs(s-c)>=margin(c),cuts),augmented)
    edges=sort!(unique(vcat(lower,cuts,upper)))
    roots=NamedTuple[];seeds=copy(probes)
    for (lo,hi) in zip(edges[1:end-1],edges[2:end])
        local_grid=filter(s->lo<=s<=hi,augmented)
        length(local_grid)>=2 || continue
        found=_regular_feature_seeds(denominator,local_grid)
        append!(seeds,found.seeds);append!(roots,found.features)
    end
    return (seeds=sort!(unique(seeds)),features=roots)
end

function scattering_feature_seeds(process,state,measure::ThermalMeasure,s0)
    rt=parentmodule(AS);qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    channels=rt.TotalPropagator.SCATTERING_MESON_MAP[process][:channels]
    haskey(channels,:s) || return (seeds=Float64[],features=NamedTuple[])
    info=channels[:s]
    thresholds=[4qp.m.u^2,4qp.m.d^2,4qp.m.s^2,(qp.m.u+qp.m.s)^2]
    grid=thermal_seed_grid(measure,s0;thresholds=thresholds,order=64)
    append!(grid,range(sqrt(first(grid)),sqrt(last(grid));length=193).^2)
    sort!(unique!(grid))
    # Do not insert discovery probes into the quantized production cache:
    # discovery must not change the subsequent cross-section evaluation order.
    # Equal-mass k=0 discovery uses the already derived endpoint-subtracted PV
    # representation. The historical kernel is nonfinite at some internal ss
    # threshold probes. This backend only proposes anchors; sigma and the
    # independent numerical comparisons retain the production kernel.
    pair_plans=Dict((m,mu)=>PairLoopPlan(m,mu,tp;n_nodes=128)
        for (m,mu) in ((qp.m.u,qp.μ.u),(qp.m.s,qp.μ.s)))
    function polarization(channel,m1,m2,mu1,mu2,a1,a2,ns,s)
        if m1==m2 && mu1==mu2
            b0=pair_b0_pv(pair_plans[(m1,mu1)],s)
            coefficient=channel===:P ? -s : 4m1^2-s
            return -Main.Constants_PNJL.N_color/(8pi^2)*(a1+a2+coefficient*b0)
        end
        return complex(rt.PolarizationAniso.polarization_aniso(channel,sqrt(s),0.0,
            m1,m2,mu1,mu2,tp.T,tp.Φ,tp.Φbar,tp.ξ,a1,a2,ns)...)
    end
    seeds=Float64[];features=NamedTuple[]
    function collect_family(family,denominator,breakpoints)
        found=real_axis_feature_seeds(denominator,grid;breakpoints=breakpoints)
        append!(seeds,found.seeds)
        append!(features,[merge((family=String(family),),r) for r in found.features])
    end
    for meson in info[:simple]
        pp=rt.TotalPropagator.get_polarization_params(meson,qp)
        coupling=meson===:pi ? kc.K123_plus : meson===:K ? kc.K4567_plus :
            meson===:sigma_pi ? kc.K123_minus : kc.K4567_minus
        denominator=s->1-4coupling*polarization(pp.channel,pp.m1,pp.m2,pp.μ1,pp.μ2,
            pp.A1,pp.A2,pp.num_s_quark,s)
        cutoff=Main.Constants_PNJL.Λ_inv_fm
        e1,e2=sqrt(cutoff^2+pp.m1^2),sqrt(cutoff^2+pp.m2^2)
        delta_mu=pp.μ1-pp.μ2
        # Include the sum/difference endpoints and the symmetrized K loop.
        cuts=[(energy+shift)^2 for energy in (e1+e2,e1-e2) for shift in (-delta_mu,delta_mu)]
        collect_family(meson,denominator,cuts)
    end
    for channel in (:P,:S)
        info[channel===:P ? :mixed_P : :mixed_S] || continue
        function denominator(s)
            pu=polarization(channel,qp.m.u,qp.m.u,qp.μ.u,qp.μ.u,qp.A.u,qp.A.u,0,s)
            ps=polarization(channel,qp.m.s,qp.m.s,qp.μ.s,qp.μ.s,qp.A.s,qp.A.s,2,s)
            mm=rt.EffectiveCouplings.mixing_matrix_elements(pu,ps,kc,channel)
            return (mm.M00*mm.M88-mm.M08^2)/(channel===:P ? kc.det_K_plus : kc.det_K_minus)
        end
        cutoff=Main.Constants_PNJL.Λ_inv_fm
        collect_family(Symbol(:mixed_,channel),denominator,
            [4(cutoff^2+qp.m.u^2),4(cutoff^2+qp.m.s^2)])
    end
    return (seeds=sort!(unique(seeds)),features=features)
end

end
