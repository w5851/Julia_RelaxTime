# Fixed-background performance candidates; no production default is changed here.
if !isdefined(Main, :ThresholdUpgradeCommon)
    include(joinpath(@__DIR__, "ThresholdUpgradeCommon.jl"))
end

module ThermalAverageCandidates

using Main.ThresholdUpgradeCommon: AS, CS

export ThermalPlan, prepare_thermal_plan, prepared_rate, shared_rates, integrate_queries,
       project_linear_rates

struct ThermalPlan
    mi::Float64
    mj::Float64
    p::Vector{Float64}
    pw::Vector{Float64}
    ei::Vector{Float64}
    ej::Vector{Float64}
    fi::Matrix{Float64}
    fj::Matrix{Float64}
    cw::Vector{Float64}
    phiw::Vector{Float64}
    cos_relative::Array{Float64,3}
    prefactor::Float64
end

function prepare_thermal_plan(process, state; p_nodes=128, angle_nodes=20, phi_nodes=36,
    p_grid=nothing, p_w=nothing, cos_grid=nothing, cos_w=nothing,
    phi_grid=nothing, phi_w=nothing, density_p_grid=nothing, density_p_w=nothing)
    qp, tp = state.quark_params, state.thermo_params
    i, j, _, _ = CS.parse_particles_from_process(process)
    mi, mj = AS.get_mass(i, qp), AS.get_mass(j, qp)
    mui, muj = AS.get_mu(i, qp), AS.get_mu(j, qp)
    p, pw = p_grid === nothing ?
        AS._cached_semi_infinite_momentum_grid(p_nodes, AS.DEFAULT_SEMI_INF_SCALE) : (p_grid, p_w)
    ct, cw = cos_grid === nothing ? AS._cached_interval_grid(-1.0, 1.0, angle_nodes) : (cos_grid, cos_w)
    az, aw = phi_grid === nothing ? AS._cached_interval_grid(0.0, 2pi, phi_nodes) : (phi_grid, phi_w)
    ei, ej = AS.energy_from_p.(p, mi), AS.energy_from_p.(p, mj)
    fi = [AS.distribution_with_anisotropy(i, v, mi, mui, tp.T, tp.Φ, tp.Φbar, tp.ξ, c) for v in p, c in ct]
    fj = [AS.distribution_with_anisotropy(j, v, mj, muj, tp.T, tp.Φ, tp.Φbar, tp.ξ, c) for v in p, c in ct]
    st = sqrt.(max.(1.0 .- ct .* ct, 0.0))
    ca = cos.(az)
    cr = [ct[a] * ct[b] + st[a] * st[b] * ca[k] for k in eachindex(az), a in eachindex(ct), b in eachindex(ct)]
    rhoi = AS.number_density(i, mi, mui, tp.T, tp.Φ, tp.Φbar, tp.ξ;
        p_nodes=AS.DEFAULT_SEMI_INF_NODES, angle_nodes=angle_nodes, p_grid=density_p_grid, p_w=density_p_w)
    rhoj = AS.number_density(j, mj, muj, tp.T, tp.Φ, tp.Φbar, tp.ξ;
        p_nodes=AS.DEFAULT_SEMI_INF_NODES, angle_nodes=angle_nodes, p_grid=density_p_grid, p_w=density_p_w)
    prefactor = rhoi == 0 || rhoj == 0 ? 0.0 : AS.DQ^2 / (32pi^5 * rhoi * rhoj)
    return ThermalPlan(mi, mj, p, pw, ei, ej, fi, fj, cw, aw, cr, prefactor)
end

@inline function relative_velocity(s, mi, mj, ::Val{:legacy})
    root = sqrt(s)
    ei = (s + mi^2 - mj^2) / (2root)
    ej = (s - mi^2 + mj^2) / (2root)
    pcm = sqrt(max(0.0, (s - (mi + mj)^2) * (s - (mi - mj)^2))) / (2root)
    numerator = (ei * ej + pcm * pcm)^2 - (mi * mj)^2
    return numerator > 0 ? sqrt(numerator) / (ei * ej) : 0.0
end

@inline function relative_velocity(s, mi, mj, ::Val{:equal_mass})
    mi == mj || return relative_velocity(s, mi, mj, Val(:legacy))
    return 2sqrt(max(0.0, (s - 4mi^2) / s))
end

@inline accumulate_queries(::Tuple{}, ::Tuple{}, s, base, ::Tuple{}) = ()

@inline function accumulate_queries(queries::Tuple{F,Vararg{Any,N}}, totals::Tuple,
    s, base, bounds::Tuple) where {F,N}
    # Heterogeneous callable tuples must be unrolled by type. Runtime indexing
    # inside ntuple can lose inference with four different query wrappers.
    lo, hi = first(bounds)
    value = first(totals)
    if lo < s < hi
        value += base * first(queries)(s)
    end
    return (value, accumulate_queries(Base.tail(queries), Base.tail(totals),
        s, base, Base.tail(bounds))...)
end

function integrate_queries(plan::ThermalPlan, queries::Tuple, bounds::Tuple, velocity::Val)
    length(queries) == length(bounds) || throw(ArgumentError("query/bound lengths differ"))
    totals = map(_ -> 0.0, queries)
    plan.prefactor == 0 && return totals
    lower = minimum(first, bounds)
    upper = maximum(last, bounds)
    mi, mj = plan.mi, plan.mj
    p, pw, ei, ej, fi, fj = plan.p, plan.pw, plan.ei, plan.ej, plan.fi, plan.fj
    cw, aw, cr = plan.cw, plan.phiw, plan.cos_relative
    @inbounds for i in eachindex(p), j in eachindex(p)
        pi, pj = p[i], p[j]
        energy_product, momentum_product = ei[i] * ej[j], pi * pj
        for a in eachindex(cw), b in eachindex(cw)
            f1, f2 = fi[i, a], fj[j, b]
            (f1 == 0 || f2 == 0) && continue
            for k in eachindex(aw)
                s = mi^2 + mj^2 + 2 * (energy_product - momentum_product * cr[k, a, b])
                (s <= lower || s >= upper) && continue
                v = relative_velocity(s, mi, mj, velocity)
                (v == 0 || v > 2) && continue
                base = pw[i] * pw[j] * cw[a] * cw[b] * aw[k] * pi^2 * pj^2 * f1 * f2 * v
                totals = accumulate_queries(queries, totals, s, base, bounds)
            end
        end
    end
    return map(value -> plan.prefactor * value, totals)
end

function process_bounds(process, state; sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
    p_grid=nothing, p_w=nothing, apply_s_domain_cut=true)
    qp = state.quark_params
    i, j, c, d = CS.parse_particles_from_process(process)
    mi, mj, mc, md = (AS.get_mass(q, qp) for q in (i, j, c, d))
    cutoff, lo, hi = AS._resolve_sigma_support_bounds(mi, mj, mc, md, sigma_cutoff, p_grid, p_w)
    return apply_s_domain_cut && !isnan(cutoff) ? (lo, hi) : ((mi + mj)^2, Inf)
end

function shared_rates(queries::Tuple, processes::Tuple, state; velocity=:equal_mass,
    sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm, apply_s_domain_cut=true, kwargs...)
    velocity in (:legacy, :equal_mass) || throw(ArgumentError("unknown velocity candidate"))
    length(queries) == length(processes) || throw(ArgumentError("query/process lengths differ"))
    isempty(queries) && throw(ArgumentError("at least one query is required"))
    incoming = CS.parse_particles_from_process(first(processes))[1:2]
    all(p -> CS.parse_particles_from_process(p)[1:2] == incoming, processes) ||
        throw(ArgumentError("shared quadrature requires identical incoming particle symbols"))
    plan = prepare_thermal_plan(first(processes), state; kwargs...)
    p_grid, p_w = get(kwargs, :p_grid, nothing), get(kwargs, :p_w, nothing)
    bounds = map(p -> process_bounds(p, state; sigma_cutoff=sigma_cutoff,
        apply_s_domain_cut=apply_s_domain_cut, p_grid=p_grid, p_w=p_w), processes)
    return integrate_queries(plan, queries, bounds, Val(velocity))
end

prepared_rate(query, process, state; kwargs...) = only(shared_rates((query,), (process,), state; kwargs...))

"""Integrate several piecewise linear/raw or x-sigma tables on their union grid.

The four nonnegative barycentric sums are exact for these interpolants at the
unchanged 5D nodes. No bin-center replacement or additional s interpolation is
introduced. The primary callable is evaluated normally, including addback and
clipping. weighted_until[k] marks the end of each query's x-sigma segment.
"""
function project_linear_rates(plan::ThermalPlan, knots::Vector{Float64}, s0::Float64,
    primary::F, queries::Tuple, weighted_until::Tuple, bounds::Tuple{Float64,Float64}) where {F}
    length(queries)==length(weighted_until) || throw(ArgumentError("representation lengths differ"))
    issorted(knots) && all(diff(knots).>0) || throw(ArgumentError("strictly increasing union grid required"))
    left=zeros(length(knots)-1); right=similar(left); fill!(right,0.0)
    inv_left=copy(left); inv_right=copy(left)
    fallback=map(_->0.0,queries); original=0.0
    mi,mj=plan.mi,plan.mj
    p,pw,ei,ej,fi,fj=plan.p,plan.pw,plan.ei,plan.ej,plan.fi,plan.fj
    cw,aw,cr=plan.cw,plan.phiw,plan.cos_relative
    lower,upper=bounds
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
                original+=base*primary(s)
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
