# Isolated, single-process comparisons of a loop backend through the unchanged
# production amplitude and t quadrature. This is not a production backend API.
if !isdefined(Main, :ThresholdUpgradeCommon)
    include(joinpath(@__DIR__, "ThresholdUpgradeCommon.jl"))
end

module CrossSectionLoopProbe

using Main.ThresholdUpgradeCommon: AS, CS
const RT=Main.RelaxTime
const PC=RT.PolarizationCache
const TP=RT.TotalPropagator
const PolValue=Tuple{Float64,Float64}

export LoopProbe, probe_cross_section, polarization_specs

struct LoopProbe{F}
    backend::F
    channels::Symbol
    cache::Dict{PC.PolarizationKey,PolValue}
end

function LoopProbe(backend;channels=:all)
    channels in (:all,:s,:finite_k) || throw(ArgumentError("channels must be :all, :s or :finite_k"))
    return LoopProbe(backend,channels,Dict{PC.PolarizationKey,PolValue}())
end

function polarization_specs(info, qp)
    specs=[TP.get_polarization_params(meson,qp) for meson in info[:simple]]
    for channel in (:P,:S)
        info[channel===:P ? :mixed_P : :mixed_S] || continue
        for (flavor,ns) in ((:u,0),(:s,2))
            m=getproperty(qp.m,flavor);mu=getproperty(qp.μ,flavor);a=getproperty(qp.A,flavor)
            push!(specs,(channel=channel,m1=m,m2=m,μ1=mu,μ2=mu,A1=a,A2=a,num_s_quark=ns))
        end
    end
    return unique(specs)
end

"""Evaluate a selected backend, restoring every overwritten cache key on exit.

The backend takes the same arguments as polarization_aniso and returns its
real/imaginary pair. Its private cache belongs to one explicit backend object;
different backends and tolerances must use different objects. With :all, every
production polarization query must hit a supplied key or this probe fails.
"""
function probe_cross_section(probe::LoopProbe,process,s,state;n_points=24)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    ctx=CS.prepare_scattering_context(process,qp)
    masses=(ctx.m1,ctx.m2,ctx.m3,ctx.m4)
    s>max((ctx.m1+ctx.m2)^2,(ctx.m3+ctx.m4)^2) || return 0.0
    tb=CS.calculate_t_bounds(s,masses...)
    tlo=ctx.particle_c==ctx.particle_d ? tb.t_min/2 : tb.t_min
    thi=tb.t_max
    thi-tlo>CS.EPS_KINEMATIC || return 0.0
    z,_=RT.GaussLegendre.standard_nodes_weights(n_points)
    ts=(thi-tlo)/2 .* z .+ (thi+tlo)/2
    previous=Dict{PC.PolarizationKey,Union{Nothing,PolValue}}()
    # Bound diagnostic memory without changing the backend or tolerance.
    length(probe.cache)>100_000 && empty!(probe.cache)
    try
        for (channel,info) in TP.SCATTERING_MESON_MAP[process][:channels]
            probe.channels===:s && channel!==:s && continue
            probe.channels===:finite_k && channel===:s && continue
            specs=polarization_specs(info,qp)
            for t in (channel===:s ? (first(ts),) : ts)
                u=ctx.m1^2+ctx.m2^2+ctx.m3^2+ctx.m4^2-s-t
                cms=TP.calculate_cms_momentum(process,s,t,channel,qp;u=u)
                for spec in specs
                    args=(spec.channel,cms.k0,cms.k,spec.m1,spec.m2,spec.μ1,spec.μ2,
                        tp.T,tp.Φ,tp.Φbar,tp.ξ,spec.A1,spec.A2,spec.num_s_quark)
                    key=PC.PolarizationKey(args...)
                    haskey(previous,key) && continue
                    previous[key]=get(PC.POLARIZATION_CACHE,key,nothing)
                    value=get!(probe.cache,key) do
                        result=probe.backend(args...)
                        all(isfinite,result) || throw(DomainError(result,"nonfinite loop probe result"))
                        result
                    end
                    PC.POLARIZATION_CACHE[key]=value
                end
            end
        end
        misses=PC.CACHE_TOTAL_CALLS[]-PC.CACHE_HIT_CALLS[]
        value=CS.total_cross_section(process,s,qp,tp,kc;n_points=n_points)
        if probe.channels===:all
            PC.CACHE_TOTAL_CALLS[]-PC.CACHE_HIT_CALLS[]==misses ||
                error("loop probe did not cover every production polarization query")
        end
        return value
    finally
        for (key,value) in previous
            value===nothing ? delete!(PC.POLARIZATION_CACHE,key) : (PC.POLARIZATION_CACHE[key]=value)
        end
    end
end

end
