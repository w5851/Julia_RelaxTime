#!/usr/bin/env julia
# quote-remote only. Independent continuum thermal weight in total-pair variables.
include(joinpath(@__DIR__, isfile(joinpath(@__DIR__,"run_sigma_mott_spectrum_v2.jl")) ?
    "run_sigma_mott_spectrum_v2.jl" : "run_sigma_mott_spectrum.jl"))

"""Return W(s)/x for equal-mass incoming u/ubar, x=s-4m^2.

The original measure is d^3p1 d^3p2 f1 f2 v_cm. With Q=p1+p2,
E=sqrt(s+Q^2), k*=sqrt(x)/2 and beta^2=x/s, its Jacobian is
d^3p1 d^3p2 = E1 E2 beta/(4E) ds d^3Q dOmega*, while v_cm=2beta.
The common azimuth gives the same Dq^2/(32pi^5 n1 n2) prefactor as AS.
Thus W(s)/x = prefactor/(2s) integral Q^2 dQ dcQ dc* dphi* E1 E2/E f1 f2.
This changes integration variables only; it does not change the collision kernel.
"""
function sigma_continuum_g(state, x, prefactor; p_nodes=64, c_nodes=16, internal_nodes=12,
                           boltzmann=false, swap=false)
    m=state.qp.m.u
    s=4m^2+x
    root_s=sqrt(s)
    kstar=sqrt(x)/2
    q,qw=AS._cached_semi_infinite_momentum_grid(p_nodes,20.0)
    cq,cqw=AS._cached_interval_grid(-1.0,1.0,c_nodes)
    c,cw=AS._cached_interval_grid(-1.0,1.0,internal_nodes)
    phi,phiw=AS._cached_interval_grid(0.0,2pi,internal_nodes)
    fi,fj=swap ? (:ubar,:u) : (:u,:ubar)
    mui,muj=AS.get_mu(fi,state.qp),AS.get_mu(fj,state.qp)
    tp=state.tp
    total=0.0
    for i in eachindex(q)
        Q=q[i]
        E=sqrt(s+Q^2)
        gamma=E/root_s
        radial=qw[i]*Q^2/E
        for a in eachindex(cq), b in eachindex(c), d in eachindex(phi)
            E1=E/2+Q*kstar/root_s*c[b]
            E2=E/2-Q*kstar/root_s*c[b]
            fproduct=if boltzmann
                exp(-E/tp.T)
            else
                transverse=kstar*sqrt(1-c[b]^2)*cos(phi[d])*sqrt(1-cq[a]^2)
                pz1=(Q/2+gamma*kstar*c[b])*cq[a]-transverse
                pz2=(Q/2-gamma*kstar*c[b])*cq[a]+transverse
                p1=sqrt(max(E1^2-m^2,0.0));p2=sqrt(max(E2^2-m^2,0.0))
                ct1=p1==0 ? 0.0 : clamp(pz1/p1,-1.0,1.0)
                ct2=p2==0 ? 0.0 : clamp(pz2/p2,-1.0,1.0)
                AS.distribution_with_anisotropy(fi,p1,m,mui,tp.T,tp.Φ,tp.Φbar,tp.ξ,ct1)*
                AS.distribution_with_anisotropy(fj,p2,m,muj,tp.T,tp.Φ,tp.Φbar,tp.ξ,ct2)
            end
            total+=radial*cqw[a]*cw[b]*phiw[d]*E1*E2*fproduct
        end
    end
    prefactor/(2s)*total
end

# Independent analytical angular reduction for Maxwell-Boltzmann input.
function sigma_boltzmann_g_oracle(state,x,prefactor;p_nodes=64)
    s=4state.qp.m.u^2+x
    q,w=AS._cached_semi_infinite_momentum_grid(p_nodes,20.0)
    integral=0.0
    for (Q,weight) in zip(q,w)
        E=sqrt(s+Q^2)
        integral+=weight*Q^2/E*(E^2-Q^2*x/(3s))*exp(-E/state.tp.T)
    end
    prefactor*pi/s*integral
end

function sigma_threshold_g_oracle(state,prefactor;p_nodes=64,c_nodes=16)
    s0=4state.qp.m.u^2
    q,w=AS._cached_semi_infinite_momentum_grid(p_nodes,20.0)
    c,cw=AS._cached_interval_grid(-1.0,1.0,c_nodes)
    tp=state.tp;m=state.qp.m.u
    total=0.0
    for (Q,weight) in zip(q,w),(ct,angular) in zip(c,cw)
        E=sqrt(s0+Q^2)
        product=AS.distribution_with_anisotropy(:u,Q/2,m,AS.get_mu(:u,state.qp),tp.T,tp.Φ,tp.Φbar,tp.ξ,ct)*
                AS.distribution_with_anisotropy(:ubar,Q/2,m,AS.get_mu(:ubar,state.qp),tp.T,tp.Φ,tp.Φbar,tp.ξ,ct)
        total+=weight*angular*Q^2*E*product
    end
    prefactor*pi/(2s0)*total
end

struct SigmaWeightTable
    nodes::Vector{Float64}
    values::Vector{Float64}
    weights::Vector{Float64}
end
function (table::SigmaWeightTable)(x)
    exact=findfirst(==(x),table.nodes)
    exact===nothing || return table.values[exact]
    numerator=0.0;denominator=0.0
    for i in eachindex(table.nodes)
        weight=table.weights[i]/(x-table.nodes[i])
        numerator+=weight*table.values[i];denominator+=weight
    end
    numerator/denominator
end

function sigma_weight_table(state,prefactor,xmax;fine=false,order=12)
    nodes=[xmax/2*(1-cos(pi*j/order)) for j in 0:order]
    values=[sigma_continuum_g(state,x,prefactor;p_nodes=fine ? 96 : 64,
        c_nodes=fine ? 24 : 16,internal_nodes=fine ? 20 : 12) for x in nodes]
    weights=[(-1.)^j*(j in (0,order) ? 0.5 : 1.0) for j in 0:order]
    SigmaWeightTable(nodes,values,weights)
end

function sigma_gl_vector(f,a,b,n)
    nodes,weights=SIG_CS.standard_nodes_weights(n)
    value=zeros(6)
    for (z,w) in zip(nodes,weights)
        value .+= (b-a)/2*w .* f((a+b)/2+(b-a)/2*z)
    end
    value
end

function sigma_integrate_vector(f,a,b;atol=1e-13,rtol=1e-8,depth=0)
    coarse=sigma_gl_vector(f,a,b,16)
    fine=sigma_gl_vector(f,a,b,32)
    defect=abs.(fine.-coarse)
    if all(defect .<= atol .+ rtol.*abs.(fine))
        return (value=fine,error=defect,leaves=1)
    end
    depth<14 || error("Continuous sigma integral did not converge on [$a,$b]")
    left=sigma_integrate_vector(f,a,(a+b)/2;atol=atol/2,rtol=rtol,depth=depth+1)
    right=sigma_integrate_vector(f,(a+b)/2,b;atol=atol/2,rtol=rtol,depth=depth+1)
    (value=left.value.+right.value,error=left.error.+right.error,leaves=left.leaves+right.leaves)
end

function run_sigma_thermal_audit()
    spectrum=joinpath(dirname(OUTPUT_DIR),"spectrum_v2")
    isfile(joinpath(spectrum,"summary.json")) || error("Completed spectrum is required")
    profile=SP_OTHER.load_profile(joinpath(INPUT_DIR,"s6.toml"))
    runtime=SP_OTHER.build_runtime(profile)
    backgrounds=collect(CSV.File(joinpath(spectrum,"background.csv")))
    allroots=collect(CSV.File(joinpath(spectrum,"real_axis_roots.csv")))
    selected=get(ENV,"MOTT_POINT_IDS","all")
    selected=="all" || filter!(r->String(r.point_id) in split(selected,','),backgrounds)
    bounds=[1e-9,1e-4,5e-4,1e-3,2e-3,5e-3,1e-2]
    xmax=last(bounds)
    results,weightrows,checks,events,peakrows,history=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    function checkpoint()
        for (name,rows) in (("thermal_band_comparison",results),("continuous_weight",weightrows),
                           ("verification",checks),("direct_peak_events",events),("thermal_peak_areas",peakrows),
                           ("discrete_cache_history",history))
            isempty(rows) || CSV.write(joinpath(OUTPUT_DIR,name*".csv"),rows)
        end
    end
    for row in backgrounds
        started=time()
        state=frozen_state(row,profile,32,16)
        info=(point_id=String(row.point_id),xi=Float64(row.xi),delta_xi=Float64(row.delta_xi))
        s0=4state.qp.m.u^2
        problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,state.qp,state.tp,state.kc;runtime.tau_kwargs...)
        measure=SIG_TS.thermal_measure(problem.plan,problem.context.bounds;relative_velocity=:cm)
        prefactor=problem.plan.prefactor
        near_indices=searchsortedfirst(measure.nodes,s0+first(bounds)):searchsortedlast(measure.nodes,s0+xmax)
        near=SIG_TS.ThermalMeasure(measure.nodes[near_indices],measure.weights[near_indices])
        roots=filter(r->String(r.point_id)==info.point_id && String(r.backend)=="native",allroots)
        gbase=sigma_weight_table(state,prefactor,xmax)
        gfine=sigma_weight_table(state,prefactor,xmax;fine=true)
        for x in (0.,1e-4,0.01)
            numerical=sigma_continuum_g(state,x,prefactor;boltzmann=true)
            analytic=sigma_boltzmann_g_oracle(state,x,prefactor)
            relative=numerical/analytic-1
            abs(relative)<2e-11 || error("Boltzmann angular/Jacobian reduction failed")
            push!(checks,merge(info,(check="Boltzmann_angular_reduction",x_inv_fm2=x,relative_error=relative)))
        end
        threshold=sigma_threshold_g_oracle(state,prefactor)
        abs(gbase(0.)/threshold-1)<2e-11 || error("Threshold thermal weight limit failed")
        push!(checks,merge(info,(check="threshold_limit",x_inv_fm2=0.,relative_error=gbase(0.)/threshold-1)))
        swapped=sigma_continuum_g(state,1e-4,prefactor;swap=true)
        normal=sigma_continuum_g(state,1e-4,prefactor)
        abs(swapped/normal-1)<2e-11 || error("Incoming-particle exchange check failed")
        push!(checks,merge(info,(check="incoming_exchange",x_inv_fm2=1e-4,relative_error=swapped/normal-1)))
        for x in (1.37e-5,3.731e-4,0.006371)
            direct=sigma_continuum_g(state,x,prefactor;p_nodes=96,c_nodes=24,internal_nodes=20)
            relative=gfine(x)/direct-1
            abs(relative)<1e-7 || error("Continuous-weight interpolation is unresolved")
            push!(checks,merge(info,(check="weight_interpolation",x_inv_fm2=x,relative_error=relative)))
        end
        for x in sort!(unique(vcat(gfine.nodes,[Float64(r.x_inv_fm2) for r in roots])))
            push!(weightrows,merge(info,(x_inv_fm2=x,g_base=gbase(x),g_fine=gfine(x),
                W_base_fm2=x*gbase(x),W_fine_fm2=x*gfine(x),quadrature_relative=gbase(x)/gfine(x)-1)))
        end
        features=Float64[]
        for r in roots
            append!(features,[s0+r.x_inv_fm2+y*r.half_width_inv_fm2 for y in
                (-64.,-32.,-16.,-8.,-4.,-2.,-1.,-.5,0.,.5,1.,2.,4.,8.,16.,32.,64.)])
        end
        for process in (:uubar_to_uubar,:uubar_to_ddbar)
            partvalues=Dict{Float64,Any}()
            function parts(s)
                get!(partvalues,s) do
                    SIG_PC.reset_cache!()
                    sigma_components(process,s,state)
                end
            end
            direct=s->parts(s).sigma_fm2
            built=SIG_TS.build_thermal_weighted_table(direct,near,s0;rtol=1e-8,atol=1e-13,
                weighted_window=xmax,seed_order=64,extra_seeds=features,max_evaluations=320000,
                max_nodes=96000,max_rounds=80)
            table=SIG_TS.validated_table(built)
            for h in built.history
                push!(history,merge(info,(channel=String(process),),h))
            end
            CSV.write(joinpath(OUTPUT_DIR,info.point_id*"__"*String(process)*"__near_table.csv"),
                [(s_inv_fm2=s,x_inv_fm2=s-s0,sigma_fm2=v) for (s,v) in zip(table.nodes,table.sigma)])
            weighted=x->begin
                p=parts(s0+x)
                f,b=x*gfine(x),x*gbase(x)
                [f*p.sigma_fm2,f*p.s_mixed_scalar_square_fm2,f*p.mixed_delta_fm2,
                 b*p.sigma_fm2,b*p.s_mixed_scalar_square_fm2,b*p.mixed_delta_fm2]
            end
            for (a,b) in zip(bounds[1:end-1],bounds[2:end])
                splits=sort!(unique(vcat(a,b,[s-s0 for s in features if a<s-s0<b])))
                value=zeros(6);error=zeros(6);leaves=0
                for (lo,hi) in zip(splits[1:end-1],splits[2:end])
                    integrated=sigma_integrate_vector(weighted,lo,hi;atol=1e-12*(hi-lo)/xmax)
                    value .+= integrated.value;error .+= integrated.error;leaves+=integrated.leaves
                end
                indices=searchsortedfirst(near.nodes,s0+a):searchsortedlast(near.nodes,s0+b)
                discrete=sum(near.weights[i]*table(near.nodes[i]) for i in indices)
                # A positive-part scalar table uses the identical dense nodes; it is diagnostic only.
                scalar_table=SIG_TS.WeightedSigmaTable(table.nodes,
                    [parts(s).s_mixed_scalar_square_fm2 for s in table.nodes],s0,s0+xmax)
                scalar_discrete=sum(near.weights[i]*scalar_table(near.nodes[i]) for i in indices)
                push!(results,merge(info,(channel=String(process),x_left_inv_fm2=a,x_right_inv_fm2=b,
                    discrete_full_rate_fm2=discrete,continuous_full_rate_fm2=value[1],
                    continuous_mixed_square_rate_fm2=value[2],continuous_mixed_delta_rate_fm2=value[3],
                    base_weight_full_rate_fm2=value[4],base_weight_mixed_square_rate_fm2=value[5],
                    base_weight_mixed_delta_rate_fm2=value[6],discrete_mixed_square_rate_fm2=scalar_discrete,
                    quadrature_error_full_fm2=error[1],quadrature_error_mixed_fm2=error[2],
                    discrete_total_estimated_error_fm2=built.estimated_error,leaves=leaves,events=length(indices))))
            end
            for r in roots
                a=max(first(bounds),r.x_inv_fm2-32r.half_width_inv_fm2)
                b=min(xmax,r.x_inv_fm2+32r.half_width_inv_fm2)
                a<b || continue
                splits=sort!(unique(vcat(a,b,[s-s0 for s in features if a<s-s0<b])))
                value=zeros(6);err=zeros(6)
                for (lo,hi) in zip(splits[1:end-1],splits[2:end])
                    got=sigma_integrate_vector(weighted,lo,hi;atol=1e-13)
                    value .+= got.value;err .+= got.error
                end
                ids=searchsortedfirst(measure.nodes,s0+a):searchsortedlast(measure.nodes,s0+b)
                length(ids)<=20000 || error("Peak event audit exceeds explicit budget")
                exact=0.0;approximate=0.0;scalar_exact=0.0;l1=0.0
                for i in ids
                    s,w=measure.nodes[i],measure.weights[i]
                    p=parts(s);v=table(s)
                    exact+=w*p.sigma_fm2;approximate+=w*v;scalar_exact+=w*p.s_mixed_scalar_square_fm2
                    l1+=w*abs(p.sigma_fm2-v)
                    push!(events,merge(info,(channel=String(process),s_inv_fm2=s,x_inv_fm2=s-s0,
                        thermal_weight=w,direct_sigma_fm2=p.sigma_fm2,interpolated_sigma_fm2=v,
                        direct_mixed_square_fm2=p.s_mixed_scalar_square_fm2)))
                end
                x=Float64(r.x_inv_fm2)
                native=parts(s0+x).sigma_fm2
                SIG_PC.reset_cache!()
                higher_t=SIG_CS.total_cross_section(process,s0+x,state.qp,state.tp,state.kc;n_points=48)
                push!(checks,merge(info,(check="t24_t48_"*String(process),x_inv_fm2=x,
                    relative_error=higher_t/native-1)))
                push!(peakrows,merge(info,(channel=String(process),x_peak_inv_fm2=x,
                    x_left_inv_fm2=a,x_right_inv_fm2=b,events=length(ids),
                    continuous_full_rate_fm2=value[1],continuous_mixed_square_rate_fm2=value[2],
                    continuous_mixed_delta_rate_fm2=value[3],base_weight_full_rate_fm2=value[4],
                    base_weight_mixed_square_rate_fm2=value[5],quadrature_error_fm2=err[1],
                    discrete_exact_full_rate_fm2=exact,discrete_exact_mixed_square_rate_fm2=scalar_exact,
                    discrete_interpolated_full_rate_fm2=approximate,discrete_interpolation_l1_fm2=l1)))
            end
            checkpoint()
            println("THERMAL_AUDIT_CHANNEL ",info.point_id," ",process," direct_evaluations=",length(partvalues));flush(stdout)
        end
        println("THERMAL_AUDIT_DONE ",info.point_id," seconds=",round(time()-started;digits=2));flush(stdout)
        SIG_PC.reset_cache!();GC.gc()
    end
    write_json(joinpath(OUTPUT_DIR,"summary.json"),(status="completed_diagnostic",points=length(backgrounds),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",production_kernel_changed=false,
        interpretation="Continuum equal-mass thermal weight versus the retained p48/c36/phi36 discrete measure. Kernel, CM relative velocity, blocking and density normalization are identical. This is a diagnostic; no production quadrature is replaced.",
        continuous_weight_quadratures=((64,16,12,12),(96,24,20,20)),x_range_inv_fm2=[first(bounds),last(bounds)],
        omitted_endpoint_interval_inv_fm2=[0.,first(bounds)],script_sha256=digest(@__FILE__),
        spectrum_sha256=digest(joinpath(@__DIR__,"run_sigma_mott_spectrum_v2.jl")),
        library_sha256=SIG_LIBRARY_SHA,
        files=Dict(f=>digest(joinpath(OUTPUT_DIR,f)) for f in readdir(OUTPUT_DIR))))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_sigma_thermal_audit()
end
