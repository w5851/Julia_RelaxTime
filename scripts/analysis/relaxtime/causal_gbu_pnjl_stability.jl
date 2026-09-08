"""Analysis-only PNJL and pole-health checks. No phase/density repair or promotion."""
module CausalGBUPNJLStability
include("causal_gbu_gap_completeness.jl")
include("causal_gbu_pereira_reference.jl")
const G=CausalGBUGapCompleteness
const R=G.R
const P=CausalGBUPereiraReference

"""Re-evaluate the three scalar gap identities; never solve or replace a background."""
function tadpole_gap_state(m,A,bare,G,K;Nc=3)
    length(m)==length(A)==length(bare)==3 && all(isfinite,(m...,A...,bare...,G,K,Nc)) &&
        Nc>0 || throw(ArgumentError("finite three-flavor gap inputs required"))
    phi=ntuple(i->Nc*m[i]*A[i]/(4pi^2),3)
    predicted=ntuple(i->bare[i]-4*G*phi[i]+2*K*phi[mod1(i+1,3)]*phi[mod1(i+2,3)],3)
    return (condensates_inv_fm3=phi,masses_inv_fm=predicted,
        K12_fm2=G-K*phi[3]/2,K45_fm2=G-K*phi[2]/2)
end

function _derivative_moment(t)
    if abs(t)<0.01
        term=t
        total=zero(t)
        for n in 1:12
            total+=n*term/(n+1)
            term*=-t
        end
        return total
    end
    return log1p(t)/t-1/(1+t)
end

"""Derivative of the SAME linear spectral interpolant, in an open gap or UHP."""
function spectral_derivative(p,z::Number)
    isfinite(z) && imag(z)>=0 || throw(ArgumentError("closed UHP required"))
    total=zero(complex(z))
    for i in 1:length(p.energy)-1
        a,b=p.energy[i],p.energy[i+1]
        u,v=p.imaginary[i],p.imaginary[i+1]
        u==v==0 && continue
        imag(z)==0 && a<=real(z)<=b && throw(ArgumentError("derivative needs an open gap"))
        t=(b-a)/(a-z)
        total+=u*(1/(a-z)-1/(b-z))+(v-u)*_derivative_moment(t)/(a-z)
    end
    return total/pi
end

"""For D=2K/F and rho_D=Im(D)/pi, a real simple pole has weight -2K/F'.

The sign condition is k0*weight>0, NOT positivity of d(delta)/domega.
"""
function pole_weight(slope,K,k0;minimum_slope=1e-9)
    all(isfinite,(slope,K,k0)) && K>0 && minimum_slope>0 ||
        throw(ArgumentError("invalid pole inputs"))
    simple=abs(slope)>minimum_slope
    weight=simple ? -2K/slope : NaN
    return (simple=simple,weight_fm=weight,sign_passed=simple && k0*weight>0)
end

"""Step-halving derivative check wholly inside an open real gap.

Two successive finite differences AND the analytic derivative must agree.
Changing the verification step does not change the spectral profile or root.
"""
function gap_derivative_check(f,root,left,right,exact;tolerance=1e-5,max_halvings=12)
    all(isfinite,(root,left,right,exact,tolerance)) && left<root<right &&
        tolerance>0 && max_halvings>=1 || throw(ArgumentError("invalid gap derivative check"))
    step=min(1e-5,min(root-left,right-root)/20)
    difference(h)=real(f(root+h)-f(root-h))/(2h)
    previous=difference(step)
    scale=max(1.,abs(exact))
    initial_error=abs(previous-exact)/scale
    numerical=previous;error=initial_error;change=Inf
    for i in 1:max_halvings
        step/=2
        root+step>root>root-step || break
        numerical=difference(step)
        error=abs(numerical-exact)/scale
        change=abs(numerical-previous)/scale
        if isfinite(error) && max(error,change)<tolerance
            return (passed=true,value=numerical,error=error,step_change=change,
                initial_error=initial_error,step_inv_fm=step,halvings=i)
        end
        previous=numerical
    end
    return (passed=false,value=numerical,error=error,step_change=change,
        initial_error=initial_error,step_inv_fm=step,halvings=max_halvings)
end

"""Check k0*rho_Pi>=0 over every linear cell, without changing any weight.

Exact passivity + F(0)>0 excludes UHP zeros by the subtracted dispersion
identity. Neither Float64 checks nor a tolerance-qualified result are an
interval-arithmetic proof for the original continuous kernel.
"""
function passivity_audit(p,K,shift;tolerance=1e-10)
    all(isfinite,(K,shift,tolerance)) && K>0 && tolerance>=0 ||
        throw(ArgumentError("invalid passivity settings"))
    minimum_product=Inf
    negative_segments=0
    for i in 1:length(p.energy)-1
        a,b=p.energy[i]-shift,p.energy[i+1]-shift
        u,v=p.imaginary[i],p.imaginary[i+1]
        slope=(v-u)/(b-a)
        intercept=u-slope*a
        # Use stored endpoint weights: slope*x+intercept can turn an exact
        # zero endpoint into a tiny negative value through cancellation.
        smallest=min(a*u,b*v)
        if slope!=0
            x=-intercept/(2slope)
            if a<x<b
                weight=((b-x)*u+(x-a)*v)/(b-a)
                smallest=min(smallest,x*weight)
            end
        end
        minimum_product=min(minimum_product,smallest)
        smallest<0 && (negative_segments+=1)
    end
    f0=1-4K*R.cauchy_transform(p,shift)
    exact=minimum_product>=0 && imag(f0)==0 && real(f0)>0
    numerical=minimum_product>=-tolerance && abs(imag(f0))<=tolerance && real(f0)>0
    return (minimum_k0_rho=minimum_product,negative_segments=negative_segments,
        static_inverse_real=real(f0),static_inverse_imag=imag(f0),
        sign_sufficient_condition_float64=exact,numerical_passivity=numerical,
        production_authorized=false)
end

"""Upper rectangle count, anchored to spectral knots and tested at two resolutions.

It excludes the strip 0<Im(lambda)<eta_floor; no all-UHP claim is made.
Same-grid quadrature poles are never counted as physical RPA poles.
"""
function upper_rectangle_count(f,radius,eta_floor;anchors=Float64[],nodes=96,max_nodes=1536)
    isfinite(radius) && isfinite(eta_floor) && 0<eta_floor<radius &&
        8<=nodes<=max_nodes && all(isfinite,anchors) || throw(ArgumentError("invalid contour"))
    edge=collect(range(-radius,radius;length=nodes+1))
    append!(edge,[x for x in anchors if -radius<x<radius])
    sort!(unique!(edge))
    zs=complex.(edge,eta_floor)
    append!(zs,[complex(radius,eta_floor+(radius-eta_floor)*i/nodes) for i in 1:nodes])
    append!(zs,[complex(radius-2radius*i/nodes,radius) for i in 1:nodes])
    append!(zs,[complex(-radius,radius-(radius-eta_floor)*i/nodes) for i in 1:nodes-1])
    function one(points)
        # Subdivide the MERGED contour, including every fixed-anchor interval.
        # Increasing only the original uniform grid can leave a steep anchor
        # interval unchanged at all nominal resolutions.
        fine=ComplexF64[]
        for i in eachindex(points)
            push!(fine,points[i],(points[i]+points[mod1(i+1,length(points))])/2)
        end
        function winding(points)
            values=f.(points)
            all(isfinite,values) || return (winding=NaN,step=Inf,minimum=0.)
            increments=[angle(values[mod1(i+1,length(values))]/values[i]) for i in eachindex(values)]
            return (winding=sum(increments)/(2pi),step=maximum(abs,increments),minimum=minimum(abs,values))
        end
        return winding(points),winding(fine),fine
    end
    n=nodes
    while true
        coarse,fine,refined=one(zs)
        finite=isfinite(coarse.winding) && isfinite(fine.winding)
        count=finite ? round(Int,fine.winding) : -1
        resolved=finite && count>=0 && min(coarse.minimum,fine.minimum)>1e-8 &&
            max(coarse.step,fine.step)<pi/4 &&
            max(abs(coarse.winding-count),abs(fine.winding-count))<1e-7
        if resolved || n>=max_nodes
            return (count=count,resolved=resolved,coarse_winding=coarse.winding,
                fine_winding=fine.winding,max_step=max(coarse.step,fine.step),
                minimum_inverse_abs=min(coarse.minimum,fine.minimum),
                contour_nodes=length(refined),eta_floor_inv_fm=eta_floor,
                unresolved_near_axis_strip=true,full_UHP_certified=false)
        end
        2n>max_nodes && return (count=count,resolved=false,coarse_winding=coarse.winding,
            fine_winding=fine.winding,max_step=max(coarse.step,fine.step),
            minimum_inverse_abs=min(coarse.minimum,fine.minimum),
            contour_nodes=length(refined),eta_floor_inv_fm=eta_floor,
            unresolved_near_axis_strip=true,full_UHP_certified=false)
        zs=refined
        n*=2
    end
end

"""Exactly the user's first-line sphere, not the historical shifted-log B0.

First routing: |p1|<L, |p2| unrestricted. No flavor averaging is imposed.
The same domain is used for BOTH A terms and B0; p2 still depends on q.
Two coordinate systems test integration error separately from routing effects.
"""
function single_sphere_loop(q,a,u,b,v,T,phi,bar,vc,tc,zs;nodes=128,coordinates=:angular)
    P.validate(q,a,u,b,v,T,phi,bar,vc,nodes)
    isfinite(tc) && tc>=vc && coordinates in (:angular,:radial) ||
        throw(ArgumentError("invalid single-sphere settings"))
    !isempty(zs) && all(z->isfinite(z) && imag(z)>0,zs) || throw(ArgumentError("UHP probes required"))
    b0=zeros(ComplexF64,length(zs));A1,A2=0.,0.
    for (component,L) in ((:vacuum,vc),(:thermal,tc))
        # Split at p=q so the radial triangle's abs(p-q) cusp is a panel edge.
        edges=sort!(unique!(vcat([0.,L],0<q<L ? [q] : Float64[])))
        for k in 1:length(edges)-1
            ps,ws=R.gauleg(edges[k],edges[k+1],nodes)
            for (p,pw) in zip(ps,ws)
                xs,xws=coordinates===:radial && q>0 ? R.gauleg(abs(p-q),p+q,nodes) :
                    (q==0 ? ([0.],[2.]) : R.gauleg(-1.,1.,nodes))
                for (x,xw) in zip(xs,xws)
                    r=coordinates===:radial && q>0 ? x : sqrt(p^2+q^2-2p*q*x)
                    measure=coordinates===:radial && q>0 ? pw*xw*p*r/q : pw*xw*p^2
                    e1,e2=hypot(p,a),hypot(r,b)
                    n1,n2,c=P.weights(e1,e2,u,v,T,phi,bar,component)
                    A1-=2measure*(c-n1.quark-n1.anti)/e1
                    A2-=2measure*(c-n2.quark-n2.anti)/e2
                    for i in eachindex(zs)
                        pair,scatter=P.pair_scattering(e1,e2,zs[i],n1,n2,c)
                        b0[i]+=measure*(pair+scatter)/(e1*e2)
                    end
                end
            end
        end
    end
    values=[3/(8pi^2)*((z^2-q^2-(a-b)^2)*b0[i]-A1-A2) for (i,z) in enumerate(zs)]
    return (values=values,b0=b0,A1=A1,A2=A2,contact=A1+A2,
        regulator="first_line_sphere_no_shift_of_boundary",production_authorized=false)
end
end
