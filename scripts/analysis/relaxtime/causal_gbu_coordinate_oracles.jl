"""Independent radial-coordinate loop oracle on an explicitly unchanged domain.

No density, equilibrium solve, phase choice or production default lives here.
For q>0 change variables from (p1,cos(theta)) to (p1,p2), with measure
p1*p2*dp1*dp2/(q*E1*E2). This tests a coordinate routing, not a new regulator.
"""
module CausalGBUCoordinateOracles
include("causal_gbu_regulator_checks.jl")
const C=CausalGBURegulatorChecks
const R=C.R
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution

function radial_domain(q,L,regulator)
    isfinite(q) && isfinite(L) && q>=0 && L>0 || throw(ArgumentError("invalid domain"))
    regulator in (:two_line,:centered) || throw(ArgumentError("unknown regulator"))
    if regulator===:two_line
        q>=2L && return Float64[]
        return sort!(unique([max(0.,q-L),abs(L-q),L]))
    end
    # The p2 upper bound changes between the triangle and the center circle here.
    lo,hi=max(0.,q/2-L),q/2+L
    edges=[lo,hi]
    lo<L-q/2<hi && push!(edges,L-q/2)
    lo<q<hi && push!(edges,q)
    return sort!(unique!(edges))
end

function radial_limits(p,q,L,regulator)
    low=abs(p-q)
    high=regulator===:two_line ? min(L,p+q) :
        min(p+q,sqrt(max(0.,2L^2+q^2/2-p^2)))
    return low,high
end

"""Return vacuum and thermal Pi/B0/contact separately at INTERNAL complex lambda.

The radial cut boundaries are analytic; no indicator function is sampled by
quadrature. `component=:vacuum` means remove occupations at fixed masses, not
a separately solved physical vacuum background.
"""
function radial_loop(q,a,mu1,b,mu2,T,phi,phibar,vc,tc,zs,regulator;np=128,nr=128)
    all(isfinite,(q,a,mu1,b,mu2,T,phi,phibar,vc,tc)) && q>=0 && min(a,b,T,vc)>0 &&
        tc>=vc && min(np,nr)>=4 && 0<=phi<=1 && 0<=phibar<=1 || throw(ArgumentError("invalid loop inputs"))
    !isempty(zs) && all(z->isfinite(z) && imag(z)>0,zs) || throw(ArgumentError("UHP lambda probes required"))
    radial_domain(q,vc,regulator) # validate even when q=0
    results=NamedTuple[]
    for (component,L) in ((:vacuum,vc),(:thermal,tc))
        values=zeros(ComplexF64,length(zs))
        bvalues=zeros(ComplexF64,length(zs))
        moment0,moment1,contact=0.,0.,0.
        function add(p,r,weight)
            e1,e2=hypot(p,a),hypot(r,b)
            n1=component===:vacuum ? (1.,0.) :
                (-antiquark_distribution(e1,mu1,T,phi,phibar),quark_distribution(e1,mu1,T,phi,phibar))
            n2=component===:vacuum ? (1.,0.) :
                (-antiquark_distribution(e2,mu2,T,phi,phibar),quark_distribution(e2,mu2,T,phi,phibar))
            measure=weight/(e1*e2)
            # Tadpole reduction is independent of the following residue moments.
            contact+=2measure*(e1*(n2[2]-n2[1])+e2*(n1[2]-n1[1]))
            for (si,s) in enumerate((-1,1)),(ti,t) in enumerate((-1,1))
                u=s*e1-t*e2
                residue=measure*s*t*(n2[ti]-n1[si])
                moment0+=residue
                moment1+=residue*u
                trace=u^2-q^2-(a-b)^2
                for i in eachindex(zs)
                    term=residue/(zs[i]-u)
                    bvalues[i]+=term
                    values[i]+=C.NC/(8pi^2)*trace*term
                end
            end
        end
        if q==0
            ps,ws=R.gauleg(0.,L,np)
            for (p,w) in zip(ps,ws)
                add(p,p,2w*p^2)
            end
        else
            edges=radial_domain(q,L,regulator)
            for i in 1:length(edges)-1
                ps,ws=R.gauleg(edges[i],edges[i+1],np)
                for (p,w) in zip(ps,ws)
                    left,right=radial_limits(p,q,L,regulator)
                    left<right || continue
                    rs,rws=R.gauleg(left,right,nr)
                    for (r,rw) in zip(rs,rws)
                        add(p,r,w*rw*p*r/q)
                    end
                end
            end
        end
        reconstructed=[C.NC/(8pi^2)*((z^2-q^2-(a-b)^2)*bvalues[i]-contact) for (i,z) in enumerate(zs)]
        push!(results,(component=component,values=values,b0=bvalues,contact=contact,
            moment0=moment0,moment1=moment1,
            contact_identity_residual=maximum(abs.(values.-reconstructed)),
            moment_identity_residual=abs(moment1-contact)))
    end
    return (vacuum=results[1],thermal=results[2],
        values=results[1].values.+results[2].values,
        contact=results[1].contact+results[2].contact)
end
end
