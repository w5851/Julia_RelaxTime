"""Numerical contour and weak-limit checks for the infinite-thermal candidate.

Finite resolution tests are not interval-arithmetic or global UV proofs.
"""
module CausalGBUInfiniteAcceptance
include("causal_gbu_infinite_yield.jl")
const Y=CausalGBUInfiniteYield
const P=Y.P
const I=Y.I
const R=Y.R

function nyquist(f,edges,holes,radius;nodes=64)
    radius>0 && nodes>=8 || throw(ArgumentError("positive radius and nodes>=8 required"))
    all(h->-radius<h[1]<h[2]<radius,holes) || throw(ArgumentError("holes must be inside contour"))
    ps=sort!(unique!(vcat([-radius,radius],filter(x->-radius<x<radius,edges),
        [x for h in holes for x in h])))
    filter!(x->!any(h->h[1]<x<h[2],holes),ps)
    path=ComplexF64[]
    for j in 1:length(ps)-1
        a,b=ps[j],ps[j+1]
        if (a,b) in holes
            center=(a+b)/2; h=(b-a)/2
            append!(path,[center+h*cis(pi*(1-i/nodes)) for i in 0:nodes-1])
        else
            append!(path,complex.([a+(b-a)*sinpi(i/(2nodes))^2 for i in 0:nodes-1]))
        end
    end
    append!(path,[radius*cis(pi*i/(4nodes)) for i in 0:4nodes-1])
    values=f.(path)
    all(isfinite,values) && minimum(abs,values)>1e-12 || throw(ArgumentError("singular Nyquist path"))
    increments=[angle(values[mod1(j+1,length(values))]/values[j]) for j in eachindex(values)]
    winding=sum(increments)/(2pi)
    count=round(Int,winding)
    return (count=count,winding=winding,max_step=maximum(abs,increments),
        minimum=minimum(abs,values),passed=abs(winding-count)<1e-7 && maximum(abs,increments)<pi/2)
end

function contour_audit(p;nodes=32,radius=28.,indent=1e-5)
    roots=vcat(Y.gap_roots(p)...)
    centers=copy(roots)
    if p.kernel.q==0
        bg=p.kernel.bg;i,j=R.charged_rpa_spec(p.kernel.ch).pair
        edge=hypot(bg.m[i],bg.vacuum)+hypot(bg.m[j],bg.vacuum)
        append!(centers,[-edge,edge])
    end
    holes=[(c-indent,c+indent) for c in sort(centers)]
    f(z)=P.inverse(p,ComplexF64(z))
    a=nyquist(f,p.kernel.edges,holes,radius;nodes=nodes)
    b=nyquist(f,p.kernel.edges,holes,radius;nodes=2nodes)
    # Contours around each isolated analytic-gap root also include the indented cap.
    rootchecks=map(roots) do root
        reflect(z)=imag(z)>=0 ? f(z) : conj(f(conj(z)))
        R.contour_count(reflect,root-2indent,root+2indent,-2indent,2indent;nodes=32,max_nodes=128)
    end
    return (count=b.count,passed=a.passed && b.passed && a.count==b.count==0 &&
        all(r->r.passed && r.count==1,rootchecks),max_step=max(a.max_step,b.max_step),
        minimum=min(a.minimum,b.minimum),root_disks=length(rootchecks),
        root_disks_passed=all(r->r.passed && r.count==1,rootchecks),
        radius=radius,indent=indent,global_UV_certified=false)
end

"""Finite-eta derivative density, expressed by parts WITH both endpoint terms.

The Bose IR cutoff is explicit. Take eta->0 at fixed lower before lower->0.
The root resonance is integrated here, not added a second time.
"""
function eta_shell(p,eta;lower=1e-3,upper=24.,nodes=64)
    eta>0 && 0<lower<upper || throw(ArgumentError("eta and Bose lower bound must be positive"))
    k=p.kernel;T=k.bg.T
    roots=[r-k.shift for r in vcat(Y.gap_roots(p)...)]
    edges=sort!(unique!(vcat([lower,upper],filter(x->lower<x<upper,k.edges.-k.shift),
        filter(x->lower<x<upper,[r+s*eta for r in roots for s in (-32.,-8.,-2.,0.,2.,8.,32.)]))))
    weight(w)=Y.O.gbu_weight(Y.cut_phase(p,Float64(w);eta=eta))
    g(w)=inv(expm1(w/T))
    part=real(I.C.mapped_integral(edges;nodes=nodes) do w
        g(w)*(1+g(w))/T*weight(w)/pi
    end)
    boundary=(g(upper)*weight(upper)-g(lower)*weight(lower))/pi
    return k.q^2/(2pi^2)*(part+boundary)
end
end
