# Independent reference integrations; QuadGK remains confined to benchmark/.
const LOOP_REFERENCE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
joinpath(LOOP_REFERENCE_ROOT, "benchmark") in LOAD_PATH ||
    push!(LOAD_PATH, joinpath(LOOP_REFERENCE_ROOT, "benchmark"))
if !isdefined(Main, :ThresholdUpgradeCommon)
    include(joinpath(LOOP_REFERENCE_ROOT,"scripts","analysis","relaxtime","ThresholdUpgradeCommon.jl"))
end

module LoopQuadratureReference

using QuadGK
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: correction_cos_theta_coefficient
using Main.Constants_PNJL: Λ_inv_fm, N_color

export loop_reference, polarization_reference, angular_pv_second_moment, angular_log

@inline function occupation_components(sign, p, m, mu, tp)
    E = sqrt(m*m+p*p)
    if sign===:plus
        d = quark_distribution(E,mu,tp.T,tp.Φ,tp.Φbar)
        c = tp.ξ==0 ? 0.0 : correction_cos_theta_coefficient(:quark,p,m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ)
    else
        d = 1-antiquark_distribution(E,mu,tp.T,tp.Φ,tp.Φbar)
        c = tp.ξ==0 ? 0.0 : -correction_cos_theta_coefficient(:antiquark,p,m,mu,tp.T,tp.Φ,tp.Φbar,tp.ξ)
    end
    return d,c
end

@inline function angular_log(A, B)
    A==0 && return 0.0
    B!=0 && abs(A/B)<0.5 && return 2atanh(A/B)
    return log(abs(B+A))-log(abs(B-A))
end

"""PV integral of c^2/(A*c+B), evaluated without small-A cancellation."""
function angular_pv_second_moment(A, B)
    if B!=0 && abs(A/B)<0.125
        r2=(A/B)^2
        term=1.0/3; total=term; power=1.0
        for n in 1:24
            power*=r2
            term=power/(2n+3)
            next=total+term
            next==total && break
            total=next
        end
        return 2total/B
    end
    B==0 && return 0.0
    return (-2B+B*B/A*angular_log(A,B))/(A*A)
end

function k0_component(sign, lambda, m, mp, mu, tp; rtol=2e-10)
    L=Λ_inv_fm
    dist(p)=let (d,c)=occupation_components(sign,p,m,mu,tp); d+c/3 end
    mass_difference=(m-mp)*(m+mp)
    c=(lambda^2+mass_difference)/2
    if lambda==0
        c!=0 || throw(DomainError((lambda,m,mp),"static degenerate k=0 reference requires a specified limiting prescription"))
        value,error=quadgk(p->2p*p*dist(p)/(sqrt(m*m+p*p)*c),0.0,L;rtol=rtol,atol=1e-12)
        return complex(value,0.0),error
    end
    e0=-0.5*(lambda+mass_difference/lambda)
    if e0<=0
        value,error=quadgk(p->2p*p*dist(p)/(sqrt(m*m+p*p)*(lambda*sqrt(m*m+p*p)+c)),
            0.0,L;rtol=rtol,atol=1e-12)
        return complex(value,0.0),error
    end
    # Factor the Kallen polynomial: evaluating E0^2-m^2 at the
    # threshold otherwise manufactures an O(eps) offset and an
    # O(sqrt(eps)) discontinuity in the reference itself.
    q2=(lambda-(m+mp))*(lambda+(m+mp))*(lambda-(m-mp))*(lambda+(m-mp))/(4lambda^2)
    G(p)=dist(p)*(sqrt(m*m+p*p)+e0)/sqrt(m*m+p*p)
    if q2<0
        kappa=sqrt(-q2)
        edges=sort!(unique(vcat([0.0,L],[p for p in (kappa,10kappa) if 0<p<L])))
        value,error=quadgk(p->2/lambda*p*p*G(p)/(p*p+kappa*kappa),edges...;rtol=rtol,atol=1e-12)
        return complex(value,0.0),error
    elseif q2==0
        value,error=quadgk(p->2/lambda*G(p),0.0,L;rtol=rtol,atol=1e-12)
        return complex(value,0.0),error
    end
    q=sqrt(q2)
    q==L && throw(DomainError(q,"sharp upper-cutoff logarithmic divergence"))
    h(p)=p*p*G(p)/(p+q)
    hq=q*G(q)/2
    regular(p)=if p==q
        step=cbrt(eps(Float64))*max(q,1.0)
        (h(q+step)-h(max(0.0,q-step)))/(q+step-max(0.0,q-step))
    else
        (h(p)-hq)/(p-q)
    end
    edges=q<L ? (0.0,q,L) : (0.0,L)
    value,error=quadgk(regular,edges...;rtol=rtol,atol=1e-12)
    value=2/lambda*(value+hq*(log(abs(L-q))-log(q)))
    imaginary=q<L ? 2pi*q*dist(q)/lambda : 0.0
    return complex(value,imaginary),2error/abs(lambda)
end

function radial_breakpoints(lambda,k,m,mp)
    # Solve (2*lambda*E+C)^2 = 4*k^2*(E^2-m^2), independently
    # of the production root/interval helper.
    a=lambda^2-k^2; C=a+m^2-mp^2
    b=lambda*C; c=C*C/4+k*k*m*m
    roots=Float64[]
    if a==0
        b!=0 && push!(roots,-c/b)
    else
        discriminant=b*b-4a*c
        if discriminant>=0
            q=-0.5*(b+copysign(sqrt(discriminant),b))
            q==0 ? push!(roots,-b/(2a)) : append!(roots,(q/a,c/q))
        end
    end
    Emax=sqrt(m*m+Λ_inv_fm^2)
    p=[sqrt((E-m)*(E+m)) for E in roots if m<E<Emax]
    return sort!(unique(vcat(0.0,p,Λ_inv_fm)))
end

function finite_k_component(sign,lambda,k,m,mp,mu,tp;rtol=2e-10)
    edges=radial_breakpoints(lambda,k,m,mp)
    function real_integrand(p)
        E=sqrt(m*m+p*p); A=2p*k
        B=lambda^2+2lambda*E+m*m-mp*mp-k*k
        d,c=occupation_components(sign,p,m,mu,tp)
        correction=c==0 ? 0.0 : 2p*c*angular_pv_second_moment(A,B)
        return p/E*(d*angular_log(A,B)/k+correction)
    end
    real_value,real_error=quadgk(real_integrand,edges...;rtol=rtol,atol=1e-11)
    imaginary=0.0;imaginary_error=0.0
    if lambda!=0
        for (a,b) in zip(edges[1:end-1],edges[2:end])
            p=(a+b)/2;E=sqrt(m*m+p*p);A=2p*k
            B=lambda^2+2lambda*E+m*m-mp*mp-k*k
            abs(B)<A || continue
            value,error=quadgk(a,b;rtol=rtol,atol=1e-11) do p
                E=sqrt(m*m+p*p);A=2p*k
                B=lambda^2+2lambda*E+m*m-mp*mp-k*k
                d,c=occupation_components(sign,p,m,mu,tp)
                pi*Base.sign(lambda)*p/E*(d/k+2p*c*B*B/(A*A*A))
            end
            imaginary+=value;imaginary_error+=error
        end
    end
    return complex(real_value,imaginary),real_error+imaginary_error
end

function loop_reference(lambda,k,m1,m2,mu1,mu2,tp;rtol=2e-10)
    k>=0 || throw(ArgumentError("nonnegative momentum required"))
    component(sign,lambda,m,mp,mu)= k==0 ?
        k0_component(sign,lambda,m,mp,mu,tp;rtol=rtol) :
        finite_k_component(sign,lambda,k,m,mp,mu,tp;rtol=rtol)
    terms=(component(:plus,-lambda,m1,m2,mu1),component(:minus,lambda,m1,m2,mu1),
        component(:plus,lambda,m2,m1,mu2),component(:minus,-lambda,m2,m1,mu2))
    return (value=terms[1][1]-terms[2][1]+terms[3][1]-terms[4][1],
        error=sum(t[2] for t in terms))
end

function polarization_reference(channel,k0,k,m1,m2,mu1,mu2,tp,a1,a2,ns;rtol=2e-10)
    lambda=k0+mu1-mu2
    reference=loop_reference(lambda,k,m1,m2,mu1,mu2,tp;rtol=rtol)
    b0=reference.value;error=reference.error
    if ns==1
        extra=loop_reference(-k0+mu1-mu2,k,m1,m2,mu1,mu2,tp;rtol=rtol)
        b0=(b0+extra.value)/2;error=(error+extra.error)/2
    end
    coefficient=k*k-lambda*lambda+(channel===:P ? (m1-m2)^2 : (m1+m2)^2)
    factor=-N_color/(8pi^2)
    return (value=factor*(a1+a2+coefficient*b0),error=abs(factor*coefficient)*error)
end

end
