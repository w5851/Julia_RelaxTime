# Analytic k=0, gamma=0 endpoint of the historical energy-symmetrized K loop.
# At w=m1+m2, D_i(+w)^(-1)+D_i(-w)^(-1)=-m_i/(w*p^2).
# Hence B0_sym=2/w * sum_i m_i*integral_0^Lambda[(1-f_i-fbar_i)/E_i] dp.
# Occupations include the same linear-xi correction as PairLoopQuadrature.
# This evaluator is for diagnostic threshold location, not transport production.
function kaon_endpoint_inverse(state; n_nodes=128)
    qp,tp,kc=state.qp,state.tp,state.kc
    abs(qp.μ.u-qp.μ.s)<1e-14 || error("kaon endpoint requires common chemical potential")
    mass=qp.m.u+qp.m.s
    plans=map((:u,:s)) do f
        PAIR_OTHER.PairLoopPlan(getproperty(qp.m,f),getproperty(qp.μ,f),tp;n_nodes=n_nodes)
    end
    b0=2/mass*(qp.m.u*plans[1].threshold_integral+qp.m.s*plans[2].threshold_integral)
    pol=-Main.Constants_PNJL.N_color/(8pi^2)*(qp.A.u+qp.A.s-4qp.m.u*qp.m.s*b0)
    return complex(1-4MP._meson_coupling(:K,kc)*pol)
end

function kaon_pair_inverse(state, mass; n_nodes=128)
    qp,tp,kc=state.qp,state.tp,state.kc
    s=mass^2
    s>abs(qp.m.u^2-qp.m.s^2) || error("diagnostic formula restricted to near pair threshold")
    b0=zero(ComplexF64)
    for (f,other) in ((:u,:s),(:s,:u))
        m,mother=getproperty(qp.m,f),getproperty(qp.m,other)
        a=s+m^2-mother^2
        plan=PAIR_OTHER.PairLoopPlan(m,getproperty(qp.μ,f),tp;n_nodes=n_nodes)
        b0+=a/(2s)*PAIR_OTHER.pair_b0_pv(plan,a^2/s)
    end
    pol=-Main.Constants_PNJL.N_color/(8pi^2)*(qp.A.u+qp.A.s+((qp.m.u-qp.m.s)^2-s)*b0)
    return 1-4MP._meson_coupling(:K,kc)*pol
end
