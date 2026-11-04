# Diagnostic expressions at k=0, common quark chemical potential, gamma=0.
# Only threshold screening/locations use these; the production kernel is intact.
const PAIR_OTHER = Main.RelaxTime.PairLoopQuadrature

function equal_pair_polarization(state, flavor, channel, mass; n_nodes=128)
    qp, tp = state.qp, state.tp
    m, mu, av = getproperty(qp.m,flavor), getproperty(qp.μ,flavor), getproperty(qp.A,flavor)
    coefficient = channel == :P ? -mass^2 : (2m-mass)*(2m+mass)
    # For a scalar exactly at its own pair threshold the B0 coefficient vanishes.
    loop_term = coefficient == 0.0 ? complex(0.0) : coefficient * PAIR_OTHER.pair_b0_pv(
        PAIR_OTHER.PairLoopPlan(m,mu,tp;n_nodes=n_nodes),mass^2)
    return -Main.Constants_PNJL.N_color/(8pi^2) * (2av + loop_term)
end

function diagnostic_threshold_inverse(state, meson, cut; n_nodes=128, mass_offset=0.0)
    qp,tp,kc = state.qp,state.tp,state.kc
    mass = meson in MIXED ? (cut == :uu ? 2qp.m.u : 2qp.m.s) : MT.mott_threshold_mass(meson,qp)
    if meson in MIXED
        channel = meson in (:eta,:eta_prime) ? :P : :S
        pus = map(f->equal_pair_polarization(state,f,channel,mass+mass_offset;n_nodes=n_nodes),(:u,:s))
        elems = EC.mixing_matrix_elements(pus...,kc,channel)
        which = meson in (:eta,:sigma) ? :light : :heavy
        return MP._mixed_inverse(channel,which,elems)
    elseif meson in (:pi,:sigma_pi)
        channel = meson == :pi ? :P : :S
        pol = equal_pair_polarization(state,:u,channel,mass+mass_offset;n_nodes=n_nodes)
        return 1-4MP._meson_coupling(meson,kc)*pol
    elseif meson == :sigma_K && mass_offset == 0.0
        pol = -Main.Constants_PNJL.N_color/(8pi^2)*(qp.A.u+qp.A.s)
        return complex(1-4MP._meson_coupling(meson,kc)*pol)
    else
        # K is unequal-mass; a finite-offset historical-kernel screen is explicit.
        offset = mass_offset == 0.0 ? -1e-6 : mass_offset
        return MP.meson_mass_equation(meson,mass+offset,0.0,0.0,qp,tp,kc)
    end
end
