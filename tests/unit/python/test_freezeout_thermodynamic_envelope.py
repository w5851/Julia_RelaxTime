"""Scientific identities for projecting frozen freezeout parameters into T-muB."""
import math
import unittest

from scripts.analysis.relaxtime.freezeout_parameter_envelope import trajectory
from scripts.analysis.relaxtime.plot_freezeout_thermodynamic_envelope import (
    andronic_temperature, fixed_temperature_band, polynomial_temperature,
)


class ThermodynamicEnvelopeTests(unittest.TestCase):
    def setUp(self):
        self.andronic = dict(T_lim_GeV=.1584,log_offset=2.6,log_width=.45,d_GeV=1.3075,e_GeV_inv1=.288)

    def test_polynomial_projection_removes_only_energy_mapping(self):
        first = [.15,.10,.05,.75,.27]
        second = [.15,.10,.05,1.1,.20]
        for mu in (0.,.2,.5):
            self.assertEqual(polynomial_temperature(first,mu),polynomial_temperature(second,mu))
        self.assertNotEqual(trajectory(first,7.7),trajectory(second,7.7))
        t,mu = trajectory(first,7.7)
        self.assertEqual(t,polynomial_temperature(first,mu))

    def test_andronic_inverse_matches_energy_form_and_zero_limit(self):
        for energy in (4.87,7.7,19.6,200.,5020.):
            expected = .1584/(1+math.exp(2.6-math.log(energy)/.45))
            mu = 1.3075/(1+.288*energy)
            self.assertAlmostEqual(andronic_temperature(self.andronic,mu),expected,places=14)
        self.assertEqual(andronic_temperature(self.andronic,0.),.1584)
        self.assertAlmostEqual(andronic_temperature(self.andronic,1e-300),.1584,places=14)

    def test_fixed_band_preserves_members_and_exact_units(self):
        candidates = {'one':dict(values=[.16,.10,.05,1.,.2]),'two':dict(values=[.15,.12,.07,.8,.3])}
        rows = fixed_temperature_band(candidates,['one','two'],[0.,.2])
        self.assertEqual(rows[0]['minimum_T_MeV'],150.)
        self.assertEqual(rows[0]['maximum_T_MeV'],160.)
        self.assertAlmostEqual(rows[1]['minimum_T_MeV'],145.088,places=11)
        self.assertAlmostEqual(rows[1]['maximum_T_MeV'],155.92,places=11)
        self.assertEqual([r['member_count'] for r in rows],[2,2])
        self.assertEqual(rows[1]['minimum_profile'],'two')

    def test_bad_domains_do_not_silently_drop_a_member(self):
        candidates = {'one':dict(values=[.15,.2,1.,1.,.2])}
        for grid in ([.1,.1],[-.1,.1],[0.,.8]):
            with self.assertRaises(ValueError):
                fixed_temperature_band(candidates,['one'],grid)
        with self.assertRaises(ValueError):
            fixed_temperature_band(candidates,['one','one'],[0.,.2])
        for mu in (-.1,1.3075,float('nan')):
            with self.assertRaises(ValueError):
                andronic_temperature(self.andronic,mu)


if __name__ == '__main__':
    unittest.main()
