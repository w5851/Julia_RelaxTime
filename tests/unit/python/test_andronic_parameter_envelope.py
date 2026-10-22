"""Family-specific formulas, deterministic coverage and exact energy remapping."""
import math
import unittest

from scripts.analysis.relaxtime.andronic_parameter_envelope import (
    DEFAULT_CONFIG, read_config, initial_candidates, rescale_energy,
)
from scripts.analysis.relaxtime.freezeout_parameter_envelope import (
    ANDRONIC_FAMILY, trajectory, trajectory_valid, in_box,
)
from scripts.analysis.relaxtime.plot_freezeout_thermodynamic_envelope import (
    andronic_temperature, fixed_temperature_band, geometry_signature,
)


class AndronicExplorationTests(unittest.TestCase):
    def test_energy_formula_and_inverse_agree_without_moving_the_division(self):
        p=[.1566,2.20,.48,1.3075,.288]
        keys=('T_lim_GeV','log_offset','log_width','d_GeV','e_GeV_inv1')
        for energy in (4.87,7.7,19.6,200.,5020.):
            t,mu=trajectory(p,energy,ANDRONIC_FAMILY)
            self.assertAlmostEqual(t,p[0]/(1+math.exp(p[1]-math.log(energy)/p[2])),places=14)
            self.assertAlmostEqual(t,andronic_temperature(dict(zip(keys,p)),mu),places=14)
        self.assertTrue(trajectory_valid(p,[4.87,200.],ANDRONIC_FAMILY))
        self.assertFalse(trajectory_valid([.15,2.6,0.,1.,.3],[4.87,200.],ANDRONIC_FAMILY))
        with self.assertRaises(ValueError):trajectory(p,7.7,'unknown_family')

    def test_design_covers_both_faces_and_uses_only_complete_nonbaseline_anchors(self):
        config,domain=read_config(DEFAULT_CONFIG)
        candidates,bounds,inner=initial_candidates(config,domain)
        self.assertEqual((candidates,bounds,inner),initial_candidates(config,domain))
        self.assertEqual(len(candidates),46)
        self.assertEqual(sum(c['role']=='source_anchor' for c in candidates),4)
        self.assertEqual(sum(c['role']=='interior_sample' for c in candidates),32)
        self.assertNotIn('thakur_hirono2022',{c['name'] for c in candidates})
        self.assertNotIn('andronic2018',{c['name'] for c in candidates})
        self.assertEqual(len({tuple(c['values']) for c in candidates}),46)
        for i,pair in enumerate(bounds):
            self.assertEqual(min(c['values'][i] for c in candidates),pair[0])
            self.assertEqual(max(c['values'][i] for c in candidates),pair[1])
        for c in candidates:
            self.assertTrue(in_box(c['values'],bounds))
            self.assertTrue(trajectory_valid(c['values'],[4.87,200.],ANDRONIC_FAMILY))

    def test_joint_A_E_remapping_preserves_the_full_thermodynamic_curve(self):
        p=[.156,2.6,.45,1.2,.29]
        mapped=rescale_energy(p,8.,7.7)
        self.assertEqual([mapped[i] for i in (0,2,3)],[p[i] for i in (0,2,3)])
        k=8./7.7
        for energy in (4.87,7.7,19.6,200.):
            before=trajectory(p,k*energy,ANDRONIC_FAMILY)
            after=trajectory(mapped,energy,ANDRONIC_FAMILY)
            for a,b in zip(before,after):self.assertAlmostEqual(a,b,places=14)
        self.assertEqual(geometry_signature(p,ANDRONIC_FAMILY),geometry_signature(mapped,ANDRONIC_FAMILY))
        with self.assertRaises(ValueError):rescale_energy(p,8.,0.)

    def test_fixed_mu_band_uses_andronic_geometry_and_keeps_every_member(self):
        p=[.1566,2.2,.48,1.3075,.288]
        q=[.1535,2.82,.4,1.3568,.326]
        candidates={'one':dict(values=p),'two':dict(values=q)}
        rows=fixed_temperature_band(candidates,['one','two'],[0.,.2,.5],ANDRONIC_FAMILY)
        self.assertEqual(rows[0]['minimum_T_MeV'],153.5)
        self.assertAlmostEqual(rows[0]['maximum_T_MeV'],156.6,places=12)
        for row,mu in zip(rows,[0.,.2,.5]):
            values=[1000*andronic_temperature(dict(zip(('T_lim_GeV','log_offset','log_width','d_GeV','e_GeV_inv1'),v)),mu) for v in (p,q)]
            self.assertAlmostEqual(row['minimum_T_MeV'],min(values),places=12)
            self.assertAlmostEqual(row['maximum_T_MeV'],max(values),places=12)
            self.assertEqual(row['member_count'],2)
        with self.assertRaises(ValueError):
            fixed_temperature_band(candidates,['one','two'],[0.,1.32],ANDRONIC_FAMILY)


if __name__=='__main__':unittest.main()
