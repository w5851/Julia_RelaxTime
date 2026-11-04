from __future__ import annotations

import unittest

from scripts.analysis.relaxtime import merge_transport_eta_refinement as merge


def point(identity="p1", variant="fine", xi="0.36"):
    return dict(point_id=identity, variant=variant, xi=xi, curve_id="a_mu0", T_MeV="200.0",
                muB_MeV="0.0", original_point_id="a0087", tau_u="1.2345678901234567")


class ResultSelectionTests(unittest.TestCase):
    def test_tighter_replaces_exact_point_without_rounding_values(self):
        fine = point()
        tighter = dict(point(variant="tighter"), tau_u="1.2345678901234568")
        additional = point("p2", "tighter", "0.361")
        selected = merge.select_results([fine], [tighter, additional])
        self.assertEqual(selected, [tighter, additional])
        self.assertEqual(fine["tau_u"], "1.2345678901234567")

    def test_changed_coordinate_or_reference_is_rejected(self):
        for field, value in (("xi", "0.361"), ("T_MeV", "201.0"), ("muB_MeV", "450.0"),
                             ("curve_id", "b_mu0"), ("original_point_id", "a0088")):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "does not match"):
                merge.select_results([point()], [dict(point(variant="tighter"), **{field: value})])

    def test_duplicates_at_either_precision_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate result"):
            merge.select_results([point(), point()], [])
        with self.assertRaisesRegex(ValueError, "duplicate physical coordinate"):
            merge.select_results([point(), point("other", xi="0.360")], [])

    def test_precision_label_is_checked(self):
        with self.assertRaisesRegex(ValueError, "precision/source"):
            merge.select_results([point(variant="tighter")], [])

    def test_boundary_audit_counts_coordinates_instead_of_observables(self):
        a = dict(mode="A", muB_MeV=900, curve_parameter=1, xi_left=-0.003, xi_right=0.003,
                 reference_first_order_boundary=True)
        b = dict(a, mode="B", curve_parameter=120, xi_left=-0.13, xi_right=-0.12)
        repeated = [dict(p, observable=field) for p in (a, b) for field in ("tau_u", "eta_over_s", "sigma_over_T")]
        self.assertEqual(len(merge.reviewed_boundaries(dict(known_boundary_pairs=repeated))), 2)
        with self.assertRaisesRegex(ValueError, "unreviewed"):
            merge.reviewed_boundaries(dict(known_boundary_pairs=[dict(a, reference_first_order_boundary=False), b]))


class RowAssemblyTests(unittest.TestCase):
    def template_and_point(self, original_id):
        template = {field: "base metadata" for field in (*merge.META_FIELDS, *merge.PHASE_FIELDS)}
        template.update(T_MeV="200.0", muq_MeV="0.0", muB_MeV="0.0", xi="0.35", iterations="7",
                        phase_prev="quark", phase_curr="hadron", phase_structure="crossover_possible",
                        phase_boundary_xi_used="0.36", mode="mode_a_fixed_muB_phase_scaled")
        selected = dict(point(), original_point_id=original_id, source_run="fine_a_mu0_v2", source_csv_line="4",
                        background_residual="1e-14", s_fm3inv="10.0", T_fm="1.0", muq_fm="0.0",
                        Phi="0.5", Phibar="0.5", m_u="1.0", m_d="1.0", m_s="2.0")
        for name in (*merge.OBSERVABLES, *("n_" + q for q in merge.SPECIES), *("tauinv_" + q for q in merge.SPECIES)):
            selected.setdefault(name, "1.2345678901234567")
        template.update({name: "old" for name in selected if name not in template})
        # The real base scan has no selection-only columns.
        for name in ("variant", "point_id", "curve_id", "original_point_id", "source_run", "source_csv_line", "background_residual"):
            template.pop(name, None)
        return template, selected

    def test_replay_keeps_original_phase_metadata_but_uses_saved_transport(self):
        template, selected = self.template_and_point("a0087")
        assembled = merge.assembled_row(selected, template, template, {"eta_over_s": "bad override"}, template["mode"])
        for name in merge.OBSERVABLES:
            self.assertEqual(assembled[name], selected[name])
        for name in merge.PHASE_FIELDS:
            self.assertEqual(assembled[name], template[name])
        self.assertEqual(assembled["point_id"], "a0087")
        self.assertEqual(assembled["source_kind"], "recomputed_original")
        self.assertEqual(assembled["iterations"], "")
        self.assertEqual(template["iterations"], "7")

    def test_insert_does_not_invent_phase_classification(self):
        template, selected = self.template_and_point("")
        assembled = merge.assembled_row(selected, None, template, {}, template["mode"])
        self.assertEqual(assembled["point_id"], "eta_p1")
        self.assertEqual(assembled["source_kind"], "xi_inserted")
        self.assertEqual(assembled["phase_curr"], "unknown")
        self.assertEqual(assembled["phase_boundary_xi_used"], "")


if __name__ == "__main__":
    unittest.main()
