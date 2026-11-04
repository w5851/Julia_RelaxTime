"""Reference selection must cover both signs and avoid duplicate xi points."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = (Path(os.environ["MOTT_PLAN_PREPARE_SCRIPT"]) if "MOTT_PLAN_PREPARE_SCRIPT" in os.environ
          else Path(__file__).resolve().parents[3]/"scripts/analysis/relaxtime/prepare_transport_other_mott_convergence.py")
FIXTURES = Path(os.environ.get("MOTT_TEST_ROOT", tempfile.gettempdir())).resolve()


def point(name, x):
    return dict(point_id=name,curve_id="toy",meson="K",cut="pair",window="toy",
                T_MeV=200.,muB_MeV=0.,xi=x,xi_mott=0.,delta_xi=x,
                kind="mott_centered",original_point_id="")


def edge(left, right, slope):
    return dict(curve_id="toy",observable="tau_s",point_left=left["point_id"],
                point_right=right["point_id"],delta_xi_midpoint=(left["xi"]+right["xi"])/2,
                log_slope=slope)


def write_csv(path, rows):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class ReferencePlanTests(unittest.TestCase):
    def run_plan(self, points, slopes):
        FIXTURES.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="mott_plan_",dir=FIXTURES) as temporary:
            case=Path(temporary).resolve()
            self.assertTrue(case.is_relative_to(FIXTURES))
            (case/"locate_v1").mkdir()
            (case/"fine_slopes_v1").mkdir()
            write_csv(case/"locate_v1/transport_plan.csv",points)
            write_csv(case/"fine_slopes_v1/slopes.csv",slopes)
            subprocess.run([sys.executable,str(SCRIPT),str(case),str(case/"reference")],
                           check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            result=json.loads((case/"reference/summary.json").read_text(encoding="utf-8"))
            with (case/"reference/transport_plan.csv").open(encoding="utf-8",newline="") as stream:
                plan=list(csv.DictReader(stream))
            return result,plan

    def test_negative_response_extremes_and_outside_window(self):
        points=[point("far",-.008),point("left",-.004),point("inner",-.002),
                point("n1",-.001),point("n2",-.0002),point("zero",0.),
                point("p2",.0002),point("p1",.001)]
        by_id={r["point_id"]:r for r in points}
        slopes=[edge(by_id["far"],by_id["left"],1000.),
                edge(by_id["left"],by_id["inner"],-5.),
                edge(by_id["p2"],by_id["p1"],-.5)]
        result,plan=self.run_plan(points,slopes)
        job=result["jobs"][0]
        self.assertTrue({"left","inner","n1","n2","zero","p2","p1"} <= set(job["point_ids"]))
        self.assertNotIn("far",job["point_ids"])
        self.assertEqual(job["new_xi_points"],2)
        self.assertEqual(len(plan),len(points)+2)
        self.assertEqual(len(plan),len({r["point_id"] for r in plan}))
        self.assertAlmostEqual(result["refinements"][0]["xi"],-.003)
        self.assertAlmostEqual(result["refinements"][1]["xi"],.0006)

    def test_existing_midpoint_is_not_recomputed_under_a_new_identity(self):
        points=[point("left",-.001),point("center",0.),point("right",.001)]
        result,plan=self.run_plan(points,[edge(points[0],points[2],2.)])
        self.assertEqual(result["new_xi_points"],0)
        self.assertEqual(result["refinements"],[])
        self.assertEqual(len(plan),3)
        self.assertEqual(set(result["jobs"][0]["point_ids"]),{"left","center","right"})


if __name__ == "__main__":
    unittest.main()
