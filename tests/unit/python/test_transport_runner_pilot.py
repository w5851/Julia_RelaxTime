import copy
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("runner_pilot", ROOT / "scripts/analysis/relaxtime/review_transport_runner_pilot.py")
pilot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pilot)


class RunnerPilotTests(unittest.TestCase):
    def test_gnu_time_handles_minutes_and_hours(self):
        template = """User time (seconds): 61.5
System time (seconds): 0.3
Elapsed (wall clock) time (h:mm:ss or m:ss): {}
Maximum resident set size (kbytes): 1048576
Exit status: 0
"""
        self.assertEqual(pilot.parse_process_time(template.format("1:02.34"))["wall_seconds"], 62.34)
        result = pilot.parse_process_time(template.format("2:01:02.34"))
        self.assertEqual(result["wall_seconds"], 7262.34)
        self.assertEqual(result["peak_rss_bytes"], 1024**3)
        for text in (template.replace("Exit status: 0", "Exit status: 1"),
                     template.replace("User time (seconds): 61.5", "User time (seconds): nan"),
                     template.replace("Maximum resident set size (kbytes): 1048576\n", "")):
            with self.subTest(text=text), self.assertRaises(ValueError):
                pilot.parse_process_time(text.format("1:02.34"))

    def test_exact_dispatch_point_set_including_phase_anchor(self):
        job = dict(expected_points=2, expected_T_MeV=200.0898076157829,
                   inputs=dict(mode="fixed-muB-phase-scaled", muB_list="0", xi_list="-0.22,0.37"))
        rows = [dict(T_MeV="200.089807616", muB_MeV="0", xi=str(xi)) for xi in (-0.22, 0.37)]
        pilot.check_point_set(rows, job)
        for bad in (rows[:1], [rows[0], rows[0]], [dict(rows[0], T_MeV="200.1"), rows[1]],
                    [dict(rows[0], xi="-0.23"), rows[1]], [dict(rows[0], T_MeV="nan"), rows[1]]):
            with self.subTest(rows=bad), self.assertRaises(ValueError):
                pilot.check_point_set(bad, job)

    def test_cost_does_not_divide_multipt_job_or_double_count_jit(self):
        runs = [dict(completed_points=2, scan=dict(elapsed_seconds=300, compile_seconds=120),
                     job_seconds=420, process=dict(peak_rss_bytes=1024**3), total_memory_bytes=8*1024**3)]
        result = pilot.calibrate_cost(runs)
        self.assertEqual(result["per_point_budget_seconds"], 600)
        self.assertEqual(result["reserved_overhead_seconds"], 240)
        self.assertEqual(result["max_points_per_shard"], 17)
        self.assertIsNone(result["steady_point_seconds"])
        runs[0]["completed_points"] = 10
        self.assertEqual(pilot.calibrate_cost(runs), result)

    def test_cost_rejects_unusable_memory_or_time(self):
        base = dict(scan=dict(elapsed_seconds=100), job_seconds=110,
                    process=dict(peak_rss_bytes=1), total_memory_bytes=10)
        for run in (dict(base, process=dict(peak_rss_bytes=9)),
                    dict(base, scan=dict(elapsed_seconds=6000), job_seconds=6010),
                    dict(base, scan=dict(elapsed_seconds=float("nan")))):
            with self.subTest(run=run), self.assertRaises(ValueError):
                pilot.calibrate_cost([run])

    def test_archive_rejects_tampering_and_traversal(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive, unpacked = root / "archive.zip", root / "extracted"
            unpacked.mkdir()
            (unpacked / "result.csv").write_bytes(b"original\n")
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.writestr("result.csv", b"original\n")
            checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
            pilot.verified_zip(archive, unpacked, checksum)
            (unpacked / "result.csv").write_bytes(b"changed\n")
            with self.assertRaisesRegex(ValueError, "differs"):
                pilot.verified_zip(archive, unpacked, checksum)
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.writestr("../outside.csv", b"data\n")
            with self.assertRaisesRegex(ValueError, "leaves"):
                pilot.verified_zip(archive, unpacked, hashlib.sha256(archive.read_bytes()).hexdigest())

    def test_timing_identity_and_interval_guards(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for relative in ("scripts/relaxtime/run_phase_guided_transport_scan.jl",
                             "benchmark/relaxtime/transport_runner_timing.jl"):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fixture\n")
            checksum = hashlib.sha256(b"fixture\n").hexdigest()
            args = dict(mode="fixed-T-sparse-muB", case_name="fixture", sampling_profile="thermal_weighted_v1",
                        xi_list="-0.5,0.5", muB_list="0", t_list="120", calculation_sha="a"*40)
            phase = dict(elapsed_seconds=20.0, gc_seconds=0.1, compile_seconds=5.0,
                         recompile_seconds=0.0, allocated_bytes=1000, compile_timing_available=True)
            timing = dict(entrypoint_returned=True, julia_version="1.12.5", julia_threads=1,
                          blas_threads_at_start=1, blas_threads_after_load=1, blas_threads_at_end=1,
                          calculation_sha="a"*40, workflow_sha="a"*40, github_run_id="123",
                          argv=pilot.expected_argv(args), entrypoint_sha256=checksum,
                          timing_script_sha256=checksum, compile_time_included_in_elapsed=True,
                          load=phase, scan=phase, measured_wrapper_seconds=40.1,
                          wrapper_initialization_seconds=0.3, total_memory_bytes=8000)
            process = dict(wall_seconds=41.0, peak_rss_bytes=2000)
            pilot.check_timing(timing, process, args, root, 123, "a"*40)
            for field, value in (("workflow_sha", "b"*40), ("blas_threads_after_load", 2),
                                 ("argv", timing["argv"] + ["--no-compute-bulk"]),
                                 ("measured_wrapper_seconds", 39.0), ("julia_version", "1.11.6")):
                altered = copy.deepcopy(timing)
                altered[field] = value
                with self.subTest(field=field), self.assertRaises(ValueError):
                    pilot.check_timing(altered, process, args, root, 123, "a"*40)


if __name__ == "__main__":
    unittest.main()
