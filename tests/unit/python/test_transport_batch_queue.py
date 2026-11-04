import importlib.util
from pathlib import Path
import tempfile
import sys
import subprocess
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("transport_batch_queue", ROOT / "scripts/analysis/relaxtime/run_transport_batch_queue.py")
queue = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(queue)


def entry(label):
    return dict(mode="A", job=dict(label=label, inputs=dict(case_name="fixture"), expected_points=14),
                request_file=label + ".json", request_sha256="hash")


def specification(entries):
    return dict(repository="fixture/repo", workflow_id=1, workflow_sha="a" * 40,
                workflow_ref="frozen", planned_runs=entries, runs=[], max_outstanding=2)


def run_record(job, **changes):
    return dict(dict(id=12, display_title=queue.title(job), head_sha="a" * 40,
                     head_branch="frozen", event="workflow_dispatch", run_attempt=1,
                     html_url="https://example.invalid/12", status="in_progress", conclusion=None), **changes)


class TransportBatchQueueTests(unittest.TestCase):
    def test_main_recovers_a_read_timeout_without_restarting_or_posting(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            spec = dict(specification([entry("job")]), frozen_files={})
            queue.write(root / "queue_plan.json", spec)
            queue.write(root / "dispatch_manifest.json", spec)
            ready = dict(all_jobs_terminal=True, ready_for_remote_review=True,
                         collection_errors=[], runs=[dict(conclusion="success")])
            failure = subprocess.CalledProcessError(1, ["gh", "api", "read-only-fixture"])
            with patch.object(sys, "argv", ["queue", str(root)]), \
                 patch.object(queue, "api", return_value=dict(object=dict(sha=spec["workflow_sha"]))), \
                 patch.object(queue, "cycle", side_effect=[failure, ready]) as cycle, \
                 patch.object(queue, "dispatch") as post, patch.object(queue.time, "sleep"):
                queue.main()
            self.assertEqual(cycle.call_count, 2)
            post.assert_not_called()
            recovered = queue.read(root / "network_retry_status.json")
            self.assertFalse(recovered["recovering"])
            self.assertIn("read-only-fixture", recovered["error"])

    def test_main_waits_for_delayed_artifact_after_successful_terminal_run(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            spec = dict(specification([entry("job")]), frozen_files={})
            queue.write(root / "queue_plan.json", spec)
            queue.write(root / "dispatch_manifest.json", spec)
            pending = dict(all_jobs_terminal=True, ready_for_remote_review=False,
                           collection_errors=[], runs=[dict(conclusion="success")])
            ready = dict(pending, ready_for_remote_review=True)
            with patch.object(sys, "argv", ["queue", str(root)]), \
                 patch.object(queue, "api", return_value=dict(object=dict(sha=spec["workflow_sha"]))), \
                 patch.object(queue, "cycle", side_effect=[pending, ready]) as cycle, \
                 patch.object(queue.time, "sleep"):
                queue.main()
            self.assertEqual(cycle.call_count, 2)

    def test_pagination_recovers_more_than_one_hundred_runs(self):
        calls = []
        def get(endpoint):
            calls.append(endpoint)
            page = int(endpoint.split("page=")[-1])
            return dict(workflow_runs=[dict(id=i) for i in (range(100) if page == 1 else range(100, 142))])
        result = queue.list_runs(specification([]), get=get)
        self.assertEqual(len(result), 142)
        self.assertEqual(len(calls), 2)
        self.assertTrue(calls[-1].endswith("page=2"))

    def test_recover_requires_unique_title_ref_sha_and_attempt(self):
        job = entry("job")
        spec = specification([job])
        run = run_record(job)
        self.assertEqual(queue.recover_entry(job, spec, [run])["run_id"], 12)
        for records in ([run, dict(run, id=13)], [dict(run, head_branch="moved")],
                        [dict(run, head_sha="b" * 40)], [dict(run, run_attempt=2)]):
            with self.subTest(records=records), self.assertRaises(ValueError):
                queue.recover_entry(job, spec, records)

    def test_saved_attempt_is_not_reposted_and_counts_against_cap(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            jobs = [entry(str(i)) for i in range(4)]
            queue.retain(root / "attempts/0.json", dict(request_may_have_been_accepted=True))
            with patch.object(queue, "list_runs", return_value=[]), patch.object(queue, "dispatch") as post, \
                 patch.object(queue.time, "sleep"):
                status = queue.cycle(root, specification(jobs))
            self.assertEqual(post.call_count, 1)
            self.assertEqual(post.call_args.args[1]["job"]["label"], "1")
            self.assertEqual(status["unresolved_attempts"], ["0"])
            self.assertFalse(status["ready_for_remote_review"])

    def test_accepted_but_unrecorded_run_is_recovered_without_post(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            job = entry("job")
            with patch.object(queue, "list_runs", return_value=[run_record(job)]), \
                 patch.object(queue, "dispatch") as post:
                status = queue.cycle(root, specification([job]))
            post.assert_not_called()
            self.assertEqual(status["recovered_jobs"], 1)
            self.assertEqual(queue.read(root / "dispatch_manifest.json")["runs"][0]["run_id"], 12)

    def test_completed_workflow_needs_successful_artifact_collection(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            job = entry("job")
            for conclusion, available, ready in (("success", True, True), ("failure", True, False),
                                                  ("success", False, False)):
                with patch.object(queue, "list_runs", return_value=[run_record(job, status="completed", conclusion=conclusion)]), \
                     patch.object(queue, "collect_run", return_value=dict(artifact_available=available)):
                    status = queue.cycle(root, specification([job]))
                self.assertEqual(status["ready_for_remote_review"], ready)

    def test_extract_preserves_bytes_and_rejects_tampering_and_traversal(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive, unpacked = root / "artifact.zip", root / "artifact"
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.writestr("data/a.csv", b"raw bytes \n")
            self.assertEqual(queue.extract_verified(archive, unpacked), 1)
            self.assertEqual(queue.extract_verified(archive, unpacked), 1)
            (unpacked / "data/a.csv").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "differs"):
                queue.extract_verified(archive, unpacked)
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.writestr("../outside.csv", b"escape")
            with self.assertRaisesRegex(ValueError, "leaves"):
                queue.extract_verified(archive, unpacked)
            self.assertFalse((root / "outside.csv").exists())

    def test_existing_attempt_prevents_post_before_any_network_call(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            job = entry("job")
            spec = specification([job])
            queue.write(root / "job.json", dict(ref="frozen", inputs=job["job"]["inputs"]))
            job["request_sha256"] = queue.digest(root / "job.json")
            queue.retain(root / "attempts/job.json", dict(request_may_have_been_accepted=True))
            with patch.object(queue.subprocess, "run") as post, self.assertRaisesRegex(ValueError, "already exists"):
                queue.dispatch(root, job, spec)
            post.assert_not_called()

    def test_prepared_api_boolean_preserves_false_and_rejects_changed_flag(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            job = entry("job")
            job["job"]["inputs"]["render_plots"] = False
            spec = specification([job])
            inputs = dict(job["job"]["inputs"], render_plots="false")
            queue.write(root / "job.json", dict(ref="frozen", inputs=inputs))
            job["request_sha256"] = queue.digest(root / "job.json")
            with patch.object(queue.subprocess, "run") as post:
                post.return_value.returncode = 0
                post.return_value.stdout = post.return_value.stderr = ""
                result = queue.dispatch(root, job, spec)
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(queue.read(root / "attempts/job.json")["request_may_have_been_accepted"])
            inputs["render_plots"] = "true"
            queue.write(root / "job.json", dict(ref="frozen", inputs=inputs))
            job["request_sha256"] = queue.digest(root / "job.json")
            with patch.object(queue.subprocess, "run") as post, self.assertRaisesRegex(ValueError, "differs"):
                queue.dispatch(root, job, spec)
            post.assert_not_called()


if __name__ == "__main__":
    unittest.main()
