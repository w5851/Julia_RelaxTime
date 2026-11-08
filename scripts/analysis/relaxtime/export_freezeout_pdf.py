"""Export accepted freezeout PNGs with their unchanged renderers and inputs."""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
import tempfile
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.analysis.relaxtime import plot_freezeout_comparison as comparison
from scripts.analysis.relaxtime import plot_freezeout_envelope as envelope
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import generator_record, input_record, sha256_file, write_manifest
from scripts.plotting.plot_provenance import validate_hash_record
from scripts.plotting.plot_quality import export_figure, placement_limits
from scripts.plotting.validate_plot_artifact import validate_manifest


class RenderCaptured(Exception):
    """Stop the review-only renderer at its export boundary, before any writes."""


def resolve(record):
    path = Path(record["path"])
    return path if path.is_absolute() else ROOT / path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_for(manifest, role, name=None):
    paths = [resolve(item) for item in manifest["inputs"] if item["role"] == role
             and (name is None or Path(item["path"]).name == name)]
    require(len(paths) == 1, f"expected one {role} input: {name}")
    return paths[0]


def frozen_directory(directory, names, archive_name, temporary):
    """Stage archived dependencies for the unchanged renderer's hash checks."""
    if all((directory / name).is_file() for name in names):
        return directory
    archive = ROOT / "data/outputs/results/relaxtime/meson_density/freezeout_review_20261022" / archive_name
    staged = temporary / directory.name
    with zipfile.ZipFile(archive) as bundle:
        for name in names:
            target = (staged / name).resolve()
            require(target.is_relative_to(staged.resolve()), "invalid calculation member path")
            target.parent.mkdir(parents=True, exist_ok=True)
            original = directory / name
            if original.is_file():
                shutil.copyfile(original, target)
            else:
                target.write_bytes(bundle.read(original.relative_to(ROOT).as_posix()))
    # The original loaders verify restored bytes against the frozen manifests.
    return staged


def deliver(acceptance_path, output):
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
    require(acceptance.get("schema") == "plot_png_acceptance_v1"
            and acceptance.get("status") == "author_accepted"
            and acceptance.get("pdf_export_authorized") is True
            and bool(acceptance.get("author_instruction", "").strip()), "actual PNG acceptance is required")
    for name in ("accepted_manifest", "accepted_png"):
        require(not validate_hash_record(acceptance[name], root=ROOT, label=name, allow_historical=False),
                f"accepted artifact changed: {name}")
    source = resolve(acceptance["accepted_manifest"])
    accepted = json.loads(source.read_text(encoding="utf-8"))
    require(not validate_manifest(source, repo_root=ROOT), f"invalid accepted PNG manifest: {source}")
    require(sha256_file(resolve(accepted["generator"])) == accepted["generator"]["sha256"],
            "accepted renderer changed")
    experiment = source_for(accepted, "experimental_observations")
    is_defaults = Path(accepted["generator"]["path"]).name == "plot_freezeout_comparison.py"
    require(is_defaults or Path(accepted["generator"]["path"]).name == "plot_freezeout_envelope.py",
            "unsupported freezeout renderer")
    if is_defaults:
        stem = "kplus_piplus_freezeout_defaults"
    else:
        stem = "kplus_piplus_" + ("andronic" if "andronic" in accepted["asset_id"] else "polynomial") + "_envelope"
    captured = {}
    with tempfile.TemporaryDirectory(prefix="freezeout-vector-") as directory:
        temporary = Path(directory).resolve()

        def capture(figure, ignored_stem, profile):
            pngs, _ = export_figure(figure, temporary / "rerender", profile, formats=["png"])
            require(pngs[0]["sha256"] == acceptance["accepted_png"]["sha256"],
                    f"PNG rerender differs from accepted bytes: {stem}")
            pdfs, quality = export_figure(figure, output / stem, profile, formats=["pdf"])
            require(quality["actual_fonts"] == accepted["rendering"]["quality"]["actual_fonts"],
                    "measured font identity changed")
            captured.update(pdfs=pdfs, quality=quality, profile=profile)
            raise RenderCaptured

        try:
            if is_defaults:
                run_dir = source_for(accepted, "direct_q0_calculations").parent
                outputs = comparison.checked_json(run_dir / "manifest.json")["outputs"]
                names = ["run.json", "run.json.sha256", "manifest.json", "manifest.json.sha256", *outputs]
                available = frozen_directory(run_dir, names, "freezeout_defaults_20261017_raw.zip", temporary)
                models, calculation, state = comparison.load_calculation(available)
                with patch.object(comparison, "export_review_png", capture):
                    comparison.render(models, comparison.load_experiment(experiment), temporary,
                        temporary, calculation, state, run_dir, experiment)
            else:
                summary_dir = source_for(accepted, "frozen_exploration", "summary.json").parent
                checksums = json.loads((summary_dir / "checksums.json").read_text(encoding="utf-8"))
                archive = ("andronic" if "andronic" in accepted["asset_id"] else "freezeout") + "_envelope_20261018_raw.zip"
                available = frozen_directory(summary_dir, ["checksums.json", *checksums], archive, temporary)
                args = [envelope.__file__, "--summary-dir", str(available), "--experiment", str(experiment),
                        "--output-dir", str(temporary)]
                with patch.object(sys, "argv", args), patch.object(envelope, "export_review_png", capture), \
                        patch.object(envelope, "prepare_review_directories", lambda *args, **kwargs: None):
                    envelope.main()
        except RenderCaptured:
            pass
        finally:
            plt.close("all")
    require(bool(captured), "renderer did not reach its export boundary")
    record = deepcopy(accepted)
    record.update(asset_id=stem, case_slug=output.name, generated_at_utc=datetime.now(timezone.utc).isoformat(),
        outputs=[*accepted["outputs"], *captured["pdfs"]], delivery_stage="vector_delivery",
        vector_delivery_pending=False, publication_scope="qualitative_meson_manuscript",
        author_review=dict(status="author_accepted_png",
                           record=input_record(acceptance_path, role="author_png_acceptance")))
    record["validation"].update(author_review="accepted", png_rerender_byte_identical=True)
    record["rendering"].update(delivery_stage="vector_delivery", vector_delivery_pending=False,
        quality=captured["quality"], output_formats=["png", "pdf"],
        placement_limits=placement_limits(captured["quality"], captured["profile"], record["outputs"]))
    return dict(stem=stem, accepted_manifest=acceptance["accepted_manifest"]["path"]), record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acceptance", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"delivery output already exists: {output}")
    output.mkdir(parents=True)
    charts, records = [], []
    for path in args.acceptance:
        chart, record = deliver(path.resolve(), output)
        charts.append(chart)
        records.append(record)
        print(f"Verified identical PNG and exported PDF: {chart['stem']}", flush=True)
    package = dict(case_slug=output.name, charts=charts, delivery_stage="vector_delivery",
        vector_delivery_pending=False, png_storage="reference_accepted_png_cases",
        vector_exporter=generator_record(__file__, command=" ".join(sys.argv)),
        delivery_validator=input_record(ROOT / "scripts/plotting/plot_delivery.py", role="vector_delivery_validation"))
    manifest = output / "plot_manifest.json"
    write_manifest(manifest, build_bundle(package, records))
    issues = validate_manifest(manifest, repo_root=ROOT)
    require(not issues, "; ".join(issues))
    print(json.dumps(dict(manifest=str(manifest), vector_figures=len(records), validation="passed")))


if __name__ == "__main__":
    main()
