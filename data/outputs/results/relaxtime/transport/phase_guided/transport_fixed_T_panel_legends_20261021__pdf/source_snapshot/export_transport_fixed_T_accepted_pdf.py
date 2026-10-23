"""Export the author-accepted Fig.4 through its unchanged PNG renderer.

The wrapper changes delivery metadata and export formats only. It requires
identical PNG pixels/bytes, plotted samples, ranges, and legend specifications.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from datetime import datetime, timezone
import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from PIL import Image
from scripts.plotting.plot_manifest import generator_record, input_record, runtime_record, sha256_file
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.validate_plot_artifact import validate_manifest


def resolve_record(record):
    path = Path(record["path"])
    return path if path.is_absolute() else ROOT / path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-png-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--author-acceptance", required=True)
    args = parser.parse_args()
    source_path = args.accepted_png_manifest.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing existing case: {output}")
    if errors := validate_manifest(source_path, repo_root=ROOT):
        raise ValueError(f"invalid accepted PNG package: {errors}")
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if source["rendering"]["delivery_stage"] != "png_review":
        raise ValueError("expected a PNG review manifest")
    renderer_path = resolve_record(source["generator"])
    if sha256_file(renderer_path) != source["generator"]["sha256"]:
        raise ValueError("the accepted renderer has changed")
    png_record = next(record for record in source["outputs"] if record["format"] == "png")
    accepted_png = resolve_record(png_record)
    spec = importlib.util.spec_from_file_location("accepted_fixed_T_renderer", renderer_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import accepted renderer")
    renderer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(renderer)
    original_build = renderer.build_manifest
    comparison = {}

    def accepted_export(figure, stem, profile, **kwargs):
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=profile.dpi, bbox_inches=None)
        buffer.seek(0)
        with Image.open(buffer) as generated, Image.open(accepted_png) as reference:
            if generated.size != reference.size or generated.convert("RGBA").tobytes() != reference.convert("RGBA").tobytes():
                raise ValueError("rendered pixels differ from the accepted PNG")
        outputs, quality = export_figure(figure, stem, profile,
                                        intended_width_inches=kwargs["intended_width_inches"], formats=("png", "pdf"))
        new_png = next(record for record in outputs if record["format"] == "png")
        if new_png["sha256"] != png_record["sha256"]:
            raise ValueError("PNG bytes differ from the accepted PNG")
        comparison.update({"png_pixels_identical": True, "png_bytes_identical": True,
                           "png_sha256": new_png["sha256"]})
        return outputs, quality

    def accepted_build(**kwargs):
        if kwargs["validation"]["plotted_state_sha256"] != source["validation"]["plotted_state_sha256"]:
            raise ValueError("samples, ranges, markers, or styles differ from accepted PNG")
        if kwargs["rendering"]["legend_layout"] != source["rendering"]["legend_layout"]:
            raise ValueError("legend layout differs from accepted PNG")
        snapshot = output / "source_snapshot" / Path(__file__).name
        shutil.copy2(__file__, snapshot)
        kwargs["inputs"].extend([
            input_record(source_path, role="author_accepted_png_manifest"),
            input_record(accepted_png, role="author_accepted_png_output"),
            input_record(renderer_path, role="unchanged_accepted_png_renderer"),
        ])
        command = subprocess.list2cmdline([sys.executable, str(Path(__file__).resolve()),
            "--accepted-png-manifest", str(source_path), "--output-dir", str(output),
            "--author-acceptance", args.author_acceptance])
        kwargs["generator"] = generator_record(__file__, command=command, source_snapshot=snapshot,
                                                runtime=runtime_record({"matplotlib": renderer.matplotlib.__version__}))
        kwargs["rendering"].update({"delivery_stage": "vector_delivery", "output_formats": ["png", "pdf"],
                                    "vector_delivery_pending": False,
                                    "author_visual_acceptance": {"recorded_at_utc": datetime.now(timezone.utc).isoformat(),
                                                                 "statement": args.author_acceptance},
                                    "accepted_png_comparison": comparison,
                                    "in_panel_exception_basis": "author accepted the source PNG layout; identical render verified"})
        manifest = original_build(**kwargs)
        manifest["vector_export_complete"] = True
        return manifest

    argv = [str(renderer_path), "--output-dir", str(output)]
    with patch.object(sys, "argv", argv), patch.object(renderer, "export_figure", accepted_export), \
            patch.object(renderer, "build_manifest", accepted_build), redirect_stdout(io.StringIO()):
        renderer.main()
    path = output / "plot_manifest.json"
    if errors := validate_manifest(path, repo_root=ROOT):
        raise ValueError(errors)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    print(json.dumps({"manifest": str(path), "accepted_png_comparison": comparison,
                      "pdf_inspection": next(x["inspection"] for x in manifest["outputs"] if x["format"] == "pdf"),
                      "artifact_validation": "passed"}, indent=2))


if __name__ == "__main__":
    main()
