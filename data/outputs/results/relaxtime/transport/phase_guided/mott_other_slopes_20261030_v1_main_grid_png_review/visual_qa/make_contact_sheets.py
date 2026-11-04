"""Create image-only visual QA sheets from the existing frozen PNG manifest."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[8]
QA = Path(__file__).resolve().parent
FIGURES = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/mott_other_slopes_20261030_v1_main_grid_png_review"
manifest_path = FIGURES / "plot_manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
font = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 18)
files = sorted(
    ROOT / output["path"]
    for figure in manifest["figures"]
    for output in figure["record"]["outputs"]
    if output["format"] == "png" and not output["path"].endswith("_gray.png")
)
assert len(files) == 39
index = dict(
    schema="image_only_contact_sheet_v1",
    created_utc=datetime.now(timezone.utc).isoformat(),
    source_manifest=str(manifest_path.relative_to(ROOT)).replace("\\", "/"),
    source_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
    source_values_modified=False,
    numerical_analysis_performed=False,
    contact_sheets=[],
    details=[],
)

for start in range(0, len(files), 6):
    selected = files[start:start + 6]
    canvas = Image.new("RGB", (2025, 1240), "white")
    draw = ImageDraw.Draw(canvas)
    for offset, path in enumerate(selected):
        col, row = offset % 3, offset // 3
        with Image.open(path) as original:
            resized = original.convert("RGB").resize((675, 600), Image.Resampling.LANCZOS)
        canvas.paste(resized, (col * 675, row * 620 + 20))
        draw.text((col * 675 + 10, row * 620), path.name, fill="black", font=font)
    target = QA / f"overview_{start // 6 + 1:02d}.png"
    if target.exists():
        raise FileExistsError(target)
    canvas.save(target)
    index["contact_sheets"].append(dict(path=target.name, originals=[p.name for p in selected]))

for name in (
    "other01_tau_sbar", "other03_tau_sbar", "other07_sigma_over_T",
    "other08_tau_ubar", "other08_sigma_over_T", "other01_tau_sbar_gray",
):
    path = FIGURES / (name + ".png")
    target = QA / ("detail_" + path.name)
    if target.exists():
        raise FileExistsError(target)
    with Image.open(path) as original:
        original.convert("RGB").resize((1350, 1200), Image.Resampling.LANCZOS).save(target)
    index["details"].append(dict(path=target.name, original=path.name))

(QA / "contact_sheets.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
print(json.dumps(dict(contact_sheets=len(index["contact_sheets"]), detail_images=len(index["details"]))))
