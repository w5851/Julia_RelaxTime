"""Retain the existing anomaly cases and add isotropic/domain/phase controls."""
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path


def main(input_dir: Path, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Use a new empty output directory")
    selection = json.loads((input_dir / "selection.json").read_text(encoding="utf-8"))
    source = next(s for s in selection["sources"] if Path(s["path"]).name == "mode_a_scan.csv")
    raw_path = Path(source["path"])
    if hashlib.sha256(raw_path.read_bytes()).hexdigest() != source["sha256"]:
        raise ValueError("Production source hash differs from the retained selection")
    with (input_dir / "review_points.csv").open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        selected = list(reader)
    with raw_path.open(encoding="utf-8-sig", newline="") as handle:
        production = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    identity = lambda r: (float(r["muB_MeV"]), float(r["alpha_T"]), float(r["xi"]))
    existing = {identity(r) for r in selected}
    controls = [r for r in production if any(abs(float(r["xi"]) - x) < 1e-10 for x in (-0.5, 0.0, 0.5))]
    groups = {}
    for row in production:
        groups.setdefault((row["muB_MeV"], row["alpha_T"]), []).append(row)
    # Include the nearest retained node on either side of each first-order
    # branch change. This selects controls; it never alters sampling rules.
    for rows in groups.values():
        rows.sort(key=lambda r: float(r["xi"]))
        for left, right in zip(rows[:-1], rows[1:]):
            if left["phase_curr"] != right["phase_curr"] and "first_order" in right["phase_structure"]:
                controls.extend((left, right))
    added = []
    for row in controls:
        key = identity(row)
        if key in existing:
            continue
        if str(row["converged"]).lower() != "true":
            raise ValueError(f"Control background did not converge: {key}")
        mu, alpha, xi = key
        row = {**row, "point_id": f"muB{mu:.1f}_alpha{alpha:.1f}_xi{xi:+.4f}",
               "review_mode": "mode_a", "window_ids": "independent_convergence_control"}
        selected.append(row)
        added.append(row["point_id"])
        existing.add(key)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "review_points.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sorted(selected, key=identity))
    shutil.copyfile(input_dir / "channel_plan.csv", output_dir / "channel_plan.csv")
    (output_dir / "selection.json").write_text(json.dumps({
        "status": "analysis_only", "source": source, "source_selection": str(input_dir.resolve()),
        "retained_anomaly_points": len(selected) - len(added), "added_controls": added,
        "total_points": len(selected), "channel_plan_scope": "original anomaly cases only",
        "purpose": "independent and joint numerical convergence; production data unchanged",
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"total_points": len(selected), "added_controls": len(added)}))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
