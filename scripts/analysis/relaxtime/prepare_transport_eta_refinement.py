#!/usr/bin/env python3
"""Freeze the three accepted eta-Mott windows for a quote-remote refinement."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import tomllib


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("survey", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source_rows = args.survey / "inputs/source_rows.csv"
    report_path = args.survey / "review_v2/review.json"
    report = json.loads(report_path.read_text(encoding="utf-8-sig"))
    assert sha(source_rows) == report["input_rows_sha256"]
    wanted = {
        "mode_a:muB0.0_alpha1.0": "a_mu0",
        "mode_a:muB450.0_alpha1.0": "a_mu450",
        "mode_b:T200.0_muB0.0": "b_mu0",
    }
    targets = []
    for event in report["reverse_light_continuous_events"]:
        if event["window"] in wanted:
            targets.append(dict(curve_id=wanted[event["window"]], **event))
    assert len(targets) == 3
    (args.output / "source_rows.csv").write_bytes(source_rows.read_bytes())
    profile_path = args.source / "config/models/pnjl/transport_sampling/thermal_weighted_v1.toml"
    base = profile_path.read_text(encoding="utf-8")
    expected = tomllib.loads(base)
    profiles = {}
    for name, rtol in (("base", 1e-3), ("fine", 1e-5), ("tighter", 1e-6)):
        text = base
        if name != "base":
            assert text.count('id = "thermal_weighted_v1"') == 1
            assert text.count("rtol = 1.0e-3") == 1
            text = text.replace('id = "thermal_weighted_v1"', f'id = "eta_xi_{name}_v1"')
            text = text.replace("rtol = 1.0e-3", f"rtol = {rtol:.1e}")
            # More room to meet a stricter error target; no tolerance relaxation.
            text = text.replace("max_evaluations = 40000", "max_evaluations = 80000")
            text = text.replace("max_nodes = 12000", "max_nodes = 24000")
            text = text.replace("max_rounds = 40", "max_rounds = 50")
        parsed = tomllib.loads(text)
        for section in ("background", "a_integral", "tau_integral", "terminal_integral"):
            assert parsed[section] == expected[section]
        assert parsed["sampling"]["rtol"] == rtol
        out = args.output / (name + ".toml")
        out.write_text(text, encoding="utf-8")
        profiles[name] = dict(path=out.name, sha256=sha(out), effective=parsed)
    config = dict(schema="transport_eta_xi_refinement_inputs_v1", host="quote-remote",
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444", targets=targets,
        source_rows_sha256=sha(source_rows), previous_review_sha256=sha(report_path),
        original_profile_sha256=sha(profile_path), profiles=profiles,
        mass_offsets_inv_fm=[1e-4, 1e-6, 1e-7], xi_bisection_tolerance=1e-9,
        location_method="exact endpoint-subtracted equal-mass loop; original finite-offset controls",
        endpoint_loop_node_controls=[128, 256],
        xi_offsets=[-0.008, -0.004, -0.002, -0.001, -0.0005, -0.0002, -0.0001,
                     0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.008],
        original_controls="bracket left, bracket right, and next original point",
        limits="diagnostic only; all new backgrounds and all 17 transport rates are solved; original grid remains frozen",
        preparation_script_sha256=sha(Path(__file__)))
    (args.output / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(targets=[x["curve_id"] for x in targets], profiles=profiles), indent=2))


if __name__ == "__main__":
    main()
