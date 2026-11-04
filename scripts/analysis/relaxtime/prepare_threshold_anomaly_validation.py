"""Bind the mode-A anomaly study to immutable original production inputs."""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def key(row):
    return tuple(round(float(row[k]), 8) for k in ("T_MeV", "muB_MeV", "xi"))


def main(plot_manifest: Path, out: Path):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new empty output directory")
    repo = Path(__file__).resolve().parents[3]
    manifest = json.loads(plot_manifest.read_text(encoding="utf-8"))
    scan_path = Path(manifest["source_csv"])
    assert sha(scan_path) == manifest["source_csv_sha256"]
    source_plot = Path(manifest["source_plot_manifest"])
    assert sha(source_plot) == manifest["source_plot_manifest_sha256"]
    diag_path = scan_path.with_name("mode_a_diag.csv")
    scan, diag = read_rows(scan_path), read_rows(diag_path)
    diag_index = collections.defaultdict(dict)
    for r in diag:
        diag_index[key(r)][r["species"], r["channel"]] = r
    audit_path = repo / "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v5/tables/display_smoothing_window_audit_v2.csv"
    windows = []
    for r in read_rows(audit_path):
        if r["mode_key"] != "mode_a" or r["action"] != "apply_all_observables":
            continue
        windows.append(dict(window_id=r["window_id"], panel=r["plot_panel"], series=r["plot_series"],
            xis=[float(r[k]) for k in ("left_xi", "xi", "right_xi")],
            selection="previously smoothed isolated/rapid-change window", control=False))
    # These are same-branch residual windows, not the first-order crossing at xi=0.
    extras = [
        ("mu900_a1_strange_and_branch", "muB900.0", "alpha1.0", [-.07,-.06,-.05,-.04,-.03,-.02,-.01]),
        ("mu900_a11_rise", "muB900.0", "alpha1.1", [.16,.17,.18,.19,.20,.21,.22]),
        ("mu900_a11_late_rise", "muB900.0", "alpha1.1", [.41,.42,.43,.44,.45,.46]),
        ("mu900_a12_endpoint", "muB900.0", "alpha1.2", [.47,.48,.49,.50]),
        ("mu0_a1_rapid_control", "muB0.0", "alpha1.0", [-.24,-.22,-.20,-.18]),
    ]
    for name, panel, series, xis in extras:
        windows.append(dict(window_id=name,panel=panel,series=series,xis=xis,
            selection="same-branch residual or rapid-change control",control=name.endswith("control")))
    chosen, memberships, selected_channels = {}, collections.defaultdict(list), collections.defaultdict(set)
    for window in windows:
        selected = [r for r in scan if r["plot_panel"] == window["panel"] and r["plot_series"] == window["series"]
                    and any(abs(float(r["xi"])-x) < 1e-8 for x in window["xis"])]
        selected.sort(key=lambda r: float(r["xi"]))
        assert len(selected) == len(window["xis"]), window
        assert not (selected[0]["phase_reference_kind"] == "first_order" and
                    len({r["phase_curr"] for r in selected}) > 1), "Do not bridge first-order phases"
        scores = collections.defaultdict(float)
        for a,b,c in zip(selected,selected[1:],selected[2:]):
            ta,tb,tc = (diag_index[key(r)] for r in (a,b,c))
            fraction = (float(b["xi"])-float(a["xi"]))/ (float(c["xi"])-float(a["xi"]))
            for q,d in tb.items():
                predicted = float(ta[q]["contribution"])*(1-fraction)+float(tc[q]["contribution"])*fraction
                score = abs(float(d["contribution"])-predicted)/float(d["total"])
                scores[q[1]] = max(scores[q[1]],score)
        ranked = sorted(scores,key=scores.get,reverse=True)
        channels = sorted(set(ranked[:2]) | {c for c in ranked if scores[c] >= 5e-4})
        window["channels"] = channels
        window["channel_residual_scores"] = dict(scores)
        window["phase_branch"] = sorted({r["phase_curr"] for r in selected})
        window["selection_threshold"] = "max absolute channel-contribution curvature / species total >= 5e-4; at least top two"
        for r in selected:
            k=key(r)
            chosen[k]=r
            memberships[k].append(window["window_id"])
            selected_channels[k].update(channels)
    points, plans = [], []
    for k,r in sorted(chosen.items()):
        point_id=f"{r['scan_group']}_xi{float(r['xi']):+.4f}"
        points.append(dict(r,point_id=point_id,review_mode="mode_a",window_ids=";".join(memberships[k])))
        for d in diag_index[k].values():
            plans.append(dict(d,point_id=point_id,recompute=str(d["channel"] in selected_channels[k]).lower()))
    out.mkdir(parents=True,exist_ok=True)
    save_rows(out/"review_points.csv",points)
    save_rows(out/"channel_plan.csv",plans)
    report=dict(status="analysis_only",production_rows=len(scan),windows=windows,
        unique_backgrounds=len(points),point_process_jobs=sum(len(v) for v in selected_channels.values()),
        first_order_exclusion="No first-order branch crossing or xi=0 coexistence point is changed",
        untouched_channels="Retained explicitly in channel_plan; partial channel replacement, not full production re-run",
        sources=[dict(path=str(p.resolve()),sha256=sha(p)) for p in (plot_manifest,source_plot,scan_path,diag_path,audit_path)],
        outputs=[dict(path=p.name,sha256=sha(p)) for p in (out/"review_points.csv",out/"channel_plan.csv")])
    (out/"selection.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("SELECTED",len(windows),"windows",len(points),"backgrounds",report["point_process_jobs"],"point/process jobs")
    for w in windows:print(w["window_id"],w["xis"],w["channels"])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plot_manifest",type=Path)
    parser.add_argument("new_output_directory",type=Path)
    args=parser.parse_args()
    main(args.plot_manifest,args.new_output_directory)
