"""Recount a manually reviewed literature-figure audit without external services.

Usage: python summarize_figure_style_audit.py AUDIT_DIRECTORY [--verify-pdfs]
The figure judgments are inputs, not classifications inferred by this script.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


LOCATIONS = (
    "inside", "top_external", "other_external", "direct_or_caption",
    "no_key_needed", "uncertain",
)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def count_figures(rows):
    counts = Counter(row["legend_location"] for row in rows)
    return {
        "figure_count": len(rows),
        "paper_count": len({row["paper_id"] for row in rows}),
        "counts": {key: counts[key] for key in LOCATIONS},
        "percent": {
            key: round(100 * counts[key] / len(rows), 4) if rows else None
            for key in LOCATIONS
        },
    }


def summarize(audit_dir, verify_pdfs=False):
    coding = read_json(audit_dir / "visual_coding.json")
    papers = read_json(audit_dir / "papers.json")
    paper_map = {paper["paper_id"]: paper for paper in papers}
    if set(paper_map) != set(coding["papers"]):
        raise ValueError("Paper coverage differs between coding and source inventory")

    rows = []
    verified = []
    for paper in papers:
        pid = paper["paper_id"]
        codes = coding["papers"][pid]
        if not codes.get("reviewed_sheets"):
            raise ValueError(f"{pid}: no visual-review record")
        numbered = sorted(int(number) for number in codes["figures"])
        expected = sorted(paper["figure_pages"], key=int)
        if numbered != [int(number) for number in expected]:
            raise ValueError(f"{pid}: incomplete numbered-figure coverage")
        if numbered != list(range(1, numbered[-1] + 1)):
            raise ValueError(f"{pid}: non-contiguous figure numbering; inspect explicitly")
        if verify_pdfs:
            digest = hashlib.sha256(Path(paper["local_pdf"]).read_bytes()).hexdigest()
            if digest != paper["pdf_sha256"]:
                raise ValueError(f"{pid}: PDF hash mismatch")
            verified.append(pid)
        for number in numbered:
            value = coding["default"] | codes["defaults"] | codes["figures"][str(number)]
            location = value["legend_location"]
            if value["eligible"] and location not in LOCATIONS:
                raise ValueError(f"{pid} F{number}: invalid eligible legend category")
            if not value["eligible"] and location != "not_applicable":
                raise ValueError(f"{pid} F{number}: exclusion must be explicit")
            if value["global_top_header"] and location != "top_external":
                raise ValueError(f"{pid} F{number}: inconsistent global top key")
            if value["top_header_rows"] and not value["global_top_header"]:
                raise ValueError(f"{pid} F{number}: header row count without top key")
            rows.append({
                "paper_id": pid,
                "stratum": paper["stratum"],
                "doi": paper["doi"],
                "arxiv": paper["arxiv"],
                "pdf_version": paper["pdf_version"],
                "pdf_sha256": paper["pdf_sha256"],
                "figure_number": number,
                "pdf_page": paper["figure_pages"][str(number)],
                **value,
                "local_pdf": paper["local_pdf"],
                "source_url": paper["source_url"],
            })

    eligible = [row for row in rows if row["eligible"]]
    by_paper = []
    for paper in papers:
        pid = paper["paper_id"]
        subset = [row for row in eligible if row["paper_id"] == pid]
        counts = count_figures(subset)
        by_paper.append({
            "paper_id": pid,
            "stratum": paper["stratum"],
            "journal": paper["journal"],
            "year": paper["year"],
            "doi": paper["doi"],
            "numbered_figures": len(coding["papers"][pid]["figures"]),
            "eligible_figures": len(subset),
            **counts["counts"],
            "uses_top_external": counts["counts"]["top_external"] > 0,
        })
    subsets = {
        "all": eligible,
        "core": [row for row in eligible if row["stratum"] == "core"],
        "adjacent": [row for row in eligible if row["stratum"] == "adjacent"],
        "publisher_only": [row for row in eligible if row["pdf_version"] == "publisher_open_aps"],
        "author_versions_only": [row for row in eligible if row["pdf_version"] != "publisher_open_aps"],
        "prd_only": [row for row in eligible if paper_map[row["paper_id"]]["journal"] == "Phys.Rev.D"],
        "multipanel": [row for row in eligible if row["panel_count"] > 1],
        "transport_specific": [row for row in eligible if row["transport_specific"]],
        "explicit_series_keys": [row for row in eligible if row["legend_location"] in LOCATIONS[:3]],
    }
    summary = {
        "schema_version": 1,
        "review_date": coding["review_date"],
        "unit": coding["unit"],
        "numbered_figures": len(rows),
        "eligible_figures": len(eligible),
        "exclusions": [
            {key: row[key] for key in ("paper_id", "figure_number", "pdf_page", "note")}
            for row in rows if not row["eligible"]
        ],
        "subsets": {name: count_figures(subset) for name, subset in subsets.items()},
        "paper_level_ever_used": {
            key: sum(paper[key] > 0 for paper in by_paper) for key in LOCATIONS
        },
        "global_top_keys": [
            {key: row[key] for key in (
                "paper_id", "figure_number", "pdf_page", "panel_count", "top_header_rows"
            )}
            for row in eligible if row["global_top_header"]
        ],
        "two_row_global_top_key_count": sum(
            row["global_top_header"] and row["top_header_rows"] == 2 for row in eligible
        ),
        "panel_count_distribution": dict(sorted(Counter(
            row["panel_count"] for row in eligible
        ).items())),
        "figures_with_insets": [
            f"{row['paper_id']} F{row['figure_number']}"
            for row in eligible if "inset" in row["detail_view"]
        ],
        "pdfs_verified_this_run": verified,
        "limitations": [
            "Fixed targeted sample; not a random or exhaustive field census.",
            "301 non-selected records lack contemporaneous per-record exclusion reasons.",
            "Figures are clustered by paper and overlapping author groups; no population confidence interval.",
            "13 publisher PDFs and 3 author versions; version sensitivity is reported separately.",
            "One reviewing agent; dimensions, exact source fonts and grayscale quality are not inferred from screenshots.",
            "No sampled figure has exactly 9 or 12 primary panels, as in the target v12 composites.",
        ],
    }
    write_csv(audit_dir / "figure_audit.csv", rows)
    write_csv(audit_dir / "paper_summary.csv", by_paper)
    (audit_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "papers": len(papers), "numbered_figures": len(rows),
        "eligible_figures": len(eligible),
        "counts": summary["subsets"]["all"]["counts"],
        "pdf_hashes_verified": len(verified),
    }, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit_dir", type=Path)
    parser.add_argument("--verify-pdfs", action="store_true")
    args = parser.parse_args()
    summarize(args.audit_dir.resolve(), args.verify_pdfs)
