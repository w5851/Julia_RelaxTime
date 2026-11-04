#!/usr/bin/env python3
"""quote-remote only: side-local differences and propagated quadrature defects."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def second_difference(rows, distance, side, fraction, field="mixed_square_rate_fm2"):
    """Use a stencil entirely on one side; distance and derivative step differ."""
    step = fraction * distance
    coordinates = [side * distance - step, side * distance, side * distance + step]
    selected = []
    for x in coordinates:
        row = min(rows, key=lambda r: abs(float(r["delta_xi"]) - x))
        if abs(float(row["delta_xi"]) - x) > 1e-12 * distance:
            raise ValueError(f"Missing side-stencil point {x}")
        selected.append(row)
    values = [float(r[field]) for r in selected]
    weights = [1.0, -2.0, 1.0]
    value = math.fsum(w * y for w, y in zip(weights, values)) / step**2
    error = math.fsum(abs(w) * float(r["quadrature_error_fm2"])
                      for w, r in zip(weights, selected)) / step**2
    return dict(curvature=value, quadrature_defect=error,
                point_ids=[r["point_id"] for r in selected], step=step)


def derive(rows, steps, fractions=(.25, .125)):
    wide_fraction, fine_fraction = fractions
    if not 0 < fine_fraction < wide_fraction < 1:
        raise ValueError("Derivative step fractions must satisfy 0 < fine < wide < 1")
    ratio_squared = (wide_fraction / fine_fraction)**2
    result = []
    for side in (-1, 1):
        previous = None
        for h in sorted(steps, reverse=True):
            wide = second_difference(rows, h, side, wide_fraction)
            fine = second_difference(rows, h, side, fine_fraction)
            extrapolated = (ratio_squared * fine["curvature"] - wide["curvature"]) / (ratio_squared - 1)
            # This fourth-order estimate is local at u != 0, never across Mott.
            defect = abs(extrapolated - fine["curvature"])
            numerical = (ratio_squared * fine["quadrature_defect"] + wide["quadrature_defect"]) / (ratio_squared - 1)
            exponent = None if previous is None or previous[1] * extrapolated <= 0 else \
                math.log(abs(extrapolated / previous[1])) / math.log(previous[0] / h)
            result.append(dict(side="left" if side < 0 else "right", h=h, delta_xi=side*h,
                curvature_wide=wide["curvature"], curvature_fine=fine["curvature"],
                curvature=extrapolated, scaled_curvature=math.sqrt(h)*extrapolated,
                step_defect=defect, quadrature_defect=numerical,
                empirical_relative_defect=(defect+numerical)/max(abs(extrapolated),1e-300),
                effective_exponent=exponent,
                point_ids=";".join(dict.fromkeys(wide["point_ids"]+fine["point_ids"]))))
            previous = h, extrapolated
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    source = args.case / "summary.json"
    summary = json.loads(source.read_text(encoding="utf-8-sig"))
    for name, expected in summary["files"].items():
        if digest(args.case/name) != expected:
            raise ValueError("Frozen input hash mismatch: " + name)
    rows = read_csv(args.case/"integrals.csv")
    result = derive(rows, summary["steps"])
    args.output.mkdir(exist_ok=False)
    target = args.output/"side_curvatures.csv"
    with target.open("w",encoding="utf-8",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(result[0]))
        writer.writeheader();writer.writerows(result)
    record=dict(status="completed_side_curvature_diagnostic",source=str(source),
        source_sha256=digest(source),scope=summary["scope"],xi_mott=summary["xi_mott"],
        derivative="Second derivative of continuous near-threshold mixed-scalar positive-square rate with respect to xi; not full tau second derivative.",
        difference_rule="u=+/-h; centered steps h/4 and h/8 on each side; (4D_h8-D_h4)/3",
        error_semantics="Propagated Gauss-rule defects and derivative-step discrepancy, not a rigorous bound or complete background/kernel quadrature certificate.",
        asymptote_test="sqrt(h)*I_second and successive log ratios, without fitting or imposing the exponent",
        script_sha256=digest(Path(__file__)),files={target.name:digest(target)})
    (args.output/"summary.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=record["status"],rows=len(result),left_inner=result[len(result)//2-3:len(result)//2],right_inner=result[-3:]),allow_nan=False))


if __name__=="__main__":
    main()
