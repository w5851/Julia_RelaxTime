#!/usr/bin/env python3
"""Generate sign-flip table for upstream detM_im decomposition terms."""

from __future__ import annotations

import csv
from pathlib import Path


IN_CSV = Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_term_switch_trace.csv")
OUT_CSV = Path(r"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_signflip_table.csv")
TARGET_COLS = ["detM_im", "c1_ReM00_ImM88", "c2_ImM00_ReM88", "c3_minus2ReM08ImM08", "ReM08"]


def sgn(x: float, eps: float = 1.0e-15) -> int:
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def main() -> None:
    rows: list[dict[str, str]] = []
    with IN_CSV.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows.extend(reader)

    by_xi: dict[float, list[dict[str, str]]] = {}
    for r in rows:
        xi = float(r["xi"])
        by_xi.setdefault(xi, []).append(r)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["scenario", "xi", "quantity", "ds_left", "value_left", "ds_right", "value_right", "sign_left", "sign_right"])

        for xi, sub in sorted(by_xi.items()):
            sub_sorted = sorted(sub, key=lambda r: float(r["ds"]))
            for q in TARGET_COLS:
                for i in range(1, len(sub_sorted)):
                    l = sub_sorted[i - 1]
                    r = sub_sorted[i]
                    lv = float(l[q])
                    rv = float(r[q])
                    ls = sgn(lv)
                    rs = sgn(rv)
                    if ls != rs:
                        writer.writerow([
                            "tauu_pos_uubaruubar",
                            f"{xi:.2f}",
                            q,
                            l["ds"],
                            f"{lv:.16e}",
                            r["ds"],
                            f"{rv:.16e}",
                            ls,
                            rs,
                        ])

    print(str(OUT_CSV))


if __name__ == "__main__":
    main()
