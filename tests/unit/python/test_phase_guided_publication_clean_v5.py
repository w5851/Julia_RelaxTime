import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v5.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v5", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_v5_label_map_separates_old_and_new_wording() -> None:
    rows = MODULE.label_map_rows()
    by_key = {(row["scope"], row["key"]): row for row in rows}
    assert by_key[("relaxation_time_y_axis", "tau_u")]["old_label"] == r"$\tau_u$"
    assert by_key[("relaxation_time_y_axis", "tau_u")]["new_label"] == r"$\tau_u\;[\mathrm{fm}]$"
    assert by_key[("first_order_endpoint_legend", "quark")]["old_label"] == "chiral-restored (quark) endpoint"
    assert by_key[("first_order_endpoint_legend", "quark")]["new_label"] == "chirally restored branch endpoint"
    assert by_key[("first_order_endpoint_legend", "hadron")]["old_label"] == "chiral-broken (hadron) endpoint"
    assert by_key[("first_order_endpoint_legend", "hadron")]["new_label"] == "chirally broken branch endpoint"


def test_v5_label_map_marks_display_only() -> None:
    assert all(row["canonical_data_modified"] == "False" for row in MODULE.label_map_rows())
    assert MODULE.phase_legend_label("quark") == "chirally restored branch endpoint"
    assert MODULE.phase_legend_label("hadron") == "chirally broken branch endpoint"
