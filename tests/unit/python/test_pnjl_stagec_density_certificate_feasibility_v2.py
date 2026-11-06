import csv
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PLOTTER = ROOT / "scripts" / "analysis" / "plot_pnjl_stagec_density_certificate_feasibility_v2.py"


def test_v2_plotter_emits_figures_and_refreshes_manifest(tmp_path):
    spec = importlib.util.spec_from_file_location("stagec_v2_plotter", PLOTTER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = tmp_path / "source" / "job"
    source.mkdir(parents=True)
    rows = [
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "0", "muq_MeV": "300"},
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "1", "muq_MeV": "301"},
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "2", "muq_MeV": "300"},
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "3", "muq_MeV": "302"},
    ]
    with (source / "curve_points.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    with (source / "slice_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["method", "xi", "T_MeV", "mu_transition_MeV"])
        writer.writeheader()
        writer.writerow({"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "mu_transition_MeV": "300.5"})
    output = tmp_path / "output"
    output.mkdir()
    (output / "manifest.json").write_text(json.dumps({"schema_version": "cep_hybrid_stagec_density_certificate_feasibility_v2"}), encoding="utf-8")
    (output / "route_comparison.csv").write_text("route,cap,xi,T_MeV,maxwell_area_gate\nstage_b_features_v1,12,-0.35,51.0,1e-6\n", encoding="utf-8")
    (output / "selected_point_index.csv").write_text("xi,T_MeV,route,cap,rank,feature,rho,muq_MeV\n", encoding="utf-8")
    result = module.render(tmp_path / "source", output, output)
    assert len(result["figures"]) == 1
    figure = output / result["figures"][0]["path"]
    assert figure.is_file() and figure.stat().st_size > 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert "plot_manifest_sha256" in manifest
