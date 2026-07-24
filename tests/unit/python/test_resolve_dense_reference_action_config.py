import importlib.util
from pathlib import Path

import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "pnjl" / "resolve_dense_reference_action_config.py"
WORKFLOW_PATH = PROJECT_ROOT / ".github" / "workflows" / "pnjl-dense-reference.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("dense_action_config", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_resolves_whitelisted_controls_in_deterministic_order():
    module = load_module()
    resolved = module.resolve_action_config(
        '{"T_refine_levels":3,"thermo_quadrature_rtol":1e-10,"crossover_T_max":240}'
    )
    assert resolved == [
        "--T-refine-levels",
        "3",
        "--thermo-quadrature-rtol",
        "1e-10",
        "--crossover-T-max",
        "240",
    ]


@pytest.mark.parametrize(
    "payload",
    [
        '[]',
        '{"unknown":1}',
        '{"T_refine_levels":1.5}',
        '{"thermo_quadrature_rtol":0}',
        '{"thermo_quadrature_maxevals":true}',
    ],
)
def test_rejects_invalid_or_unsupported_controls(payload):
    module = load_module()
    with pytest.raises(ValueError):
        module.resolve_action_config(payload)


def test_workflow_dispatch_stays_within_github_property_limit():
    workflow = yaml.load(WORKFLOW_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    inputs = workflow["on"]["workflow_dispatch"]["inputs"]
    assert len(inputs) == 25
    assert len(inputs) <= 25
    assert "advanced_config_json" in inputs
    assert "calculation_ref" in inputs
