import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
POINTER_PATH = REPO_ROOT / "docs" / "analysis" / "pnjl" / "raw_curve_archive_v1" / "external_archive_pointer.json"
README_PATH = REPO_ROOT / "docs" / "analysis" / "pnjl" / "raw_curve_archive_v1" / "README.md"
RESTORE_WORKFLOW_PATH = REPO_ROOT / ".github" / "workflows" / "pnjl-raw-curve-archive-zenodo-restore.yml"


def test_external_pointer_is_published_and_matches_archive_contract():
    pointer = json.loads(POINTER_PATH.read_text(encoding="utf-8"))
    archive = pointer["archive"]
    provenance = pointer["provenance"]

    assert pointer["schema"] == "raw_curve_archive_v1.external_archive_pointer_v1"
    assert pointer["status"] == "canonical_external_archive"
    assert pointer["record"]["doi"] == "10.5281/zenodo.21980679"
    assert pointer["record"]["record_id"] == "21980679"
    assert archive["size_bytes"] == 1_879_467_478
    assert archive["sha256"] == "467be7fb1075d1a5f0de3dd0d8afe29d9206a156c0ca7135a1e50967a4f18ccc"
    assert archive["sha256_sidecar"]["size_bytes"] == 119
    assert archive["sha256_sidecar"]["content"].startswith(archive["sha256"] + "  ")
    assert archive["inner_manifest"]["archive_manifest_sha256"] == "514d9a7dd4cf537e8b209ed7df1cb996f52da48ab0b3672f27c3437d0cba4e52"
    assert archive["inner_manifest"]["representative_count"] == 279
    assert provenance["method"] == "independent_oracle"
    assert provenance["validation_status"] == "full_domain_verified"
    assert provenance["curve_count"] == 10_458
    assert provenance["calculation_sha"] == "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
    assert provenance["postprocess_workflow_sha"] == "67d73f871578e35759c08b3c75200c51646cf6cd"
    assert "representative-only data copy" in pointer["retention"]["representative_curves"]


def test_readme_and_restore_workflow_reference_canonical_pointer():
    readme = README_PATH.read_text(encoding="utf-8")
    workflow = RESTORE_WORKFLOW_PATH.read_text(encoding="utf-8")

    assert "10.5281/zenodo.21980679" in readme
    assert "external_archive_pointer.json" in readme
    assert "--require-full-domain" in workflow
    assert "sha256sum --strict --check" in workflow
    assert "restore-curve" in workflow
