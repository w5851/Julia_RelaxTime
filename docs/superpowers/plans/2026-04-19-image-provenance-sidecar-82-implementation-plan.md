# Image Provenance Sidecar #82 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 `scripts/plot_scan_csv.py` 图像输出链路引入默认 sidecar 溯源（`*.provenance.json`）与哈希校验入口，并形成可复用 Python 通用模块。

**Architecture:** 新增 `scripts/common/provenance_image.py` 作为通用能力层，集中处理图像哈希、v1 sidecar 构建、写入与校验。`plot_scan_csv.py` 只负责收集运行上下文并在每次 `savefig` 后调用该模块，避免在脚本内部散落重复实现。通过 Python 单测覆盖“写入成功 + 篡改检测失败 + 缺字段失败”，并补充文档明确 sidecar 方案边界。

**Tech Stack:** Python 3、matplotlib 现有绘图链、JSON/sha256 标准库、Julia smoke 回归命令（integration/unit）。

---

## File Structure

- Create: `scripts/common/provenance_image.py`
  - 责任：图像 sidecar 的通用哈希/构建/写入/校验能力。
- Modify: `scripts/plot_scan_csv.py`
  - 责任：在图像保存流程接入 sidecar 默认写入，并提供校验 CLI 入口。
- Create: `tests/unit/python/test_image_provenance_sidecar.py`
  - 责任：覆盖写入与校验（通过/篡改失败/缺字段失败）契约。
- Modify: `docs/guides/scripts/README.md`
  - 责任：补充 sidecar 使用方式、校验命令与方案边界说明。

---

### Task 1: 先写失败测试（RED）锁定 v1 sidecar 契约

**Files:**
- Create: `tests/unit/python/test_image_provenance_sidecar.py`

- [ ] **Step 1: Write failing tests for sidecar contract**

```python
from pathlib import Path
import json

from scripts.common.provenance_image import (
    write_image_provenance_sidecar,
    verify_image_provenance,
)


def test_write_sidecar_contains_required_v1_fields(tmp_path: Path):
    image = tmp_path / "demo.png"
    image.write_bytes(b"fake_png_bytes")

    metadata = {
        "schema_version": "v1",
        "generated_at_utc": "2026-04-19T00:00:00Z",
        "script_path": "scripts/plot_scan_csv.py",
        "command": "python scripts/plot_scan_csv.py --mode lines ...",
        "git_commit": "deadbeef",
        "julia_version": "python-3.x",
        "input_data_hashes": [{"path": "data/a.csv", "sha256": "abc"}],
    }

    sidecar = write_image_provenance_sidecar(image, metadata)
    payload = json.loads(sidecar.read_text(encoding="utf-8"))

    for key in [
        "schema_version",
        "generated_at_utc",
        "image_path",
        "image_sha256",
        "script_path",
        "command",
        "git_commit",
        "julia_version",
        "input_data_hashes",
    ]:
        assert key in payload


def test_verify_provenance_detects_tampered_image(tmp_path: Path):
    image = tmp_path / "tamper.png"
    image.write_bytes(b"origin")

    metadata = {
        "schema_version": "v1",
        "generated_at_utc": "2026-04-19T00:00:00Z",
        "script_path": "scripts/plot_scan_csv.py",
        "command": "python scripts/plot_scan_csv.py --mode lines ...",
        "git_commit": "deadbeef",
        "julia_version": "python-3.x",
        "input_data_hashes": [],
    }
    write_image_provenance_sidecar(image, metadata)

    image.write_bytes(b"mutated")
    ok, reason = verify_image_provenance(image)
    assert ok is False
    assert "sha256" in reason.lower()
```

- [ ] **Step 2: Run tests and confirm RED**

Run: `python -m pytest tests/unit/python/test_image_provenance_sidecar.py -q`  
Expected: FAIL（`scripts.common.provenance_image` 不存在或函数未实现）。

- [ ] **Step 3: Commit RED tests**

```bash
git add tests/unit/python/test_image_provenance_sidecar.py
git commit -m "test(scripts): add failing image provenance sidecar contract tests"
```

---

### Task 2: 实现通用模块最小可用能力（GREEN）

**Files:**
- Create: `scripts/common/provenance_image.py`
- Test: `tests/unit/python/test_image_provenance_sidecar.py`

- [ ] **Step 1: Implement minimal provenance module**

```python
from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
from typing import Dict, Tuple


REQUIRED_V1_FIELDS = [
    "schema_version",
    "generated_at_utc",
    "image_path",
    "image_sha256",
    "script_path",
    "command",
    "git_commit",
    "julia_version",
    "input_data_hashes",
]


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_image_provenance(image_path: Path, metadata: Dict) -> Dict:
    p = dict(metadata)
    p.setdefault("schema_version", "v1")
    p.setdefault("generated_at_utc", _now_utc_iso())
    p["image_path"] = str(image_path)
    p["image_sha256"] = compute_sha256(image_path)
    return p


def write_image_provenance_sidecar(image_path: Path, metadata: Dict) -> Path:
    image_path = Path(image_path)
    payload = build_image_provenance(image_path, metadata)
    sidecar = image_path.with_name(f"{image_path.name}.provenance.json")
    sidecar.write_text(json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return sidecar


def verify_image_provenance(image_path: Path) -> Tuple[bool, str]:
    image_path = Path(image_path)
    sidecar = image_path.with_name(f"{image_path.name}.provenance.json")
    if not sidecar.exists():
        return False, "sidecar not found"
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except Exception as e:
        return False, f"invalid sidecar json: {e}"

    for key in REQUIRED_V1_FIELDS:
        if key not in payload:
            return False, f"missing required field: {key}"

    actual = compute_sha256(image_path)
    expected = str(payload.get("image_sha256", ""))
    if actual != expected:
        return False, "image sha256 mismatch"
    return True, "ok"
```

- [ ] **Step 2: Re-run tests to confirm GREEN**

Run: `python -m pytest tests/unit/python/test_image_provenance_sidecar.py -q`  
Expected: PASS。

- [ ] **Step 3: Commit module implementation**

```bash
git add scripts/common/provenance_image.py tests/unit/python/test_image_provenance_sidecar.py
git commit -m "refactor(scripts): add reusable image provenance sidecar module"
```

---

### Task 3: 扩展失败路径测试（缺字段/校验入口）

**Files:**
- Modify: `tests/unit/python/test_image_provenance_sidecar.py`

- [ ] **Step 1: Add failing tests for missing fields and verify pass path**

```python
def test_verify_passes_for_untampered_image(tmp_path: Path):
    image = tmp_path / "ok.png"
    image.write_bytes(b"stable")
    metadata = {
        "schema_version": "v1",
        "generated_at_utc": "2026-04-19T00:00:00Z",
        "script_path": "scripts/plot_scan_csv.py",
        "command": "python scripts/plot_scan_csv.py --mode lines ...",
        "git_commit": "deadbeef",
        "julia_version": "python-3.x",
        "input_data_hashes": [],
    }
    write_image_provenance_sidecar(image, metadata)
    ok, reason = verify_image_provenance(image)
    assert ok is True
    assert reason == "ok"


def test_verify_fails_on_missing_required_field(tmp_path: Path):
    image = tmp_path / "missing.png"
    image.write_bytes(b"stable")
    sidecar = image.with_name(f"{image.name}.provenance.json")
    sidecar.write_text('{"schema_version":"v1"}', encoding="utf-8")

    ok, reason = verify_image_provenance(image)
    assert ok is False
    assert "missing required field" in reason
```

- [ ] **Step 2: Run tests to verify RED then GREEN**

Run: `python -m pytest tests/unit/python/test_image_provenance_sidecar.py -q`  
Expected: 新增测试先 FAIL；修复后 PASS。

- [ ] **Step 3: Commit test hardening**

```bash
git add tests/unit/python/test_image_provenance_sidecar.py scripts/common/provenance_image.py
git commit -m "test(scripts): harden image provenance verification failure paths"
```

---

### Task 4: 接入 `plot_scan_csv.py` 默认写 sidecar

**Files:**
- Modify: `scripts/plot_scan_csv.py`

- [ ] **Step 1: Add helper to collect run context**

在脚本中新增 helper（最小实现）：

```python
import subprocess
import sys
from scripts.common.provenance_image import write_image_provenance_sidecar, verify_image_provenance, compute_sha256


def _git_commit_or_unknown() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(PROJECT_ROOT), text=True).strip()
    except Exception:
        return "unknown"


def _build_input_hashes(csv_paths: list[Path]) -> list[dict]:
    out = []
    for p in csv_paths:
        if p.exists():
            out.append({"path": str(p), "sha256": compute_sha256(p)})
    return out
```

- [ ] **Step 2: Hook sidecar writing after each savefig**

在所有 `fig.savefig(...)` 后调用：

```python
metadata = {
    "schema_version": "v1",
    "script_path": "scripts/plot_scan_csv.py",
    "command": " ".join(sys.argv),
    "git_commit": _git_commit_or_unknown(),
    "julia_version": f"python-{sys.version.split()[0]}",
    "input_data_hashes": _build_input_hashes([csv_path]),
}
write_image_provenance_sidecar(out, metadata)
```

- [ ] **Step 3: Commit integration**

```bash
git add scripts/plot_scan_csv.py
git commit -m "refactor(scripts): emit image provenance sidecars from plot_scan_csv outputs"
```

---

### Task 5: 增加 CLI 校验入口

**Files:**
- Modify: `scripts/plot_scan_csv.py`

- [ ] **Step 1: Add CLI argument and short-circuit verify mode**

```python
ap.add_argument("--verify-provenance", type=str, default=None, help="Verify image sidecar provenance and exit")
```

在 `main()` 前部处理：

```python
if args.verify_provenance:
    ok, reason = verify_image_provenance(Path(args.verify_provenance))
    print("PASS" if ok else "FAIL", reason)
    return 0 if ok else 2
```

- [ ] **Step 2: Add unit test for CLI verify exit behavior**

在 Python 测试中新增调用脚本子进程断言：

```python
proc = subprocess.run([sys.executable, "scripts/plot_scan_csv.py", "--verify-provenance", str(image)], capture_output=True, text=True)
assert proc.returncode == 0
assert "PASS" in proc.stdout
```

并增加篡改后 `returncode != 0` 断言。

- [ ] **Step 3: Commit verify CLI**

```bash
git add scripts/plot_scan_csv.py tests/unit/python/test_image_provenance_sidecar.py
git commit -m "feat(scripts): add provenance verification CLI for image sidecars"
```

---

### Task 6: 文档边界说明与使用示例

**Files:**
- Modify: `docs/guides/scripts/README.md`

- [ ] **Step 1: Add sidecar usage and boundary section**

补充文档段落：

```markdown
### Image Provenance Sidecar (v1)

- `plot_scan_csv.py` 默认在每张输出图旁生成 `<image>.provenance.json`。
- sidecar 提供受控环境可复现所需最小元数据与 `image_sha256` 绑定校验。
- 本方案不保证图片脱离 sidecar 后仍可验证来源，不覆盖开放传播溯源。

验证示例：

```bash
python scripts/plot_scan_csv.py --verify-provenance data/outputs/figures/example.png
```
```

- [ ] **Step 2: Commit docs update**

```bash
git add docs/guides/scripts/README.md
git commit -m "docs(guides): document image provenance sidecar scope and verification"
```

---

### Task 7: 运行 #82 验证命令与收口检查

**Files:**
- Modify (if needed): `docs/guides/scripts/README.md`

- [ ] **Step 1: Run Python unit tests**

Run: `python -m pytest tests/unit/python/test_image_provenance_sidecar.py -q`  
Expected: PASS。

- [ ] **Step 2: Run issue-required integration smoke**

Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`  
Expected: PASS。

- [ ] **Step 3: Run issue-required unit smoke**

Run: `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`  
Expected: PASS。

- [ ] **Step 4: Final commit (if verification evidence/doc tweaks needed)**

```bash
git add docs/guides/scripts/README.md
git commit -m "test(ci): confirm image provenance sidecar smoke verification for issue 82"
```

---

## Plan Self-Review

- Spec coverage: 已覆盖通用模块、`plot_scan_csv.py` 链路落地、校验入口、测试、文档边界。
- Placeholder scan: 无 `TBD/TODO/implement later` 占位内容。
- Type consistency: 全流程统一使用 `write_image_provenance_sidecar` 与 `verify_image_provenance` 命名，sidecar 必填字段与 issue v1 对齐。
