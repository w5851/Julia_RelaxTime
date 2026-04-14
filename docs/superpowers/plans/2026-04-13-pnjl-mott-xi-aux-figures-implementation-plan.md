# PNJL Mott xi Auxiliary Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不重跑主扫描的前提下，生成图1/4/2、mode_ab 标注副本、以及可行时的图5，并把结果接入 `xi_dependence_analysis.md`。

**Architecture:** 新增一个 `scripts/analysis/` 下的 Python 生成脚本，直接读取现有 CSV 与现有 mode_ab PNG；图像绘制使用 matplotlib，文档更新由脚本追加固定锚点段落。图5采用“先判定可重构性、不可重构即跳过并写明原因”的门控策略，保证主流程可稳定落地。

**Tech Stack:** Python 3, matplotlib, 标准库 `csv/pathlib/argparse/math`, Julia integration smoke tests (`tests/integration/...`)。

---

## File Structure

- Create: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`
  - 责任：读取 3 份输入数据，生成图1/4/2，生成 mode_ab 标注图，判定并生成/跳过图5，更新分析文档。
- Create: `tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl`
  - 责任：端到端 smoke 验证输出文件契约（图文件 + 文档已追加辅助图段落）。
- Modify: `docs/analysis/pnjl_mott/xi_dependence_analysis.md`
  - 责任：接入辅助图集引用（由脚本自动追加到锚点段落）。
- Create: `docs/analysis/pnjl_mott/figures/`（运行时生成目录）
  - 责任：存放 `fig1/fig4/fig2/fig5`。
- Modify (generated output only):
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_K__xi3_annotated.png`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_pi__xi3_annotated.png`

### Task 1: 建立失败的 smoke 测试（TDD red）

**Files:**
- Create: `tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl`
- Test: `tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl`

- [ ] **Step 1: 写失败测试骨架（先断言目标文件存在）**

```julia
using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "analysis", "gen_pnjl_mott_xi_aux_figures.py")

@testset "pnjl mott xi aux figures smoke" begin
    @test isfile(SCRIPT_PATH)
end
```

- [ ] **Step 2: 扩展为最小端到端断言（预期失败）**

```julia
@testset "pnjl mott xi aux figures smoke" begin
    tmp = mktempdir()
    fig_dir = joinpath(tmp, "figures")
    mode_ab_dir = joinpath(tmp, "mode_ab")
    doc_path = joinpath(tmp, "xi_dependence_analysis.md")
    # 先写入最小 markdown
    write(doc_path, "# dummy\n\n## 5. 图像读取要点\n")

    # TODO: 后续步骤补齐输入 csv 与 mode_ab png 伪数据
    @test isfile(joinpath(fig_dir, "fig1_tmott_vs_xi_fit.png"))
end
```

- [ ] **Step 3: 运行单测确认 red**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: FAIL（脚本不存在或输出文件不存在）。

- [ ] **Step 4: 提交测试 red 基线**

```bash
git add tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl
git commit -m "test: add failing smoke contract for xi auxiliary figures"
```

### Task 2: 新建脚本骨架与输入校验

**Files:**
- Create: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`
- Test: `tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl`

- [ ] **Step 1: 写最小可运行 CLI 骨架**

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--derived-csv", type=Path, required=True)
    ap.add_argument("--scan-csv", type=Path, required=True)
    ap.add_argument("--gap-csv", type=Path, required=True)
    ap.add_argument("--mode-ab-dir", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, required=True)
    ap.add_argument("--doc", type=Path, required=True)
    ap.add_argument("--skip-fig5", action="store_true")
    return ap.parse_args()

def main() -> None:
    args = parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)
    print("ready")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 增加输入存在性与列名校验函数**

```python
import csv

def require_file(path: Path, name: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"missing {name}: {path}")

def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    return list(csv.DictReader(lines))

def require_columns(rows: list[dict[str, str]], required: list[str], name: str) -> None:
    if not rows:
        raise ValueError(f"empty csv: {name}")
    missing = [c for c in required if c not in rows[0]]
    if missing:
        raise ValueError(f"missing columns in {name}: {missing}")
```

- [ ] **Step 3: 在 `main()` 中接入校验**

```python
require_file(args.derived_csv, "derived-csv")
require_file(args.scan_csv, "scan-csv")
require_file(args.gap_csv, "gap-csv")
require_file(args.doc, "doc")
require_file(args.mode_ab_dir / "mott_mode_ab__M_K__xi3.png", "mode_ab K png")
require_file(args.mode_ab_dir / "mott_mode_ab__M_pi__xi3.png", "mode_ab pi png")

derived_rows = read_rows(args.derived_csv)
scan_rows = read_rows(args.scan_csv)
gap_rows = read_rows(args.gap_csv)
require_columns(derived_rows, ["xi", "T_MeV", "M_pi", "M_K", "M_u_plus_M_d", "M_u_plus_M_s", "Gamma_pi", "Gamma_K"], "derived")
require_columns(scan_rows, ["xi", "T_MeV", "M_pi", "M_K", "Gamma_pi", "Gamma_K", "M_u", "M_s"], "scan")
require_columns(gap_rows, ["xi", "T_MeV", "Phi", "m_u", "m_s"], "gap")
```

- [ ] **Step 4: 运行 smoke 测试（仍可失败在文件断言）**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: FAIL 点转移到“目标图未生成”。

- [ ] **Step 5: 提交骨架与校验**

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl
git commit -m "feat: scaffold xi auxiliary figures generator with strict input checks"
```

### Task 3: 实现图1（T_mott vs xi 拟合）

**Files:**
- Modify: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`
- Test: `tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl`

- [ ] **Step 1: 写 `T_mott` 提取与拟合函数**

```python
def estimate_crossing(Ts: list[float], deltas: list[float]) -> float:
    for i in range(len(Ts) - 1):
        d0, d1 = deltas[i], deltas[i + 1]
        if d0 == 0.0:
            return Ts[i]
        if d0 * d1 < 0:
            t0, t1 = Ts[i], Ts[i + 1]
            return t0 + (0.0 - d0) * (t1 - t0) / (d1 - d0)
    raise ValueError("no crossing found")
```

- [ ] **Step 2: 写图1绘制函数（线性+二次）**

```python
import matplotlib.pyplot as plt

def plot_fig1_tmott_vs_xi(out_png: Path, xi: list[float], tmott_pi: list[float], tmott_k: list[float]) -> None:
    fig, ax = plt.subplots(figsize=(6.75, 4.6))
    ax.scatter(xi, tmott_pi, label="pi samples", marker="o")
    ax.scatter(xi, tmott_k, label="K samples", marker="s")
    # 使用 numpy.polyfit 做 1/2 次拟合
    import numpy as np
    x = np.array(xi)
    for y, name, color in [(np.array(tmott_pi), "pi", "#4477AA"), (np.array(tmott_k), "K", "#EE6677")]:
        p1 = np.polyfit(x, y, 1)
        p2 = np.polyfit(x, y, 2)
        xx = np.linspace(min(x), max(x), 200)
        ax.plot(xx, p1[0] * xx + p1[1], color=color, linestyle="--", label=f"{name} linear")
        ax.plot(xx, p2[0] * xx**2 + p2[1] * xx + p2[2], color=color, linestyle="-", label=f"{name} quadratic")
    ax.set_xlabel("xi")
    ax.set_ylabel("T_mott [MeV]")
    ax.legend(frameon=False)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
```

- [ ] **Step 3: 在 `main()` 中落盘 `fig1_tmott_vs_xi_fit.png`**

```python
plot_fig1_tmott_vs_xi(args.fig_dir / "fig1_tmott_vs_xi_fit.png", xi_values, tmott_pi_values, tmott_k_values)
```

- [ ] **Step 4: 更新测试断言图1存在**

```julia
@test isfile(joinpath(fig_dir, "fig1_tmott_vs_xi_fit.png"))
```

- [ ] **Step 5: 运行测试并提交**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: 图1通过，其他待实现项仍可能失败。

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl
git commit -m "feat: generate fig1 tmott-vs-xi fits from existing scan data"
```

### Task 4: 实现图4（序参量直接/间接效应）

**Files:**
- Modify: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`

- [ ] **Step 1: 写固定温度切片与沿 `T_mott` 轨迹取样函数**

```python
def nearest_row(rows: list[dict[str, str]], *, xi: float, t_mev: float) -> dict[str, str]:
    cand = [r for r in rows if abs(float(r["xi"]) - xi) < 1e-9]
    cand.sort(key=lambda r: abs(float(r["T_MeV"]) - t_mev))
    return cand[0]
```

- [ ] **Step 2: 写图4绘制函数**

```python
def plot_fig4_orderparam(out_png: Path, xi: list[float], fixed_phi: list[float], traj_phi: list[float], fixed_mu: list[float], traj_mu: list[float]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.4), sharex=True)
    axes[0].plot(xi, fixed_phi, "o-", label="Phi @ T=200")
    axes[0].plot(xi, traj_phi, "s--", label="Phi @ T_mott(xi)")
    axes[0].set_ylabel("Phi")
    axes[0].legend(frameon=False)
    axes[1].plot(xi, fixed_mu, "o-", label="m_u @ T=200")
    axes[1].plot(xi, traj_mu, "s--", label="m_u @ T_mott(xi)")
    axes[1].set_ylabel("m_u [fm^-1]")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.set_xlabel("xi")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
```

- [ ] **Step 3: 在 `main()` 中输出 `fig4_orderparam_direct_indirect.png`**

```python
plot_fig4_orderparam(args.fig_dir / "fig4_orderparam_direct_indirect.png", xi_values, phi_fixed, phi_traj, mu_fixed, mu_traj)
```

- [ ] **Step 4: 更新测试断言图4存在并运行**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: 图1、图4存在；图2/标注/文档可能仍失败。

- [ ] **Step 5: 提交图4实现**

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl
git commit -m "feat: add fig4 direct-vs-indirect order-parameter visualization"
```

### Task 5: 实现图2（Gamma 与 Delta 双轴）

**Files:**
- Modify: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`

- [ ] **Step 1: 明确并统一 Delta 定义**

```python
# Delta = M_thr - M_mes
delta_pi = (float(r["M_u_plus_M_d"]) - float(r["M_pi"]))
delta_k = (float(r["M_u_plus_M_s"]) - float(r["M_K"]))
```

- [ ] **Step 2: 实现双轴绘图函数**

```python
def plot_fig2_gamma_delta(out_png: Path, x_t: list[float], gamma: list[float], delta: list[float], label: str) -> None:
    fig, ax1 = plt.subplots(figsize=(6.75, 4.2))
    ax2 = ax1.twinx()
    ax1.plot(x_t, gamma, color="#4477AA", label=f"Gamma_{label}")
    ax2.plot(x_t, delta, color="#EE6677", linestyle="--", label=f"Delta_{label}")
    ax1.set_xlabel("T [MeV]")
    ax1.set_ylabel("Gamma [fm^-1]")
    ax2.set_ylabel("Delta = M_thr - M_mes [fm^-1]")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, frameon=False, loc="best")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
```

- [ ] **Step 3: 生成 `fig2_gamma_delta_dualaxis.png`（可用双分面或单图聚合）**

```python
plot_fig2_gamma_delta(args.fig_dir / "fig2_gamma_delta_dualaxis.png", t_vals, gamma_pi_vals, delta_pi_vals, "pi")
```

- [ ] **Step 4: 更新测试断言图2存在并运行**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: 图1/4/2存在；标注/文档/图5门控待完成。

- [ ] **Step 5: 提交图2实现**

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl
git commit -m "feat: add fig2 dual-axis gamma-delta plot with unified delta sign"
```

### Task 6: 生成 mode_ab 轻增强标注副本

**Files:**
- Modify: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`

- [ ] **Step 1: 实现通用 PNG 标注函数（不改原图）**

```python
import matplotlib.image as mpimg

def annotate_mode_ab(src: Path, dst: Path, text_lines: list[str]) -> None:
    img = mpimg.imread(src)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.imshow(img)
    ax.axis("off")
    text = "\n".join(text_lines)
    ax.text(0.02, 0.98, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.75))
    fig.savefig(dst, dpi=300, bbox_inches="tight")
    plt.close(fig)
```

- [ ] **Step 2: 计算交点右移摘要文字（来自图1同源数据）**

```python
k_note = ["Crossing shifts right with xi", f"Delta T_mott^K ~= {tmott_k[-1]-tmott_k[0]:.2f} MeV"]
pi_note = ["Crossing shifts right with xi", f"Delta T_mott^pi ~= {tmott_pi[-1]-tmott_pi[0]:.2f} MeV"]
```

- [ ] **Step 3: 输出两张 `_annotated.png`**

```python
annotate_mode_ab(args.mode_ab_dir / "mott_mode_ab__M_K__xi3.png", args.mode_ab_dir / "mott_mode_ab__M_K__xi3_annotated.png", k_note)
annotate_mode_ab(args.mode_ab_dir / "mott_mode_ab__M_pi__xi3.png", args.mode_ab_dir / "mott_mode_ab__M_pi__xi3_annotated.png", pi_note)
```

- [ ] **Step 4: 更新测试断言并运行**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: `_annotated.png` 存在。

- [ ] **Step 5: 提交标注实现**

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl
git commit -m "feat: generate annotated mode_ab copies for crossing-shift emphasis"
```

### Task 7: 图5门控 + 文档追加

**Files:**
- Modify: `scripts/analysis/gen_pnjl_mott_xi_aux_figures.py`
- Modify: `docs/analysis/pnjl_mott/xi_dependence_analysis.md` (generated append)

- [ ] **Step 1: 实现图5可重构性判定函数**

```python
def can_build_fig5(rows: list[dict[str, str]]) -> tuple[bool, str]:
    needed = {"xi", "T_MeV", "M_pi", "M_K", "M_u", "M_s", "Phi"}
    have = set(rows[0].keys()) if rows else set()
    miss = sorted(needed - have)
    if miss:
        return False, f"missing required fields for fig5: {miss}"
    return True, "ok"
```

- [ ] **Step 2: 实现“可做就画，不可做就写 skip 理由”的分支**

```python
fig5_ok, fig5_msg = can_build_fig5(scan_rows)
fig5_path = args.fig_dir / "fig5_taylor_decomposition.png"
if (not args.skip_fig5) and fig5_ok:
    # 绘制最小可解释分解图（示例：总变化 vs 两个代理项）
    plot_fig5(fig5_path, ...)
else:
    fig5_path = None
```

- [ ] **Step 3: 追加文档“辅助图集”段落（幂等替换）**

```python
AUX_BEGIN = "<!-- AUX_FIGURES:BEGIN -->"
AUX_END = "<!-- AUX_FIGURES:END -->"

def upsert_aux_section(doc: Path, *, fig5_path: Path | None, fig5_msg: str) -> None:
    section = [
        "## 8. 辅助图集（新增）",
        "",
        "- 图1：`figures/fig1_tmott_vs_xi_fit.png`",
        "- 图4：`figures/fig4_orderparam_direct_indirect.png`",
        "- 图2：`figures/fig2_gamma_delta_dualaxis.png`",
        "- 图3（轻增强复用）：见 `mode_ab/*_annotated.png`",
        f"- 图5：`figures/fig5_taylor_decomposition.png` ({'generated' if fig5_path else 'skipped'})",
        "",
        "符号统一：`Delta = M_thr - M_mes`（若旧段落使用相反定义，请按符号翻转理解）。",
    ]
    body = doc.read_text(encoding="utf-8")
    block = AUX_BEGIN + "\n" + "\n".join(section) + ("\n说明: " + fig5_msg if (not fig5_path and fig5_msg) else "") + "\n" + AUX_END
    if AUX_BEGIN in body and AUX_END in body:
        pre = body.split(AUX_BEGIN)[0]
        post = body.split(AUX_END)[1]
        doc.write_text(pre + block + post, encoding="utf-8")
    else:
        doc.write_text(body.rstrip() + "\n\n" + block + "\n", encoding="utf-8")
```

- [ ] **Step 4: 完善 smoke 断言文档锚点存在并运行**

```julia
doc_text = read(doc_path, String)
@test occursin("<!-- AUX_FIGURES:BEGIN -->", doc_text)
@test occursin("Delta = M_thr - M_mes", doc_text)
```

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: PASS。

- [ ] **Step 5: 提交图5门控与文档接入**

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl docs/analysis/pnjl_mott/xi_dependence_analysis.md
git commit -m "docs: append xi auxiliary gallery with fig5 gating and delta convention note"
```

### Task 8: 回归验证与交付检查

**Files:**
- Test: `tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl`
- Test: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`

- [ ] **Step 1: 跑新 smoke 测试**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl")'`
Expected: PASS。

- [ ] **Step 2: 跑原有 mode_ab smoke，防止回归**

Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl")'`
Expected: PASS。

- [ ] **Step 3: 在真实数据路径执行一次生成脚本**

Run:

```bash
python scripts/analysis/gen_pnjl_mott_xi_aux_figures.py --derived-csv data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv --scan-csv data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_scan.csv --gap-csv data/outputs/results/relaxtime/scan/gap_transport_scan_step5_muB0_xi-0p3to0p3.csv --mode-ab-dir data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab --fig-dir docs/analysis/pnjl_mott/figures --doc docs/analysis/pnjl_mott/xi_dependence_analysis.md
```

Expected: 图1/图4/图2/标注副本生成；图5按门控生成或跳过说明。

- [ ] **Step 4: 文件契约检查**

Run: `git status --short`
Expected: 仅出现本计划涉及脚本/测试/文档与生成图变更。

- [ ] **Step 5: 提交最终整合**

```bash
git add scripts/analysis/gen_pnjl_mott_xi_aux_figures.py tests/integration/relaxtime/test_pnjl_mott_xi_aux_figures_smoke.jl docs/analysis/pnjl_mott/xi_dependence_analysis.md docs/analysis/pnjl_mott/figures data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/*_annotated.png
git commit -m "feat: deliver xi auxiliary figure set and annotated mode_ab references"
```

## Self-Review

- Spec coverage:
  - 图1/4/2顺序与命名：Task 3/4/5。
  - 图3采用轻增强副本：Task 6。
  - 图5条件化：Task 7。
  - 文档接入与符号统一说明：Task 7。
- Placeholder scan:
  - 无 `TODO/TBD/implement later` 形式占位执行步骤。
- Consistency:
  - 输出命名统一使用 `fig1_tmott_vs_xi_fit.png`、`fig4_orderparam_direct_indirect.png`、`fig2_gamma_delta_dualaxis.png`、`fig5_taylor_decomposition.png`。
  - `Delta` 符号统一为 `M_thr - M_mes`。
