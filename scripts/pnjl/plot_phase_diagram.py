#!/usr/bin/env python3
"""
PNJL 相图绘制脚本

从 data/reference/pnjl/ 读取相结构数据，绘制：
1. T-μ 相图（一阶相变线 + CEP）
2. T-ρ 相图（共存区）

用法：
    python scripts/pnjl/plot_phase_diagram.py [options]

选项：
    --boundary PATH   相变线数据文件（显式旧 schema 输入）
    --cep PATH        CEP 数据文件（显式旧 schema 输入）
    --xi VALUE        要绘制的 ξ 值，可多次指定 (默认绘制所有)
    --output-dir DIR  输出目录 (默认 data/outputs/figures/pnjl)
    --format FMT      输出格式 png/pdf/svg (默认 png)
    --dpi VALUE       分辨率 (默认 200)
    --no-show         不显示图形窗口

审阅交付模式（Issue 130 C2，PNG + manifest）：
    python scripts/pnjl/plot_phase_diagram.py --review-c2 \
        [--output-dir DIR --results-dir DIR] [--phase-reference-root PATH]

    --review-c2 只接受 accepted 层的原生认证行，产出 171.45 mm 审阅 PNG、
    灰度派生图与哈希绑定 manifest；--output-dir/--results-dir 必须成对给出，
    省略时使用 transport_issue130_c2_v1 的默认目录。该模式与 --format、
    --xi、--dpi、--type、--phase-reference-layer 等快速绘图参数互斥。
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


def _find_project_root() -> Path:
    """Find project root by looking for Project.toml or .git."""
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent, script_dir.parent.parent]
    candidates.append(Path.cwd())
    
    for start in candidates:
        current = start
        for _ in range(5):
            if (current / "Project.toml").exists() or (current / ".git").exists():
                return current
            parent = current.parent
            if parent == current:
                break
            current = parent
    return Path.cwd()


PROJECT_ROOT = _find_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.pnjl.phase_reference_adapter import (
    default_downstream_layer, estimate_c2_endpoint_density, load_c2_review_sources, load_phase_reference,
)


DEFAULT_PHASE_REFERENCE_ROOT = PROJECT_ROOT / "data" / "reference" / "pnjl" / "issue130_phase_reference_v2"
DEFAULT_PHASE_REFERENCE_LAYER = "accepted"
DEFAULT_ACCEPTED_TABLE_ROOT = DEFAULT_PHASE_REFERENCE_ROOT / "accepted" / "tables"
DEFAULT_BOUNDARY_PATH = DEFAULT_ACCEPTED_TABLE_ROOT / "maxwell_surface_accepted_phase_map_v1.csv"
DEFAULT_CEP_PATH = DEFAULT_ACCEPTED_TABLE_ROOT / "cep_boundary_accepted_phase_map_v1.csv"
DEFAULT_SPINODAL_PATH = DEFAULT_ACCEPTED_TABLE_ROOT / "spinodal_surface_accepted_phase_map_v1.csv"
DEFAULT_CROSSOVER_PATH = DEFAULT_ACCEPTED_TABLE_ROOT / "crossover_surface_accepted_phase_map_v1.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "outputs" / "figures" / "pnjl"


@dataclass
class BoundaryPoint:
    """一阶相变线上的一个点"""
    xi: float
    T_MeV: float
    mu_coex_MeV: float
    rho_hadron: float  # 强子相密度（低密度侧）
    rho_quark: float   # 夸克相密度（高密度侧）


@dataclass
class CEPPoint:
    """临界终点"""
    xi: float
    T_CEP_MeV: float
    mu_CEP_MeV: float


@dataclass
class SpinodalPoint:
    """Spinodal（亚稳态边界）点"""
    xi: float
    T_MeV: float
    mu_spinodal_low_MeV: float
    mu_spinodal_high_MeV: float
    rho_spinodal_hadron: float  # 强子相侧 spinodal
    rho_spinodal_quark: float   # 夸克相侧 spinodal


@dataclass
class CrossoverPoint:
    """Crossover 点"""
    xi: float
    mu_MeV: float
    T_chiral_MeV: Optional[float]  # 手征 crossover 温度
    T_deconf_MeV: Optional[float]  # 退禁闭 crossover 温度
    rho_chiral: Optional[float]    # 手征 crossover 密度
    rho_deconf: Optional[float]    # 退禁闭 crossover 密度


def load_candidate_phase_data(reference_root: Path, layer: str | None = None) -> tuple[
    List[BoundaryPoint], List[CEPPoint], List[SpinodalPoint], List[CrossoverPoint], dict
]:
    """Load an explicit Issue #130 layer for diagnostic plotting.

    This path keeps the candidate's unresolved/interpolated status in the
    returned diagnostics.  The v2 ``accepted`` layer is the default for
    downstream phase-map consumers; ``strict`` remains an explicit audit view.
    """

    selected_layer = layer or default_downstream_layer(reference_root)
    bundle = load_phase_reference(reference_root, layer=selected_layer)
    boundary_points = [
        BoundaryPoint(
            xi=row["xi"],
            T_MeV=row["T_MeV"],
            mu_coex_MeV=row["muq_MeV"],
            rho_hadron=row["rho_hadron"],
            rho_quark=row["rho_quark"],
        )
        for row in bundle.tables.get("boundary", ())
    ]
    cep_points = [
        CEPPoint(
            xi=row["xi"],
            T_CEP_MeV=row["T_midpoint_MeV"],
            mu_CEP_MeV=row["muq_CEP_proxy_MeV"],
        )
        for row in bundle.tables.get("cep", ())
    ]
    spinodal_points = [
        SpinodalPoint(
            xi=row["xi"],
            T_MeV=row["T_MeV"],
            mu_spinodal_low_MeV=row["muq_spinodal_hadron_MeV"],
            mu_spinodal_high_MeV=row["muq_spinodal_quark_MeV"],
            rho_spinodal_hadron=row["rho_spinodal_hadron"] or 0.0,
            rho_spinodal_quark=row["rho_spinodal_quark"] or 0.0,
        )
        for row in bundle.tables.get("spinodals", ())
    ]
    # The candidate contains one physical crossover family.  Map it to the
    # historical chiral slot and leave the separate deconfinement slot empty;
    # no second physical surface is invented by the adapter.
    crossover_points = [
        CrossoverPoint(
            xi=row["xi"],
            mu_MeV=row["muq_MeV"],
            T_chiral_MeV=row["T_MeV"],
            T_deconf_MeV=None,
            rho_chiral=row["rho"],
            rho_deconf=None,
        )
        for row in bundle.tables.get("crossover", ())
        if row.get("physical_region", "") in ("", "crossover_below_CEP")
    ]
    return boundary_points, cep_points, spinodal_points, crossover_points, dict(bundle.diagnostics)


def load_boundary_data(path: Path) -> List[BoundaryPoint]:
    """加载相变线数据"""
    if not path.exists():
        raise FileNotFoundError(f"Boundary file not found: {path}")
    
    points = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # 支持新旧命名
                rho_hadron = float(row.get("rho_hadron") or row.get("rho_gas", 0))
                rho_quark = float(row.get("rho_quark") or row.get("rho_liquid", 0))
                mu = float(row.get("mu_transition_MeV") or row.get("mu_MeV") or row.get("mu_coex_MeV", 0))
                points.append(BoundaryPoint(
                    xi=float(row["xi"]),
                    T_MeV=float(row["T_MeV"]),
                    mu_coex_MeV=mu,
                    rho_hadron=rho_hadron,
                    rho_quark=rho_quark,
                ))
            except (KeyError, ValueError):
                continue
    return points


def load_cep_data(path: Path) -> List[CEPPoint]:
    """加载 CEP 数据"""
    if not path.exists():
        return []
    
    points = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                points.append(CEPPoint(
                    xi=float(row["xi"]),
                    T_CEP_MeV=float(row["T_CEP_MeV"]),
                    mu_CEP_MeV=float(row["mu_CEP_MeV"]),
                ))
            except (KeyError, ValueError):
                continue
    return points


def load_spinodal_data(path: Path) -> List[SpinodalPoint]:
    """加载 spinodal 数据"""
    if not path.exists():
        return []
    
    points = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # 支持新旧命名
                rho_hadron = float(row.get("rho_spinodal_hadron") or row.get("rho_spinodal_low", 0))
                rho_quark = float(row.get("rho_spinodal_quark") or row.get("rho_spinodal_high", 0))
                mu_hadron = float(row.get("mu_spinodal_hadron_MeV") or row.get("mu_spinodal_low_MeV", 0))
                mu_quark = float(row.get("mu_spinodal_quark_MeV") or row.get("mu_spinodal_high_MeV", 0))
                points.append(SpinodalPoint(
                    xi=float(row["xi"]),
                    T_MeV=float(row["T_MeV"]),
                    mu_spinodal_low_MeV=mu_hadron,
                    mu_spinodal_high_MeV=mu_quark,
                    rho_spinodal_hadron=rho_hadron,
                    rho_spinodal_quark=rho_quark,
                ))
            except (KeyError, ValueError):
                continue
    return points


def load_crossover_data(path: Path) -> List[CrossoverPoint]:
    """加载 crossover 数据"""
    if not path.exists():
        return []
    
    points = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                T_chiral = row.get("T_crossover_chiral_MeV", "")
                T_deconf = row.get("T_crossover_deconf_MeV", "")
                rho_chiral = row.get("rho_chiral", "")
                rho_deconf = row.get("rho_deconf", "")
                points.append(CrossoverPoint(
                    xi=float(row["xi"]),
                    mu_MeV=float(row["mu_MeV"]),
                    T_chiral_MeV=float(T_chiral) if T_chiral else None,
                    T_deconf_MeV=float(T_deconf) if T_deconf else None,
                    rho_chiral=float(rho_chiral) if rho_chiral else None,
                    rho_deconf=float(rho_deconf) if rho_deconf else None,
                ))
            except (KeyError, ValueError):
                continue
    return points


def group_by_xi(points: List[BoundaryPoint]) -> Dict[float, List[BoundaryPoint]]:
    """按 ξ 值分组"""
    groups: Dict[float, List[BoundaryPoint]] = {}
    for p in points:
        groups.setdefault(p.xi, []).append(p)
    # 按温度排序
    for xi in groups:
        groups[xi].sort(key=lambda p: p.T_MeV)
    return groups


def get_cep_for_xi(ceps: List[CEPPoint], xi: float, tol: float = 1e-6) -> Optional[CEPPoint]:
    """获取指定 ξ 的 CEP"""
    for cep in ceps:
        if abs(cep.xi - xi) < tol:
            return cep
    return None


def group_spinodals_by_xi(points: List[SpinodalPoint]) -> Dict[float, List[SpinodalPoint]]:
    """按 ξ 值分组 spinodal 数据"""
    groups: Dict[float, List[SpinodalPoint]] = {}
    for p in points:
        groups.setdefault(p.xi, []).append(p)
    # 按温度排序
    for xi in groups:
        groups[xi].sort(key=lambda p: p.T_MeV)
    return groups


def get_spinodals_for_xi(spinodals: List[SpinodalPoint], xi: float, tol: float = 1e-6) -> List[SpinodalPoint]:
    """获取指定 ξ 的 spinodal 数据"""
    return [s for s in spinodals if abs(s.xi - xi) < tol]


def group_crossover_by_xi(points: List[CrossoverPoint]) -> Dict[float, List[CrossoverPoint]]:
    """按 ξ 值分组 crossover 数据"""
    groups: Dict[float, List[CrossoverPoint]] = {}
    for p in points:
        groups.setdefault(p.xi, []).append(p)
    # 按 μ 排序
    for xi in groups:
        groups[xi].sort(key=lambda p: p.mu_MeV)
    return groups


def get_crossover_for_xi(crossovers: List[CrossoverPoint], xi: float, tol: float = 1e-6) -> List[CrossoverPoint]:
    """获取指定 ξ 的 crossover 数据"""
    return sorted([c for c in crossovers if abs(c.xi - xi) < tol], key=lambda c: c.mu_MeV)


# 颜色方案
XI_COLORS = {
    0.0: "#E41A1C",   # 红色
    0.2: "#377EB8",   # 蓝色
    0.4: "#4DAF4A",   # 绿色
    0.6: "#984EA3",   # 紫色
    0.8: "#FF7F00",   # 橙色
    1.0: "#A65628",   # 棕色
}


def get_color_for_xi(xi: float) -> str:
    """获取 ξ 对应的颜色"""
    if xi in XI_COLORS:
        return XI_COLORS[xi]
    # 默认颜色循环
    colors = list(XI_COLORS.values())
    idx = int(xi * 10) % len(colors)
    return colors[idx]


def plot_T_mu_phase_diagram(
    boundary_groups: Dict[float, List[BoundaryPoint]],
    ceps: List[CEPPoint],
    spinodals: Optional[List[SpinodalPoint]] = None,
    crossovers: Optional[List[CrossoverPoint]] = None,
    xi_filter: Optional[List[float]] = None,
    output_path: Optional[Path] = None,
    show: bool = True,
    dpi: int = 200,
) -> None:
    """绘制 T-μ 相图"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    xi_values = sorted(boundary_groups.keys())
    if xi_filter:
        xi_values = [xi for xi in xi_values if xi in xi_filter]
    
    spinodal_groups = group_spinodals_by_xi(spinodals) if spinodals else {}
    crossover_groups = group_crossover_by_xi(crossovers) if crossovers else {}
    
    for xi in xi_values:
        points = boundary_groups.get(xi, [])
        color = get_color_for_xi(xi)
        
        # 绘制 spinodal 线（亚稳态边界）
        xi_spinodals = spinodal_groups.get(xi, [])
        if xi_spinodals:
            T_vals = [s.T_MeV for s in xi_spinodals]
            mu_low = [s.mu_spinodal_low_MeV for s in xi_spinodals]
            mu_high = [s.mu_spinodal_high_MeV for s in xi_spinodals]
            ax.plot(mu_low, T_vals, '--', color=color, linewidth=1, alpha=0.6)
            ax.plot(mu_high, T_vals, '--', color=color, linewidth=1, alpha=0.6)
        
        # 绘制一阶相变线
        if points:
            mu_vals = [p.mu_coex_MeV for p in points]
            T_vals = [p.T_MeV for p in points]
            ax.plot(mu_vals, T_vals, '-', color=color, linewidth=2, 
                    label=f"ξ = {xi:.1f}" if len(xi_values) > 1 else "First-order transition")
        
        # 绘制 CEP
        cep = get_cep_for_xi(ceps, xi)
        if cep:
            ax.scatter([cep.mu_CEP_MeV], [cep.T_CEP_MeV], 
                      s=150, c=color, marker='*', edgecolors='black', 
                      linewidths=0.5, zorder=10,
                      label=f"CEP (ξ={xi:.1f})" if len(xi_values) > 1 else "CEP")
        
        # 绘制 crossover 线
        xi_crossovers = crossover_groups.get(xi, [])
        if xi_crossovers:
            # 手征 crossover（虚线）
            chiral_mu = [c.mu_MeV for c in xi_crossovers if c.T_chiral_MeV is not None]
            chiral_T = [c.T_chiral_MeV for c in xi_crossovers if c.T_chiral_MeV is not None]
            if chiral_mu:
                ax.plot(chiral_mu, chiral_T, ':', color=color, linewidth=2, alpha=0.8,
                       label="Chiral crossover" if xi == xi_values[0] else None)
            
            # 退禁闭 crossover（点划线）
            deconf_mu = [c.mu_MeV for c in xi_crossovers if c.T_deconf_MeV is not None]
            deconf_T = [c.T_deconf_MeV for c in xi_crossovers if c.T_deconf_MeV is not None]
            if deconf_mu:
                ax.plot(deconf_mu, deconf_T, '-.', color=color, linewidth=1.5, alpha=0.7,
                       label="Deconf. crossover" if xi == xi_values[0] else None)
    
    ax.set_xlabel(r"$\mu$ (MeV)", fontsize=12)
    ax.set_ylabel(r"$T$ (MeV)", fontsize=12)
    ax.set_title("PNJL Phase Diagram (T-μ)", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    
    # 设置合理的坐标范围
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_path}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_T_rho_phase_diagram(
    boundary_groups: Dict[float, List[BoundaryPoint]],
    ceps: List[CEPPoint],
    spinodals: Optional[List[SpinodalPoint]] = None,
    crossovers: Optional[List[CrossoverPoint]] = None,
    xi_filter: Optional[List[float]] = None,
    output_path: Optional[Path] = None,
    show: bool = True,
    dpi: int = 200,
) -> None:
    """绘制 T-ρ 相图（共存区和 crossover 线）"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    xi_values = sorted(boundary_groups.keys())
    if xi_filter:
        xi_values = [xi for xi in xi_values if xi in xi_filter]
    
    spinodal_groups = group_spinodals_by_xi(spinodals) if spinodals else {}
    crossover_groups = group_crossover_by_xi(crossovers) if crossovers else {}
    
    for xi in xi_values:
        points = boundary_groups.get(xi, [])
        if not points:
            continue
        
        T_vals = np.array([p.T_MeV for p in points])
        rho_hadron = np.array([p.rho_hadron for p in points])
        rho_quark = np.array([p.rho_quark for p in points])
        color = get_color_for_xi(xi)
        
        # 绘制 spinodal 线（亚稳态边界）
        xi_spinodals = spinodal_groups.get(xi, [])
        if xi_spinodals:
            T_sp = [s.T_MeV for s in xi_spinodals]
            rho_sp_hadron = [s.rho_spinodal_hadron for s in xi_spinodals]
            rho_sp_quark = [s.rho_spinodal_quark for s in xi_spinodals]
            ax.plot(rho_sp_hadron, T_sp, '--', color=color, linewidth=1, alpha=0.6)
            ax.plot(rho_sp_quark, T_sp, '--', color=color, linewidth=1, alpha=0.6)
        
        # 绘制共存区边界
        ax.plot(rho_hadron, T_vals, '-', color=color, linewidth=2)
        ax.plot(rho_quark, T_vals, '-', color=color, linewidth=2)
        
        # 填充共存区
        ax.fill_betweenx(T_vals, rho_hadron, rho_quark, 
                         color=color, alpha=0.2,
                         label=f"ξ = {xi:.1f}" if len(xi_values) > 1 else "Coexistence region")
        
        # 绘制 crossover 线（按 T 排序，因为 ρ 不是单调的）
        xi_crossovers = crossover_groups.get(xi, [])
        if xi_crossovers:
            # 手征 crossover (T-ρ)，按 T 降序排列（高温到低温）
            chiral_data = [(c.rho_chiral, c.T_chiral_MeV) for c in xi_crossovers 
                          if c.T_chiral_MeV is not None and c.rho_chiral is not None]
            if chiral_data:
                chiral_sorted = sorted(chiral_data, key=lambda x: x[1], reverse=True)
                rho_chiral, T_chiral = zip(*chiral_sorted)
                ax.plot(rho_chiral, T_chiral, '-.', color=color, linewidth=2, 
                       label=f"Chiral crossover (ξ={xi:.1f})" if len(xi_values) > 1 else "Chiral crossover")
            
            # 退禁闭 crossover (T-ρ)，按 T 降序排列
            deconf_data = [(c.rho_deconf, c.T_deconf_MeV) for c in xi_crossovers 
                          if c.T_deconf_MeV is not None and c.rho_deconf is not None]
            if deconf_data:
                deconf_sorted = sorted(deconf_data, key=lambda x: x[1], reverse=True)
                rho_deconf, T_deconf = zip(*deconf_sorted)
                ax.plot(rho_deconf, T_deconf, ':', color=color, linewidth=2,
                       label=f"Deconf crossover (ξ={xi:.1f})" if len(xi_values) > 1 else "Deconf crossover")
        
        # 绘制 CEP（在 ρ 图上用 (ρ_hadron + ρ_quark)/2 作为 x 坐标）
        cep = get_cep_for_xi(ceps, xi)
        if cep:
            # 找到最接近 CEP 温度的点
            idx = np.argmin(np.abs(T_vals - cep.T_CEP_MeV))
            rho_cep = (rho_hadron[idx] + rho_quark[idx]) / 2
            ax.scatter([rho_cep], [cep.T_CEP_MeV], 
                      s=150, c=color, marker='*', edgecolors='black',
                      linewidths=0.5, zorder=10)
    
    ax.set_xlabel(r"$\rho / \rho_0$", fontsize=12)
    ax.set_ylabel(r"$T$ (MeV)", fontsize=12)
    ax.set_title("PNJL Phase Diagram (T-ρ)", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_path}")
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_combined_phase_diagram(
    boundary_groups: Dict[float, List[BoundaryPoint]],
    ceps: List[CEPPoint],
    spinodals: Optional[List[SpinodalPoint]] = None,
    crossovers: Optional[List[CrossoverPoint]] = None,
    xi_filter: Optional[List[float]] = None,
    output_path: Optional[Path] = None,
    show: bool = True,
    dpi: int = 200,
) -> None:
    """绘制组合相图（T-μ 和 T-ρ 并排）"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    xi_values = sorted(boundary_groups.keys())
    if xi_filter:
        xi_values = [xi for xi in xi_values if xi in xi_filter]
    
    spinodal_groups = group_spinodals_by_xi(spinodals) if spinodals else {}
    crossover_groups = group_crossover_by_xi(crossovers) if crossovers else {}
    
    for xi in xi_values:
        points = boundary_groups.get(xi, [])
        if not points:
            continue
        
        mu_vals = [p.mu_coex_MeV for p in points]
        T_vals = np.array([p.T_MeV for p in points])
        rho_hadron = np.array([p.rho_hadron for p in points])
        rho_quark = np.array([p.rho_quark for p in points])
        color = get_color_for_xi(xi)
        label = f"ξ = {xi:.1f}" if len(xi_values) > 1 else "First-order"
        
        # T-μ 图: spinodal 线
        xi_spinodals = spinodal_groups.get(xi, [])
        if xi_spinodals:
            T_sp = [s.T_MeV for s in xi_spinodals]
            mu_low = [s.mu_spinodal_low_MeV for s in xi_spinodals]
            mu_high = [s.mu_spinodal_high_MeV for s in xi_spinodals]
            ax1.plot(mu_low, T_sp, '--', color=color, linewidth=1, alpha=0.6)
            ax1.plot(mu_high, T_sp, '--', color=color, linewidth=1, alpha=0.6)
            # T-ρ 图: spinodal 线
            rho_sp_hadron = [s.rho_spinodal_hadron for s in xi_spinodals]
            rho_sp_quark = [s.rho_spinodal_quark for s in xi_spinodals]
            ax2.plot(rho_sp_hadron, T_sp, '--', color=color, linewidth=1, alpha=0.6)
            ax2.plot(rho_sp_quark, T_sp, '--', color=color, linewidth=1, alpha=0.6)
        
        # T-μ 图: 一阶相变线
        ax1.plot(mu_vals, T_vals, '-', color=color, linewidth=2, label=label)
        
        # T-ρ 图: 共存区
        ax2.plot(rho_hadron, T_vals, '-', color=color, linewidth=2)
        ax2.plot(rho_quark, T_vals, '-', color=color, linewidth=2)
        ax2.fill_betweenx(T_vals, rho_hadron, rho_quark, color=color, alpha=0.2, label=label)
        
        # CEP
        cep = get_cep_for_xi(ceps, xi)
        if cep:
            ax1.scatter([cep.mu_CEP_MeV], [cep.T_CEP_MeV], 
                       s=150, c=color, marker='*', edgecolors='black', 
                       linewidths=0.5, zorder=10)
            idx = np.argmin(np.abs(T_vals - cep.T_CEP_MeV))
            rho_cep = (rho_hadron[idx] + rho_quark[idx]) / 2
            ax2.scatter([rho_cep], [cep.T_CEP_MeV], 
                       s=150, c=color, marker='*', edgecolors='black',
                       linewidths=0.5, zorder=10)
        
        # T-μ 图: crossover 线
        xi_crossovers = crossover_groups.get(xi, [])
        if xi_crossovers:
            # 手征 crossover
            chiral_mu = [c.mu_MeV for c in xi_crossovers if c.T_chiral_MeV is not None]
            chiral_T = [c.T_chiral_MeV for c in xi_crossovers if c.T_chiral_MeV is not None]
            if chiral_mu:
                ax1.plot(chiral_mu, chiral_T, ':', color=color, linewidth=2, alpha=0.8,
                        label="Chiral" if xi == xi_values[0] else None)
            
            # 退禁闭 crossover
            deconf_mu = [c.mu_MeV for c in xi_crossovers if c.T_deconf_MeV is not None]
            deconf_T = [c.T_deconf_MeV for c in xi_crossovers if c.T_deconf_MeV is not None]
            if deconf_mu:
                ax1.plot(deconf_mu, deconf_T, '-.', color=color, linewidth=1.5, alpha=0.7,
                        label="Deconf." if xi == xi_values[0] else None)
    
    # 设置 T-μ 图
    ax1.set_xlabel(r"$\mu$ (MeV)", fontsize=12)
    ax1.set_ylabel(r"$T$ (MeV)", fontsize=12)
    ax1.set_title("T-μ Phase Diagram", fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right")
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    
    # 设置 T-ρ 图
    ax2.set_xlabel(r"$\rho / \rho_0$", fontsize=12)
    ax2.set_ylabel(r"$T$ (MeV)", fontsize=12)
    ax2.set_title("T-ρ Phase Diagram", fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right")
    ax2.set_xlim(left=0)
    ax2.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_path}")
    
    if show:
        plt.show()
    else:
        plt.close()


_C2_REVIEW_XI = (-0.5, -0.25, 0.0, 0.25, 0.5)
_C2_REVIEW_MARKERS = ("o", "s", "^", "D", "v")
_C2_REVIEW_LABELS = (r"$-0.5$", r"$-0.25$", r"$0$", r"$+0.25$", r"$+0.5$")
_C2_REVIEW_SIZE = (6.75, 3.50)
DEFAULT_C2_REVIEW_OUTPUT = (
    PROJECT_ROOT / "data" / "outputs" / "figures" / "pnjl" / "phase_diagram" / "transport_issue130_c2_v1"
)
DEFAULT_C2_REVIEW_RESULTS = (
    PROJECT_ROOT / "data" / "outputs" / "results" / "pnjl" / "phase_diagram" / "transport_issue130_c2_v1"
)


def _c2_review_field(rows, key):
    values = np.asarray([float(r[key]) for r in rows])
    if not np.isfinite(values).all():
        raise ValueError(f"nonfinite values in {key}")
    return values


def _c2_review_axis(name, label, unit, transform, limits):
    return dict(field=name, label=label, source_unit=unit, display_unit=unit,
                transform=transform, scale="linear", limits=limits)


def render_c2_review_figure(*, output_dir: Path = DEFAULT_C2_REVIEW_OUTPUT,
                            results_dir: Path = DEFAULT_C2_REVIEW_RESULTS,
                            source_root: Path = DEFAULT_PHASE_REFERENCE_ROOT,
                            acceptance_path: Path | None = None) -> None:
    """Render the C2 review or export the PDF matching an accepted review."""
    import json
    from pathlib import Path as _Path

    import matplotlib
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator
    from scripts.plotting.publication_review import (
        LABEL_PT, LEGEND_PT, REVIEW_COLORS, TICK_PT, figure_context,
        grayscale_review_record, register_shared_series_styles, write_sop_snapshot,
    )
    from scripts.plotting.plot_accessibility import accessibility_record
    from scripts.plotting.plot_manifest import (
        build_manifest, generator_record, input_record, runtime_record, write_manifest,
    )
    from scripts.plotting.plot_quality import export_figure, placement_limits
    from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
    from scripts.plotting.plot_style import configure_axis_ticks, load_profile
    from scripts.plotting.review_vector_delivery import (
        COLOR_ROUTE, add_vector_pdf, read_accepted_review, require_fresh_case,
    )
    from scripts.plotting.validate_plot_artifact import validate_manifest

    output_dir, results_dir = _Path(output_dir), _Path(results_dir)

    def _repo_relative(path: _Path) -> str:
        try:
            return _Path(path).resolve().relative_to(PROJECT_ROOT).as_posix()
        except ValueError:
            return str(path)
    profile = load_profile("candidate_aps_v2")
    palette = (*REVIEW_COLORS, *profile.colors[3:])
    registry = SeriesStyleRegistry("c2_phase_slice_encoding_v1", tuple(
        SeriesStyle(f"c2.{xi:g}.{field}", palette[index], linestyle, mark)
        for index, (xi, marker) in enumerate(zip(_C2_REVIEW_XI, _C2_REVIEW_MARKERS))
        for field, linestyle, mark in (
            ("first_order", "-", ""), ("crossover", "--", ""), ("endpoint", "None", marker))
    ))
    styles = {style.series_id: style for style in registry.styles}
    require_fresh_case(output_dir, results_dir)
    asset_id = "pnjl.transport.phase.native_issue130_c2"
    acceptance = None
    if acceptance_path is not None:
        acceptance = read_accepted_review(acceptance_path, root=PROJECT_ROOT, generator=Path(__file__))
        if acceptance.record["asset_id"] != asset_id:
            raise ValueError("acceptance must select the C2 phase figure")
    rows, sources, meta = load_c2_review_sources(source_root, _C2_REVIEW_XI)
    sop = write_sop_snapshot(results_dir)
    inputs = [input_record(path, role=role) for path, role in (
        *sources, (profile.path, "APS plot profile"), (sop, "production figure SOP snapshot, LF as in Git"),
        (PROJECT_ROOT / "scripts/pnjl/phase_reference_adapter.py", "accepted C2 row and density adapter"),
    )]
    inputs.extend(input_record(PROJECT_ROOT / f"scripts/plotting/{helper}.py", role="plotting helper")
                  for helper in ("plot_style", "plot_series", "plot_quality", "plot_accessibility", "plot_manifest",
                                 "plot_provenance", "publication_review", "plot_bundle", "plot_delivery",
                                 "review_vector_delivery", "validate_plot_artifact"))
    output_dir.mkdir(parents=True, exist_ok=True)
    series, landmarks, slices = [], [], []
    with figure_context(profile, compact=True) as ctx:
        fig, axes = plt.subplots(1, 2, figsize=_C2_REVIEW_SIZE, sharey=True)
        ctx.bind(fig)
        fig.subplots_adjust(left=.111, right=.974, bottom=.195, top=.915, wspace=.13)
        for index, xi in enumerate(_C2_REVIEW_XI):
            def take(kind, key):
                picked = sorted((r for r in rows[kind] if float(r["xi"]) == xi),
                                key=lambda r: float(r[key]))
                if not picked or len(picked) != len(set(_c2_review_field(picked, key))):
                    raise ValueError(f"C2 {kind} slice at xi={xi} is empty or has duplicate {key}")
                return picked

            boundary = take("boundary", "T_MeV")
            crossover = take("crossover", "mu_MeV")
            endpoint_rows = take("cep", "xi")
            if len(endpoint_rows) != 1:
                raise ValueError(f"C2 CEP slice at xi={xi} must contain exactly one row")
            endpoint = endpoint_rows[0]
            low, high, t_crit = (float(endpoint[k]) for k in ("T_low_MeV", "T_high_MeV", "T_midpoint_MeV"))
            mu_crit = 3 * float(endpoint["mu_CEP_proxy_MeV"])
            if not np.isclose((low + high) / 2, t_crit) or high <= low:
                raise ValueError(f"C2 CEP bracket at xi={xi} is inconsistent")
            temperature, mu, rho_low, rho_high = (
                _c2_review_field(boundary, k) for k in ("T_MeV", "mu_MeV", "rho_hadron", "rho_quark")
            )
            crossover_t, crossover_mu, crossover_rho = (
                _c2_review_field(crossover, k) for k in ("T_MeV", "mu_MeV", "rho")
            )
            if (temperature[-1] != low or max(crossover_mu) * 3 >= mu_crit
                    or min(crossover_t) <= t_crit or not np.all(rho_low < rho_high)):
                raise ValueError(f"C2 slice at xi={xi} violates the declared support rules")
            density = estimate_c2_endpoint_density(temperature, rho_low, rho_high, t_crit)
            color = palette[index]
            for panel, thing, xx, yy, support, max_gap, linestyle in (
                (0, "Maxwell", 3 * mu, temperature, temperature, 2.01, "-"),
                (0, "crossover", 3 * crossover_mu, crossover_t, crossover_mu, np.inf, "--"),
                (1, "low-density", rho_low, temperature, temperature, 2.01, "-"),
                (1, "high-density", rho_high, temperature, temperature, 2.01, "-"),
                (1, "crossover", crossover_rho, crossover_t, crossover_mu, np.inf, "--"),
            ):
                segments = np.split(np.arange(len(yy)), np.flatnonzero(np.diff(support) > max_gap) + 1)
                series_id = f"{panel}_{thing}_{xi:g}"
                style_id = f"c2.{xi:g}.{'crossover' if linestyle == '--' else 'first_order'}"
                for chunk in segments:
                    line, = axes[panel].plot(xx[chunk], yy[chunk], lw=1.2, **styles[style_id].plot_kwargs())
                    if not np.array_equal(line.get_xydata(), np.column_stack((xx[chunk], yy[chunk]))):
                        raise ValueError("rendered vertices differ from the selected native rows")
                series.append(dict(
                    series_id=series_id, style_ids=[style_id], state="native_C2", xi=xi, panel=panel,
                    support_rule="exact xi, native and certified; order by T for coexistence and quark mu for crossover",
                    mask_rule="no omitted vertices; break native Maxwell gaps over 2 MeV; view T >= 60 MeV",
                    input_rows=len(yy), visible_rows=int(np.count_nonzero(yy >= 60)),
                    segments=len(segments), x=xx.tolist(), y=yy.tolist(),
                ))
            style_id = f"c2.{xi:g}.endpoint"
            for panel, x_cep, face in ((0, mu_crit, color), (1, density["rho_over_rho0"], "white")):
                axes[panel].plot(x_cep, t_crit, markerfacecolor=face, markeredgewidth=1.05,
                                 markersize=6.1, zorder=5, **styles[style_id].plot_kwargs())
            landmarks.append(dict(xi=xi, T_low=low, T_high=high, T_midpoint_MeV=t_crit,
                                  muB_proxy_MeV=mu_crit, crossover_last_T_MeV=crossover_t[-1], **density))
            series.extend([
                dict(series_id=f"endpoint_mu_{xi:g}", style_ids=[style_id], state="estimated_midpoint", xi=xi, T_MeV=t_crit,
                     muB_proxy_MeV=mu_crit,
                     support_rule="C2 bracket midpoint and the table's quark mu proxy times three",
                     mask_rule="none"),
                dict(series_id=f"endpoint_density_{xi:g}", style_ids=[style_id], xi=xi,
                     support_rule="C2 nearest subcritical coexistence pair", mask_rule="none", **density),
            ])
            slices.append(dict(xi=xi, maxwell_rows=len(temperature), crossover_rows=len(crossover_t),
                               maxwell_rows_in_view=int(sum(temperature >= 60))))
        xlabels = (r"$\mu_B$ (MeV)", r"$\rho/\rho_0$")
        for index, axis in enumerate(axes):
            axis.set_xlim(0, (1250, 3.15)[index])
            axis.set_ylim(60, 240)
            axis.xaxis.set_major_locator(FixedLocator(([0, 400, 800, 1200], [0, 1, 2, 3])[index]))
            axis.yaxis.set_major_locator(FixedLocator([60, 100, 140, 180, 220]))
            axis.set_xlabel(xlabels[index])
            configure_axis_ticks(axis, profile)
            axis.text(0, 1.025, "(a)" if index == 0 else "(b)", transform=axis.transAxes,
                      va="bottom", ha="left", fontsize=TICK_PT)
        axes[0].set_ylabel(r"$T$ (MeV)")
        handles = [Line2D([], [], marker=m, color=palette[i], ls="none", ms=6)
                   for i, m in enumerate(_C2_REVIEW_MARKERS)]
        common = dict(fontsize=LEGEND_PT, title_fontsize=LABEL_PT, frameon=False, handletextpad=.32,
                      labelspacing=.12, borderaxespad=.35, borderpad=.1, alignment="left")
        cep_title = r"CEP at each $\xi$"
        cep_key = axes[0].legend(handles, _C2_REVIEW_LABELS, title=cep_title, ncols=2,
                                 loc="lower left", columnspacing=.45, handlelength=.8, **common)
        axes[0].add_artist(cep_key)
        line_types = ["First order", "Crossover"]
        axes[0].legend([Line2D([], [], ls=style, color="#383838", lw=1.2) for style in ("-", "--")],
                       line_types, loc="upper right", handlelength=1.35, **common)
        placements = [
            dict(legend_index=0, host_axes_index=0, title=cep_title, labels=list(_C2_REVIEW_LABELS),
                 placement="in_axes", location="lower left",
                 scope="CEP symbols for both panels; fill style is panel-specific"),
            dict(legend_index=1, host_axes_index=0, title="", labels=line_types, placement="in_axes",
                 location="upper right", scope="curve line styles in both panels"),
        ]
        stem = output_dir / "phase_diagram_TmuB_Trho"
        outputs, quality = export_figure(fig, stem, profile, formats=["png"])
    outputs[0]["role"] = "color_review"
    outputs.append(grayscale_review_record(stem.with_suffix(".png")))
    points = results_dir / "selected_slices.json"
    write_manifest(points, dict(source_id=meta["source_run_id"], counts=slices,
                                endpoint_estimates=landmarks,
                                chemical_potential="muB = 3 * mu_quark, including the endpoint proxy",
                                density="source rho, rho_hadron, rho_quark already express nB/rho0",
                                axis_order="never sort by density"), overwrite=True, newline="\n")
    inputs.append(input_record(points, role="slice and endpoint audit"))
    rendering = dict(delivery_stage="png_review", vector_delivery_pending=True, review_contract="prd_review_v1",
                     column="double_column", figure_size_inches=list(_C2_REVIEW_SIZE),
                     size_override_reason="two side-by-side phase planes with in-panel keys", quality=quality,
                     series_styles=register_shared_series_styles(series, registry),
                     color_route=COLOR_ROUTE,
                     accessibility=accessibility_record(
                         [dict(value=str(v), color=palette[i], noncolor_encoding=f"endpoint marker {_C2_REVIEW_MARKERS[i]}")
                          for i, v in enumerate(_C2_REVIEW_XI)],
                         backgrounds=["#ffffff"], background_scope="white axes"),
                     legend_policy="declared_geometry_checked", legend_outside=False,
                     legend_placements=placements,
                     composite_style_reference="Times/STIX, 13 pt axis labels, 11 pt ticks/keys, 1.2 pt curves",
                     palette=list(palette[:5]), placement_limits=placement_limits(quality, profile, outputs))
    record = build_manifest(
        asset_id=asset_id, figure_family="phase_diagram",
        case_slug=output_dir.name, figure_mode="audit",
        semantic_status="native C2 slices with distinct estimated endpoint markers",
        style_profile=profile.profile_id, publication_scope="internal_review",
        calculation_sha=meta["calculation_sha"], source_run_id=meta["source_run_id"],
        generator=generator_record(
            __file__,
            command=("python scripts/pnjl/plot_phase_diagram.py --review-c2 "
                     f"--output-dir {_repo_relative(output_dir)} --results-dir {_repo_relative(results_dir)}"
                     + (f" --deliver-vector --acceptance {_repo_relative(acceptance_path)}"
                        if acceptance_path else "")),
            runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=inputs, outputs=outputs, series=series,
        axes=[_c2_review_axis("muB", xlabels[0], "MeV", "muB=3*muq", [0, 1250]),
              _c2_review_axis("rho_over_rho0", xlabels[1], "dimensionless", "identity", [0, 3.15]),
              _c2_review_axis("T_MeV", r"$T$ (MeV)", "MeV", "identity", [60, 240])],
        selection_rule="same five native xi slices as the previous C1 figure, from the C2-based accepted map",
        interpolation_policy="none; straight lines between adjacent native samples only",
        connector_policy="forbidden",
        missing_value_policy="no gap closure, extrapolation, or replacement of native points",
        validation=dict(finite=True, duplicate_keys=True, support=True,
                        all_selected_source_rows_native_certified=True, solver_called=False),
        rendering=rendering)
    record.update(delivery_stage="png_review", vector_delivery_pending=True,
                  manuscript_eligible=False, current_publication_layer=False)
    try:
        if acceptance is not None:
            add_vector_pdf(record, acceptance, figure=fig, stem=stem, profile=profile, root=PROJECT_ROOT)
    finally:
        plt.close(fig)
    manifest = output_dir / "plot_manifest.json"
    write_manifest(manifest, record, newline="\n")
    errors = validate_manifest(manifest)
    if errors:
        raise ValueError(errors)
    print(json.dumps({"slices": slices, "minimum_capital_numeral_height_mm":
                      quality["minimum_capital_numeral_height_mm"]}, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot PNJL phase diagrams",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--boundary", type=Path, default=DEFAULT_BOUNDARY_PATH,
                       help="Boundary data file")
    parser.add_argument("--cep", type=Path, default=DEFAULT_CEP_PATH,
                       help="CEP data file")
    parser.add_argument("--spinodal", type=Path, default=DEFAULT_SPINODAL_PATH,
                       help="Spinodal data file")
    parser.add_argument("--crossover", type=Path, default=DEFAULT_CROSSOVER_PATH,
                       help="Crossover data file")
    parser.add_argument(
        "--phase-reference-root",
        type=Path,
        default=None,
        help="Issue #130 candidate root；未指定时使用仓库内 v2 accepted layer",
    )
    parser.add_argument(
        "--phase-reference-layer",
        choices=["strict", "derived", "render", "accepted"],
        default=None,
        help="candidate layer; omitted selects accepted for v2 and strict for v1",
    )
    parser.add_argument("--xi", type=float, action="append", default=None,
                       help="ξ values to plot (can specify multiple times)")
    parser.add_argument("--output-dir", type=Path, default=None,
                       help="Output directory (quick plot defaults to data/outputs/figures/pnjl)")
    parser.add_argument("--format", type=str, default="png", choices=["png", "pdf", "svg"],
                       help="Output format")
    parser.add_argument("--dpi", type=int, default=200, help="DPI for raster formats")
    parser.add_argument("--no-show", action="store_true", help="Don't show plot window")
    parser.add_argument("--type", type=str, default="all", 
                       choices=["T-mu", "T-rho", "combined", "all"],
                       help="Type of phase diagram to plot")
    parser.add_argument("--no-spinodal", action="store_true", 
                       help="Don't plot spinodal lines")
    parser.add_argument("--no-crossover", action="store_true",
                       help="Don't plot crossover lines")
    parser.add_argument("--review-c2", action="store_true",
                       help="Render the accepted C2 review figure with an SOP manifest instead of a quick plot")
    parser.add_argument("--results-dir", type=Path, default=None,
                       help="Metadata directory used by --review-c2 (must accompany --output-dir)")
    
    from scripts.plotting.review_vector_delivery import add_vector_arguments, require_paired_vector_arguments

    add_vector_arguments(parser)
    args = parser.parse_args()
    require_paired_vector_arguments(parser, args)
    if args.deliver_vector and not args.review_c2:
        parser.error("--deliver-vector requires --review-c2")
    if args.deliver_vector and (args.output_dir is None or args.results_dir is None):
        parser.error("vector delivery requires explicit --output-dir and --results-dir")

    if args.review_c2:
        ignored = [name for name, value, default in (
            ("--format", args.format, "png"),
            ("--dpi", args.dpi, 200),
            ("--type", args.type, "all"),
            ("--xi", args.xi, None),
            ("--phase-reference-layer", args.phase_reference_layer, None),
        ) if value != default]
        if args.no_spinodal or args.no_crossover:
            ignored.append("--no-spinodal/--no-crossover")
        if ignored:
            parser.error("--review-c2 does not accept " + ", ".join(ignored))
        if (args.output_dir is None) != (args.results_dir is None):
            parser.error("--review-c2 requires --output-dir and --results-dir together")
        render_c2_review_figure(
            output_dir=(args.output_dir or DEFAULT_C2_REVIEW_OUTPUT).resolve(),
            results_dir=(args.results_dir or DEFAULT_C2_REVIEW_RESULTS).resolve(),
            source_root=(args.phase_reference_root or DEFAULT_PHASE_REFERENCE_ROOT),
            acceptance_path=args.acceptance,
        )
        return
    
    # 默认加载作者接受的 v2 phase map。旧 CSV 仍可通过显式路径回放。
    candidate_diagnostics = None
    use_default_candidate = (
        args.phase_reference_root is None
        and DEFAULT_PHASE_REFERENCE_ROOT.is_dir()
        and args.boundary == DEFAULT_BOUNDARY_PATH
        and args.cep == DEFAULT_CEP_PATH
        and args.spinodal == DEFAULT_SPINODAL_PATH
        and args.crossover == DEFAULT_CROSSOVER_PATH
    )
    if args.phase_reference_root is not None or use_default_candidate:
        reference_root = args.phase_reference_root or DEFAULT_PHASE_REFERENCE_ROOT
        (
            boundary_points,
            cep_points,
            spinodal_points,
            crossover_points,
            candidate_diagnostics,
        ) = load_candidate_phase_data(
            reference_root,
            args.phase_reference_layer or DEFAULT_PHASE_REFERENCE_LAYER,
        )
        if args.no_spinodal:
            spinodal_points = []
        if args.no_crossover:
            crossover_points = []
    else:
        boundary_points = load_boundary_data(args.boundary)
        cep_points = load_cep_data(args.cep)
        spinodal_points = [] if args.no_spinodal else load_spinodal_data(args.spinodal)
        crossover_points = [] if args.no_crossover else load_crossover_data(args.crossover)
    
    if not boundary_points:
        print(f"No boundary data found in {args.boundary}")
        return
    
    boundary_groups = group_by_xi(boundary_points)
    xi_filter = args.xi
    
    print(f"Loaded {len(boundary_points)} boundary points")
    print(f"Available ξ values: {sorted(boundary_groups.keys())}")
    print(f"CEP points: {len(cep_points)}")
    print(f"Spinodal points: {len(spinodal_points)}")
    print(f"Crossover points: {len(crossover_points)}")
    if candidate_diagnostics is not None:
        print(
            "Candidate adapter: "
            f"layer={candidate_diagnostics['layer']}, "
            f"uncertified_rows={candidate_diagnostics['uncertified_rows']}, "
            "runtime_consumption=false"
        )
    
    show = not args.no_show
    output_dir = args.output_dir or DEFAULT_OUTPUT_DIR
    fmt = args.format
    dpi = args.dpi
    
    # 生成文件名后缀
    if xi_filter and len(xi_filter) == 1:
        suffix = f"_xi{xi_filter[0]:.1f}"
    else:
        suffix = ""
    
    # 绘制相图
    if args.type in ["T-mu", "all"]:
        plot_T_mu_phase_diagram(
            boundary_groups, cep_points, spinodal_points, crossover_points, xi_filter,
            output_path=output_dir / f"phase_diagram_T_mu{suffix}.{fmt}",
            show=show, dpi=dpi,
        )
    
    if args.type in ["T-rho", "all"]:
        plot_T_rho_phase_diagram(
            boundary_groups, cep_points, spinodal_points, crossover_points, xi_filter,
            output_path=output_dir / f"phase_diagram_T_rho{suffix}.{fmt}",
            show=show, dpi=dpi,
        )
    
    if args.type in ["combined", "all"]:
        plot_combined_phase_diagram(
            boundary_groups, cep_points, spinodal_points, crossover_points, xi_filter,
            output_path=output_dir / f"phase_diagram_combined{suffix}.{fmt}",
            show=show, dpi=dpi,
        )


if __name__ == "__main__":
    main()
