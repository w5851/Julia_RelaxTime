#!/usr/bin/env python3
"""
PNJL 相图绘制脚本

从 data/reference/pnjl/ 读取相结构数据，绘制：
1. T-μ 相图（一阶相变线 + CEP）
2. T-ρ 相图（共存区）

用法：
    python scripts/pnjl/plot_phase_diagram.py [options]

选项：
    --boundary PATH   相变线数据文件 (默认 data/reference/pnjl/boundary.csv)
    --cep PATH        CEP 数据文件 (默认 data/reference/pnjl/cep.csv)
    --xi VALUE        要绘制的 ξ 值，可多次指定 (默认绘制所有)
    --output-dir DIR  输出目录 (默认 data/outputs/figures/pnjl)
    --format FMT      输出格式 png/pdf/svg (默认 png)
    --dpi VALUE       分辨率 (默认 200)
    --no-show         不显示图形窗口
"""

from __future__ import annotations

import argparse
import csv
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
DEFAULT_BOUNDARY_PATH = PROJECT_ROOT / "data" / "reference" / "pnjl" / "boundary.csv"
DEFAULT_CEP_PATH = PROJECT_ROOT / "data" / "reference" / "pnjl" / "cep.csv"
DEFAULT_SPINODAL_PATH = PROJECT_ROOT / "data" / "reference" / "pnjl" / "spinodals.csv"
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
    xi_filter: Optional[List[float]] = None,
    output_path: Optional[Path] = None,
    show: bool = True,
    dpi: int = 200,
) -> None:
    """绘制 T-ρ 相图（共存区）"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    xi_values = sorted(boundary_groups.keys())
    if xi_filter:
        xi_values = [xi for xi in xi_values if xi in xi_filter]
    
    spinodal_groups = group_spinodals_by_xi(spinodals) if spinodals else {}
    
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
    parser.add_argument("--xi", type=float, action="append", default=None,
                       help="ξ values to plot (can specify multiple times)")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                       help="Output directory")
    parser.add_argument("--format", type=str, default="png", choices=["png", "pdf", "svg"],
                       help="Output format")
    parser.add_argument("--dpi", type=int, default=200, help="DPI for raster formats")
    parser.add_argument("--no-show", action="store_true", help="Don't show plot window")
    parser.add_argument("--type", type=str, default="all", 
                       choices=["T-mu", "T-rho", "combined", "all"],
                       help="Type of phase diagram to plot")
    parser.add_argument("--no-spinodal", action="store_true", 
                       help="Don't plot spinodal lines")
    
    args = parser.parse_args()
    
    # 加载数据
    boundary_points = load_boundary_data(args.boundary)
    cep_points = load_cep_data(args.cep)
    spinodal_points = [] if args.no_spinodal else load_spinodal_data(args.spinodal)
    
    if not boundary_points:
        print(f"No boundary data found in {args.boundary}")
        return
    
    boundary_groups = group_by_xi(boundary_points)
    xi_filter = args.xi
    
    print(f"Loaded {len(boundary_points)} boundary points")
    print(f"Available ξ values: {sorted(boundary_groups.keys())}")
    print(f"CEP points: {len(cep_points)}")
    print(f"Spinodal points: {len(spinodal_points)}")
    
    show = not args.no_show
    output_dir = args.output_dir
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
            boundary_groups, cep_points, spinodal_points, xi_filter,
            output_path=output_dir / f"phase_diagram_T_mu{suffix}.{fmt}",
            show=show, dpi=dpi,
        )
    
    if args.type in ["T-rho", "all"]:
        plot_T_rho_phase_diagram(
            boundary_groups, cep_points, spinodal_points, xi_filter,
            output_path=output_dir / f"phase_diagram_T_rho{suffix}.{fmt}",
            show=show, dpi=dpi,
        )
    
    if args.type in ["combined", "all"]:
        plot_combined_phase_diagram(
            boundary_groups, cep_points, spinodal_points, xi_filter,
            output_path=output_dir / f"phase_diagram_combined{suffix}.{fmt}",
            show=show, dpi=dpi,
        )


if __name__ == "__main__":
    main()
