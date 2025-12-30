#!/usr/bin/env python3
"""
PNJL 相结构数据验证脚本

检查生成的相图数据是否光滑、合理：
1. 检测数据跳变（相邻点变化过大）
2. 检测单调性异常
3. 输出统计信息和问题报告

用法：
    python scripts/pnjl/validate_phase_data.py [options]
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np


def _find_project_root() -> Path:
    """Find project root by looking for Project.toml or .git."""
    script_dir = Path(__file__).resolve().parent
    for start in [script_dir, script_dir.parent, script_dir.parent.parent, Path.cwd()]:
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
DEFAULT_SPINODAL_PATH = PROJECT_ROOT / "data" / "reference" / "pnjl" / "spinodals.csv"


@dataclass
class ValidationIssue:
    """数据验证问题"""
    severity: str  # "warning" or "error"
    file: str
    column: str
    xi: float
    T_MeV: float
    message: str
    value: Optional[float] = None
    expected_range: Optional[Tuple[float, float]] = None


def load_csv_data(path: Path) -> Tuple[List[str], List[Dict[str, float]]]:
    """加载 CSV 数据"""
    if not path.exists():
        return [], []
    
    rows = []
    headers = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        for row in reader:
            try:
                parsed = {k: float(v) for k, v in row.items()}
                rows.append(parsed)
            except ValueError:
                continue
    return headers, rows


def group_by_xi(rows: List[Dict[str, float]]) -> Dict[float, List[Dict[str, float]]]:
    """按 xi 分组并按温度排序"""
    groups: Dict[float, List[Dict[str, float]]] = {}
    for row in rows:
        xi = row.get("xi", 0.0)
        groups.setdefault(xi, []).append(row)
    for xi in groups:
        groups[xi].sort(key=lambda r: r.get("T_MeV", 0))
    return groups


def check_smoothness(
    values: List[float], 
    T_values: List[float],
    threshold_ratio: float = 0.5
) -> List[Tuple[int, float, float, float]]:
    """
    检查数据光滑性
    
    返回跳变点列表: [(index, T, value, jump_ratio), ...]
    jump_ratio = |Δvalue| / |value_range|
    """
    if len(values) < 3:
        return []
    
    issues = []
    value_range = max(values) - min(values)
    if value_range < 1e-9:
        return []
    
    for i in range(1, len(values) - 1):
        # 计算局部变化
        prev_diff = values[i] - values[i-1]
        next_diff = values[i+1] - values[i]
        
        # 检测方向突变（符号变化且幅度大）
        if prev_diff * next_diff < 0:
            jump = abs(prev_diff) + abs(next_diff)
            jump_ratio = jump / value_range
            if jump_ratio > threshold_ratio:
                issues.append((i, T_values[i], values[i], jump_ratio))
    
    return issues


def check_monotonicity(
    values: List[float],
    T_values: List[float],
    expected_direction: str = "any"  # "increasing", "decreasing", "any"
) -> List[Tuple[int, float, str]]:
    """
    检查单调性
    
    返回单调性异常点: [(index, T, direction_change), ...]
    """
    if len(values) < 3:
        return []
    
    issues = []
    
    # 确定主要趋势
    diffs = [values[i+1] - values[i] for i in range(len(values)-1)]
    pos_count = sum(1 for d in diffs if d > 0)
    neg_count = sum(1 for d in diffs if d < 0)
    
    if expected_direction == "any":
        # 自动检测主要方向
        if pos_count > neg_count * 2:
            expected_direction = "increasing"
        elif neg_count > pos_count * 2:
            expected_direction = "decreasing"
        else:
            return []  # 无明显趋势
    
    for i, diff in enumerate(diffs):
        if expected_direction == "increasing" and diff < -1e-9:
            issues.append((i+1, T_values[i+1], "unexpected decrease"))
        elif expected_direction == "decreasing" and diff > 1e-9:
            issues.append((i+1, T_values[i+1], "unexpected increase"))
    
    return issues


def validate_boundary_data(path: Path) -> List[ValidationIssue]:
    """验证 boundary.csv 数据"""
    issues = []
    headers, rows = load_csv_data(path)
    
    if not rows:
        issues.append(ValidationIssue(
            severity="error", file=str(path), column="",
            xi=0, T_MeV=0, message="No data found"
        ))
        return issues
    
    groups = group_by_xi(rows)
    
    for xi, data in groups.items():
        T_vals = [r["T_MeV"] for r in data]
        
        # 检查各列
        columns_to_check = [
            ("mu_transition_MeV", "decreasing"),  # μ 应随 T 增加而减小
            ("mu_MeV", "decreasing"),             # 兼容旧命名
            ("mu_coex_MeV", "decreasing"),        # 兼容旧命名
            ("rho_hadron", "increasing"),         # ρ_hadron 应随 T 增加而增加
            ("rho_quark", "decreasing"),          # ρ_quark 应随 T 增加而减小
            # 兼容旧命名
            ("rho_gas", "increasing"),
            ("rho_liquid", "decreasing"),
        ]
        
        for col, expected_mono in columns_to_check:
            if col not in headers:
                continue
            
            values = [r.get(col, float('nan')) for r in data]
            if any(np.isnan(v) for v in values):
                continue
            
            # 光滑性检查
            jumps = check_smoothness(values, T_vals, threshold_ratio=0.3)
            for idx, T, val, ratio in jumps:
                issues.append(ValidationIssue(
                    severity="warning",
                    file="boundary.csv",
                    column=col,
                    xi=xi,
                    T_MeV=T,
                    message=f"Large jump detected (ratio={ratio:.2f})",
                    value=val
                ))
            
            # 单调性检查
            mono_issues = check_monotonicity(values, T_vals, expected_mono)
            for idx, T, msg in mono_issues:
                issues.append(ValidationIssue(
                    severity="warning",
                    file="boundary.csv",
                    column=col,
                    xi=xi,
                    T_MeV=T,
                    message=f"Monotonicity issue: {msg}",
                    value=values[idx]
                ))
    
    return issues


def validate_spinodal_data(path: Path) -> List[ValidationIssue]:
    """验证 spinodals.csv 数据"""
    issues = []
    headers, rows = load_csv_data(path)
    
    if not rows:
        issues.append(ValidationIssue(
            severity="error", file=str(path), column="",
            xi=0, T_MeV=0, message="No data found"
        ))
        return issues
    
    groups = group_by_xi(rows)
    
    for xi, data in groups.items():
        T_vals = [r["T_MeV"] for r in data]
        
        # 检查各列
        columns_to_check = [
            ("mu_spinodal_hadron_MeV", "decreasing"),
            ("mu_spinodal_quark_MeV", "decreasing"),
            ("mu_spinodal_low_MeV", "decreasing"),  # 兼容旧命名
            ("mu_spinodal_high_MeV", "decreasing"),
            ("rho_spinodal_hadron", "any"),
            ("rho_spinodal_quark", "any"),
            # 兼容旧命名
            ("rho_spinodal_low", "any"),
            ("rho_spinodal_high", "any"),
        ]
        
        for col, expected_mono in columns_to_check:
            if col not in headers:
                continue
            
            values = [r.get(col, float('nan')) for r in data]
            if any(np.isnan(v) for v in values):
                continue
            
            # 光滑性检查（spinodal 数据更严格）
            jumps = check_smoothness(values, T_vals, threshold_ratio=0.25)
            for idx, T, val, ratio in jumps:
                issues.append(ValidationIssue(
                    severity="warning" if ratio < 0.5 else "error",
                    file="spinodals.csv",
                    column=col,
                    xi=xi,
                    T_MeV=T,
                    message=f"Large jump detected (ratio={ratio:.2f})",
                    value=val
                ))
    
    return issues


def print_statistics(path: Path, name: str) -> None:
    """打印数据统计信息"""
    headers, rows = load_csv_data(path)
    if not rows:
        print(f"\n{name}: No data")
        return
    
    groups = group_by_xi(rows)
    
    print(f"\n{name}:")
    print(f"  Total rows: {len(rows)}")
    print(f"  ξ values: {sorted(groups.keys())}")
    
    for xi in sorted(groups.keys()):
        data = groups[xi]
        T_min = min(r["T_MeV"] for r in data)
        T_max = max(r["T_MeV"] for r in data)
        print(f"  ξ={xi}: {len(data)} points, T=[{T_min:.1f}, {T_max:.1f}] MeV")
        
        # 打印各列的范围
        for col in headers:
            if col in ("xi", "T_MeV"):
                continue
            values = [r.get(col, float('nan')) for r in data]
            values = [v for v in values if not np.isnan(v)]
            if values:
                print(f"    {col}: [{min(values):.3f}, {max(values):.3f}]")


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate PNJL phase data")
    parser.add_argument("--boundary", type=Path, default=DEFAULT_BOUNDARY_PATH)
    parser.add_argument("--spinodal", type=Path, default=DEFAULT_SPINODAL_PATH)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    
    print("=" * 60)
    print("PNJL Phase Data Validation")
    print("=" * 60)
    
    # 统计信息
    print_statistics(args.boundary, "boundary.csv")
    print_statistics(args.spinodal, "spinodals.csv")
    
    # 验证
    print("\n" + "-" * 60)
    print("Validation Results")
    print("-" * 60)
    
    all_issues = []
    all_issues.extend(validate_boundary_data(args.boundary))
    all_issues.extend(validate_spinodal_data(args.spinodal))
    
    if not all_issues:
        print("\n✓ All data passed validation!")
    else:
        errors = [i for i in all_issues if i.severity == "error"]
        warnings = [i for i in all_issues if i.severity == "warning"]
        
        print(f"\nFound {len(errors)} errors, {len(warnings)} warnings")
        
        for issue in all_issues:
            symbol = "✗" if issue.severity == "error" else "⚠"
            val_str = f" (value={issue.value:.4f})" if issue.value is not None else ""
            print(f"  {symbol} [{issue.file}] ξ={issue.xi}, T={issue.T_MeV:.1f} MeV, "
                  f"{issue.column}: {issue.message}{val_str}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
