"""Load and apply the repository plotting style profiles.

This module owns visual defaults only. It does not select physical rows,
interpolate data, call solvers, or modify numerical inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import tomllib


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROFILE_ROOT = PROJECT_ROOT / "config" / "plotting"
PROFILE_FILES = {
    "audit_v1": PROFILE_ROOT / "audit_v1.toml",
    "candidate_origin_like_v1": PROFILE_ROOT / "candidate_origin_like_v1.toml",
    "strict_origin_like_v1": PROFILE_ROOT / "strict_origin_like_v1.toml",
}
ALLOWED_COLUMNS = {"single_column", "double_column"}
ALLOWED_PROFILES = frozenset(PROFILE_FILES)


@dataclass(frozen=True)
class PlotProfile:
    """Validated, immutable view of a TOML plotting profile."""

    profile_id: str
    data: Mapping[str, Any]
    path: Path

    @property
    def default_column(self) -> str:
        return str(self.data["default_column"])

    @property
    def formats(self) -> tuple[str, ...]:
        return tuple(str(item).lower() for item in self.data["output"]["formats"])

    @property
    def optional_formats(self) -> tuple[str, ...]:
        return tuple(str(item).lower() for item in self.data["output"].get("optional_formats", []))

    @property
    def dpi(self) -> int:
        return int(self.data["output"]["dpi"])

    @property
    def colors(self) -> tuple[str, ...]:
        return tuple(str(item) for item in self.data["palette"]["colors"])

    @property
    def draw_support_points(self) -> bool:
        return bool(self.data.get("support", {}).get("draw_points", True))

    @property
    def support_visual_policy(self) -> str:
        return str(self.data.get("support", {}).get("visual_policy", "unspecified"))

    @property
    def legend_policy(self) -> str:
        return str(self.data.get("layout", {}).get("legend_policy", "unspecified"))

    @property
    def max_in_axes_legend_entries(self) -> int:
        return int(self.data.get("layout", {}).get("max_in_axes_entries", 0))

    @property
    def allow_external_legend(self) -> bool:
        return bool(self.data.get("layout", {}).get("allow_external_legend", False))


def _profile_path(profile: str | Path) -> Path:
    if isinstance(profile, Path):
        return profile if profile.is_absolute() else PROJECT_ROOT / profile
    if profile in PROFILE_FILES:
        return PROFILE_FILES[profile]
    candidate = Path(profile)
    if candidate.suffix == ".toml":
        return candidate if candidate.is_absolute() else PROJECT_ROOT / candidate
    raise ValueError(f"unknown plotting profile: {profile}")


def _validate_profile(data: Mapping[str, Any], path: Path) -> str:
    required = {"schema_version", "profile_id", "default_column", "figure_size_in", "output", "semantics", "palette"}
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"{path}: missing profile fields: {', '.join(missing)}")
    if data["schema_version"] != "plot_style_profile_v1":
        raise ValueError(f"{path}: unsupported schema_version={data['schema_version']!r}")
    profile_id = str(data["profile_id"])
    if profile_id not in ALLOWED_PROFILES:
        raise ValueError(f"{path}: unsupported profile_id={profile_id!r}")
    if str(data["default_column"]) not in ALLOWED_COLUMNS:
        raise ValueError(f"{path}: invalid default_column")
    sizes = data["figure_size_in"]
    for column in ALLOWED_COLUMNS:
        value = sizes.get(column)
        if not isinstance(value, list) or len(value) != 2 or any(float(item) <= 0 for item in value):
            raise ValueError(f"{path}: figure_size_in.{column} must be a positive [width, height] pair")
    output = data["output"]
    formats = [str(item).lower() for item in output.get("formats", [])]
    if not formats or len(formats) != len(set(formats)):
        raise ValueError(f"{path}: output.formats must be a non-empty unique list")
    if int(output.get("dpi", 0)) <= 0:
        raise ValueError(f"{path}: output.dpi must be positive")
    colors = data["palette"].get("colors", [])
    if not colors:
        raise ValueError(f"{path}: palette.colors must not be empty")
    support = data.get("support")
    if not isinstance(support, Mapping) or not isinstance(support.get("draw_points"), bool):
        raise ValueError(f"{path}: support.draw_points must be a boolean")
    if not str(support.get("visual_policy", "")):
        raise ValueError(f"{path}: support.visual_policy must be non-empty")
    layout = data.get("layout")
    if not isinstance(layout, Mapping):
        raise ValueError(f"{path}: layout must be a table")
    if not str(layout.get("legend_policy", "")):
        raise ValueError(f"{path}: layout.legend_policy must be non-empty")
    if not isinstance(layout.get("max_in_axes_entries"), int) or int(layout["max_in_axes_entries"]) <= 0:
        raise ValueError(f"{path}: layout.max_in_axes_entries must be a positive integer")
    if not isinstance(layout.get("allow_external_legend"), bool):
        raise ValueError(f"{path}: layout.allow_external_legend must be a boolean")
    return profile_id


def load_profile(profile: str | Path) -> PlotProfile:
    """Load one repository profile and validate its public shape."""

    path = _profile_path(profile)
    if not path.is_file():
        raise FileNotFoundError(f"plotting profile not found: {path}")
    with path.open("rb") as handle:
        data = tomllib.load(handle)
    profile_id = _validate_profile(data, path)
    return PlotProfile(profile_id=profile_id, data=data, path=path)


def figure_size_inches(profile: PlotProfile, column: str | None = None) -> tuple[float, float]:
    """Return the configured physical size in inches."""

    selected = column or profile.default_column
    if selected not in ALLOWED_COLUMNS:
        raise ValueError(f"invalid plotting column: {selected}")
    width, height = profile.data["figure_size_in"][selected]
    return float(width), float(height)


def phase_style(profile: PlotProfile, state: str) -> dict[str, Any]:
    """Return the line/marker contract for a physical state."""

    styles = profile.data.get("phase_styles", {})
    if state not in styles:
        raise ValueError(f"unknown phase state: {state}")
    return dict(styles[state])


def resolve_font(profile: PlotProfile) -> dict[str, str]:
    """Resolve the requested font for provenance without changing the profile."""

    from matplotlib import font_manager

    candidates = [str(item) for item in profile.data.get("font_candidates", [])]
    requested = candidates[0] if candidates else "serif"
    path = font_manager.findfont(font_manager.FontProperties(family=requested), fallback_to_default=True)
    resolved = font_manager.FontProperties(fname=path).get_name()
    return {"requested": requested, "resolved": resolved, "path": str(Path(path).resolve())}


def configure_matplotlib(profile: PlotProfile) -> dict[str, str]:
    """Apply a profile to matplotlib and return resolved font provenance."""

    import matplotlib
    from cycler import cycler

    font_candidates = [str(item) for item in profile.data.get("font_candidates", [])]
    font_size = float(profile.data.get("font_size_pt", 9.0))
    line_width = float(profile.data.get("line_width_pt", 1.0))
    axes_linewidth = float(profile.data.get("axes_linewidth_pt", 0.6))
    matplotlib.rcParams.update(
        {
            "font.family": profile.data.get("font_family", "serif"),
            "font.serif": font_candidates,
            "font.size": font_size,
            "mathtext.fontset": profile.data.get("mathtext_fontset", "stix"),
            "axes.prop_cycle": cycler(color=list(profile.colors)),
            "axes.linewidth": axes_linewidth,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size,
            "legend.fontsize": max(7.0, font_size - 1.0),
            "legend.frameon": False,
            "xtick.labelsize": max(7.0, font_size - 1.0),
            "ytick.labelsize": max(7.0, font_size - 1.0),
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.width": max(0.5, axes_linewidth),
            "ytick.major.width": max(0.5, axes_linewidth),
            "lines.linewidth": line_width,
            "lines.markersize": float(profile.data.get("marker_size_pt", 4.0)),
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": profile.dpi,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )
    return resolve_font(profile)
