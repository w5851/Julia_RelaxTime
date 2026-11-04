"""Load and apply the repository plotting style profiles.

This module owns visual defaults only. It does not select physical rows,
interpolate data, call solvers, or modify numerical inputs.
"""

from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any, Iterator, Mapping
import tomllib
import warnings


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROFILE_ROOT = PROJECT_ROOT / "config" / "plotting"
PROFILE_FILES = {
    "audit_v1": PROFILE_ROOT / "audit_v1.toml",
    "candidate_origin_like_v1": PROFILE_ROOT / "candidate_origin_like_v1.toml",
    "strict_origin_like_v1": PROFILE_ROOT / "strict_origin_like_v1.toml",
    "candidate_aps_v2": PROFILE_ROOT / "candidate_aps_v2.toml",
    "strict_aps_v2": PROFILE_ROOT / "strict_aps_v2.toml",
}
ACTIVE_PROFILE_IDS = frozenset({"candidate_aps_v2", "strict_aps_v2"})
LEGACY_PROFILE_IDS = frozenset(set(PROFILE_FILES) - set(ACTIVE_PROFILE_IDS))
ALLOWED_COLUMNS = {"single_column", "double_column"}
ALLOWED_PROFILES = frozenset(PROFILE_FILES)
FONT_POLICIES = frozenset({"error", "candidate_fallback", "legacy_default"})


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
    def lifecycle_status(self) -> str:
        return str(self.data.get("status", "active"))

    @property
    def is_legacy(self) -> bool:
        return self.profile_id in LEGACY_PROFILE_IDS or self.lifecycle_status != "active"

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
    status = str(data.get("status", "active"))
    if status not in {"active", "deprecated"}:
        raise ValueError(f"{path}: status must be active or deprecated")
    if profile_id in ACTIVE_PROFILE_IDS and status != "active":
        raise ValueError(f"{path}: current profile {profile_id!r} must remain active")
    if profile_id in LEGACY_PROFILE_IDS and status != "deprecated":
        raise ValueError(f"{path}: legacy profile {profile_id!r} must be deprecated")
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
    policy = data.get("font_policy", "legacy_default")
    current_path = path.resolve() == PROFILE_FILES[profile_id].resolve()
    if policy not in FONT_POLICIES or (profile_id in ACTIVE_PROFILE_IDS and policy == "legacy_default"
            and ("font_policy" in data or current_path)):
        raise ValueError(f"{path}: active profiles require font_policy=error or candidate_fallback")
    candidates = data.get("font_candidates", [])
    if policy != "legacy_default" and (not isinstance(candidates, list) or not candidates
            or any(not isinstance(item, str) or not item.strip() for item in candidates)
            or len(set(candidates)) != len(candidates)):
        raise ValueError(f"{path}: explicit font policy requires non-empty unique font_candidates")
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


def resolve_font(profile: PlotProfile) -> dict[str, Any]:
    """Resolve declared fonts; active profiles never use an undeclared default."""

    from matplotlib import font_manager

    candidates = [str(item) for item in profile.data.get("font_candidates", [])]
    requested = candidates[0] if candidates else "serif"
    policy = str(profile.data.get("font_policy", "legacy_default"))
    if policy not in FONT_POLICIES or (profile.profile_id in ACTIVE_PROFILE_IDS and policy == "legacy_default"):
        raise ValueError("active profiles require an explicit missing-font policy")
    if policy != "legacy_default" and (not candidates or len(set(candidates)) != len(candidates)
            or any(not name.strip() for name in candidates)):
        raise ValueError("explicit font policy requires non-empty unique font_candidates")
    missing = []
    path = None
    selected = requested
    if policy == "legacy_default":
        path = font_manager.findfont(font_manager.FontProperties(family=[requested]), fallback_to_default=True)
        return {"requested": requested,
                "resolved": font_manager.FontProperties(fname=path).get_name(),
                "path": str(Path(path).resolve())}
    else:
        choices = candidates[:1] if policy == "error" else candidates
        for name in choices:
            try:
                path = font_manager.findfont(font_manager.FontProperties(family=[name]), fallback_to_default=False)
                selected = name
                break
            except ValueError:
                missing.append(name)
        if path is None:
            raise ValueError(f"font_policy={policy}: required font not found: {', '.join(choices)}")
    actual = Path(path).resolve()
    resolved = font_manager.FontProperties(fname=str(actual)).get_name()
    fallback = bool(missing)
    if fallback and policy != "legacy_default":
        warnings.warn(f"{requested} is unavailable; declared font substitute is {resolved}",
                      UserWarning, stacklevel=2)
    return {"schema": "plot_font_resolution_v1", "requested": requested,
            "candidates": candidates, "policy": policy, "resolved": resolved,
            "selected": selected,
            "path": str(actual), "bytes": actual.stat().st_size,
            "sha256": hashlib.sha256(actual.read_bytes()).hexdigest(),
            "fallback_used": fallback, "missing_candidates": missing,
            "mathtext_fontset": str(profile.data.get("mathtext_fontset", "stix"))}


def _matplotlib_settings(profile: PlotProfile, font: Mapping[str, Any]) -> dict[str, Any]:
    """Defaults shared by the global CLI compatibility path and local contexts."""

    from cycler import cycler

    font_size = float(profile.data.get("font_size_pt", 9.0))
    line_width = float(profile.data.get("line_width_pt", 1.0))
    axes_linewidth = float(profile.data.get("axes_linewidth_pt", 0.6))
    settings = {
            "font.family": [font["resolved"]],
            "font.serif": [font["resolved"]],
            "font.size": font_size,
            "font.weight": "normal",
            "font.style": "normal",
            "font.variant": "normal",
            "font.stretch": "normal",
            "mathtext.fontset": profile.data.get("mathtext_fontset", "stix"),
            "mathtext.fallback": "stix",
            "mathtext.default": "it",
            "text.usetex": False,
            "text.parse_math": True,
            "text.color": "black",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "axes.edgecolor": "black",
            "axes.labelcolor": "black",
            "axes.grid": False,
            "axes.unicode_minus": True,
            "axes.formatter.use_locale": False,
            "figure.facecolor": "white",
            "figure.edgecolor": "white",
            "axes.facecolor": "white",
            "figure.autolayout": False,
            "figure.constrained_layout.use": False,
            "axes.prop_cycle": cycler(color=list(profile.colors)),
            "axes.linewidth": axes_linewidth,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size,
            "legend.fontsize": float(profile.data.get("legend_font_size_pt", max(7.0, font_size - 1.0))),
            "legend.title_fontsize": font_size,
            "legend.frameon": False,
            "xtick.color": "black",
            "ytick.color": "black",
            "xtick.labelsize": float(profile.data.get("tick_font_size_pt", max(7.0, font_size - 1.0))),
            "ytick.labelsize": float(profile.data.get("tick_font_size_pt", max(7.0, font_size - 1.0))),
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.width": max(0.5, axes_linewidth),
            "ytick.major.width": max(0.5, axes_linewidth),
            "lines.linewidth": line_width,
            "lines.linestyle": "-",
            "lines.marker": "None",
            "lines.markersize": float(profile.data.get("marker_size_pt", 4.0)),
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": profile.dpi,
            "savefig.bbox": None if profile.data["output"].get("bbox_inches") == "fixed" else "tight",
            "savefig.pad_inches": 0.08,
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
            "savefig.transparent": False,
        }
    if profile.is_legacy:
        # Keep the exact global defaults of frozen v1 renderers. The new local
        # interface is opt-in and must not restyle historical compatibility.
        previous_keys = {
            "font.family", "font.serif", "font.size", "mathtext.fontset", "axes.prop_cycle",
            "axes.linewidth", "axes.labelsize", "axes.titlesize", "legend.fontsize", "legend.frameon",
            "xtick.labelsize", "ytick.labelsize", "xtick.direction", "ytick.direction", "xtick.top",
            "ytick.right", "xtick.major.width", "ytick.major.width", "lines.linewidth", "lines.markersize",
            "pdf.fonttype", "ps.fonttype", "svg.fonttype", "savefig.dpi", "savefig.bbox", "savefig.pad_inches",
        }
        settings = {key: value for key, value in settings.items() if key in previous_keys}
        settings["font.family"] = profile.data.get("font_family", "serif")
        settings["font.serif"] = list(profile.data.get("font_candidates", []))
    return settings


def configure_matplotlib(profile: PlotProfile) -> dict[str, Any]:
    """Apply global defaults for existing standalone CLIs; prefer style_context."""
    import matplotlib

    font = resolve_font(profile)  # Fail before changing the caller's rcParams.
    matplotlib.rcParams.update(_matplotlib_settings(profile, font))
    return font


@dataclass(frozen=True)
class PlotStyleContext:
    """Resolved settings retained for creation and later measurement/export."""

    profile: PlotProfile
    font: Mapping[str, Any]
    settings: Mapping[str, Any]

    def bind(self, figure: Any) -> None:
        """Bind explicitly; figure families remain responsible for their artists."""
        existing = getattr(figure, "_plot_style_context", None)
        if existing is not None and existing is not self:
            raise ValueError("figure already has a different style context")
        figure._plot_style_context = self

    @contextmanager
    def activate(self) -> Iterator[PlotStyleContext]:
        import matplotlib

        with matplotlib.rc_context(rc=deepcopy(dict(self.settings))):
            yield self


@contextmanager
def style_context(profile: PlotProfile) -> Iterator[PlotStyleContext]:
    """Restore caller defaults on exit, including exceptional exits.

    Bind figures created inside the context. export_figure and measure_figure
    reactivate that exact context during deferred rendering and saving.
    """
    font = resolve_font(profile)
    frozen_profile = PlotProfile(profile.profile_id, deepcopy(profile.data), profile.path)
    context = PlotStyleContext(frozen_profile, font, _matplotlib_settings(frozen_profile, font))
    with context.activate():
        yield context


@contextmanager
def figure_style_context(figure: Any) -> Iterator[PlotStyleContext | None]:
    """Reactivate a bound style; unbound historical CLI figures keep their behavior."""
    context = getattr(figure, "_plot_style_context", None)
    if context is None:
        yield None
    else:
        with context.activate():
            yield context


def configure_axis_ticks(ax: Any, profile: PlotProfile) -> None:
    """Set scale-aware minor ticks after the figure family chooses its scales."""
    from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter

    settings = profile.data.get("ticks", {})
    for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
        if scale == "linear":
            axis.set_minor_locator(AutoMinorLocator(int(settings.get("linear_minor_subdivisions", 5))))
        elif scale == "log":
            axis.set_minor_locator(LogLocator(base=10, subs=settings.get("log_minor_subs", range(2, 10))))
            axis.set_minor_formatter(NullFormatter())
        else:
            raise ValueError(f"minor tick policy is undefined for {scale!r}")
    ax.tick_params(axis="both", which="both", direction="in", top=True, bottom=True, left=True, right=True)
    ax.tick_params(which="major", length=3.5, width=0.7)
    ax.tick_params(which="minor", length=2, width=0.5)
    ax.tick_params(axis="x", which="major", pad=6)
