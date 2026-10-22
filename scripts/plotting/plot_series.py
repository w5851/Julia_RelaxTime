"""Explicit curve identities, independent of selection and drawing order."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

from scripts.plotting.plot_style import PlotProfile


@dataclass(frozen=True)
class SeriesStyle:
    """One registered scientific identity and its line/marker encoding."""

    series_id: str
    color: str
    linestyle: str = "-"
    marker: str = ""

    def __post_init__(self) -> None:
        from matplotlib.colors import to_hex, to_rgba
        from matplotlib.lines import Line2D

        if not isinstance(self.series_id, str) or not self.series_id.strip() or self.series_id != self.series_id.strip():
            raise ValueError("series_id must be a non-empty string without surrounding whitespace")
        if not isinstance(self.color, str) or to_rgba(self.color)[3] != 1:
            raise ValueError(f"{self.series_id}: color must be an opaque named or hex color")
        if not isinstance(self.linestyle, str) or not isinstance(self.marker, str):
            raise ValueError(f"{self.series_id}: linestyle and marker must be strings")
        line = Line2D([], [], color=self.color, linestyle=self.linestyle, marker=self.marker)
        object.__setattr__(self, "color", to_hex(self.color))
        object.__setattr__(self, "linestyle", line.get_linestyle())

    @classmethod
    def from_profile(cls, profile: PlotProfile, series_id: str, style_index: int,
                     *, linestyle: str | None = None, marker: str = "") -> SeriesStyle:
        """Use a fixed palette slot; exhaustion requires an explicit new style."""
        if type(style_index) is not int or not 0 <= style_index < len(profile.colors):
            raise ValueError("style_index must identify an existing color slot; palette cycling is forbidden")
        lines = profile.data["palette"].get("parameter_linestyles", [])
        if linestyle is None:
            if style_index >= len(lines):
                raise ValueError("style_index has no registered linestyle; declare linestyle or a new SeriesStyle")
            linestyle = str(lines[style_index])
        return cls(series_id, profile.colors[style_index], linestyle, marker)

    def record(self) -> dict[str, str]:
        return {"series_id": self.series_id, "color": self.color,
                "linestyle": self.linestyle, "marker": self.marker}

    def plot_kwargs(self) -> dict[str, str]:
        """Include an artist identity so quality measurement checks actual lines."""
        return {"color": self.color, "linestyle": self.linestyle,
                "marker": self.marker, "gid": self.series_id}


@dataclass(frozen=True)
class SeriesStyleRegistry:
    """A complete identity table, registered before any subset is selected."""

    registry_id: str
    styles: tuple[SeriesStyle, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.registry_id, str) or not self.registry_id.strip():
            raise ValueError("registry_id must be a non-empty string")
        object.__setattr__(self, "styles", tuple(self.styles))
        if not self.styles or any(not isinstance(style, SeriesStyle) for style in self.styles):
            raise ValueError("styles must contain SeriesStyle entries")
        identities = [style.series_id for style in self.styles]
        if len(set(identities)) != len(identities):
            raise ValueError("duplicate registered series_id")
        encodings = [(style.color, style.linestyle, style.marker) for style in self.styles]
        if len(set(encodings)) != len(encodings):
            raise ValueError("registered identities must have distinct color/linestyle/marker encodings")

    def select(self, series_ids: Iterable[str]) -> tuple[SeriesStyle, ...]:
        if isinstance(series_ids, str):
            raise ValueError("series_ids must be a sequence of identities, not one string")
        selected = tuple(series_ids)
        if any(not isinstance(item, str) or not item.strip() for item in selected):
            raise ValueError("selected series_id must be a non-empty string")
        if len(set(selected)) != len(selected):
            raise ValueError("duplicate selected series_id")
        registered = {style.series_id: style for style in self.styles}
        unknown = [item for item in selected if item not in registered]
        if unknown:
            raise ValueError(f"unregistered series_id: {', '.join(unknown)}; declare its style explicitly")
        return tuple(registered[item] for item in selected)

    def record(self, series_ids: Iterable[str]) -> dict[str, Any]:
        selected = self.select(series_ids)
        return {"schema": "plot_series_styles_v1", "registry_id": self.registry_id,
                "styles": [style.record() for style in self.styles],
                "selected_ids": [style.series_id for style in selected]}


def validate_series_styles(record: Mapping[str, Any], series: list[dict[str, Any]],
                           measured: list[dict[str, Any]]) -> list[str]:
    """Compare registered identities with manifest series and measured artists."""
    try:
        if record.get("schema") != "plot_series_styles_v1":
            raise ValueError("unsupported series_styles schema")
        registry = SeriesStyleRegistry(record["registry_id"],
            tuple(SeriesStyle(**item) for item in record["styles"]))
        selected = registry.select(record["selected_ids"])
        if not selected:
            raise ValueError("series_styles.selected_ids must not be empty")
        declared = {item.series_id: item.record() for item in selected}
        if not set(declared).issubset({item.get("series_id") for item in series if isinstance(item, dict)}):
            raise ValueError("selected style identities must be present in manifest series")
        observed = set()
        for item in measured:
            identity = item["series_id"]
            if identity not in declared or item != declared[identity]:
                raise ValueError(f"measured curve style differs from registered series_id: {identity}")
            observed.add(identity)
        if observed != set(declared):
            raise ValueError("selected style identities lack measured curve evidence")
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        return [f"series_styles: {exc}"]
    return []
