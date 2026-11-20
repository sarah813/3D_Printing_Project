"""Automatically save Matplotlib figures produced by the DoE notebooks."""

from __future__ import annotations

from collections import defaultdict
from functools import wraps
from pathlib import Path
import re
import unicodedata


SERIES_FOLDERS = (
    (
        "08_Copper_Markforged",
        ("copper markforged", "markforged copper", "markforged", "copper mf"),
    ),
    (
        "06_Copper_VirtualFoundry2",
        ("virtualfoundry2", "virtual foundry 2", "virtual foundry ii", "vf2"),
    ),
    ("05_Channel_geometry", ("channel geometry", "channel_geometry")),
    ("03_Optimization", ("box-behnken", "optimization", "extra run", "extra-run")),
    ("04_Ironing", ("ironing",)),
    ("02_FD_II", ("fd_ii", "fd ii", "rr-only", "rr-e1", "rr rectilinear")),
    ("01_FD_I", ("fd_i", "fd i", "ffd i")),
)

ANALYSIS_FOLDERS = (
    ("00_design", ("doe matrix", "design matrix", "factor-level coverage")),
    ("07_geometry", ("channel geometry", "channel width", "channel height")),
    ("04_correlations", ("correlation",)),
    ("05_vacuum", ("vacuum", "pressure", "desorption", "outgassing", "leakage")),
    ("06_quality", ("quality", "blob", "stringing", "visual")),
    ("08_combined_responses", ("response factor effects",)),
    ("02_volume", ("volume",)),
    ("03_density", ("density", "porosity")),
    ("01_mass", ("mass", "weight")),
)


def _plain_text(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    return normalized.encode("ascii", "ignore").decode("ascii")


def _slug(value: str, fallback: str = "figure") -> str:
    value = _plain_text(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value).strip("_")
    return value[:110] or fallback


def _figure_title(figure) -> str:
    if figure._suptitle is not None:
        title = figure._suptitle.get_text().strip()
        if title:
            return title

    axes_titles = [axis.get_title().strip() for axis in figure.axes if axis.get_title().strip()]
    if axes_titles:
        return " - ".join(axes_titles[:3])
    return "Untitled figure"


def _figure_signature(figure, title: str) -> tuple:
    axes_signature = tuple(
        (axis.get_title(), axis.get_xlabel(), axis.get_ylabel())
        for axis in figure.axes
    )
    return title, axes_signature


def _matching_folder(title: str, mappings, fallback: str) -> str:
    normalized = _plain_text(title).lower()
    for folder, tokens in mappings:
        if any(token in normalized for token in tokens):
            return folder
    return fallback


def configure_notebook_plot_saving(
    output_root: str | Path = "notebook_figures",
    *,
    dpi: int = 220,
    default_series: str = "00_General",
):
    """Save every figure when ``plt.show`` is called and keep inline display unchanged."""
    import matplotlib.pyplot as plt

    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    current_show = plt.show
    original_show = getattr(current_show, "_doe_original_show", current_show)
    assigned_names: dict[tuple, Path] = {}
    folder_counts = defaultdict(int)

    @wraps(original_show)
    def save_then_show(*args, **kwargs):
        figure_numbers = list(plt.get_fignums())
        try:
            for figure_number in figure_numbers:
                figure = plt.figure(figure_number)
                title = _figure_title(figure)
                signature = _figure_signature(figure, title)
                if getattr(figure, "_doe_saved_signature", None) == signature:
                    continue

                series_folder = _matching_folder(title, SERIES_FOLDERS, default_series)
                analysis_folder = _matching_folder(title, ANALYSIS_FOLDERS, "07_other")
                if analysis_folder == "07_other" and series_folder != "00_General":
                    analysis_folder = "01_mass"
                folder = output_root / series_folder / analysis_folder
                folder.mkdir(parents=True, exist_ok=True)

                if signature not in assigned_names:
                    folder_counts[folder] += 1
                    sequence = folder_counts[folder]
                    assigned_names[signature] = folder / f"{sequence:02d}_{_slug(title)}.png"

                figure.savefig(
                    assigned_names[signature],
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )
                figure._doe_saved_signature = signature

            return original_show(*args, **kwargs)
        finally:
            for figure_number in figure_numbers:
                plt.close(figure_number)

    save_then_show._doe_original_show = original_show
    plt.show = save_then_show
    return output_root.resolve()
