"""Analyse maximum Feret diameters exported by Fiji's Analyze Particles tool."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FERET_COLUMN = "Feret"
AREA_COLUMN = "Area"
MIN_FERET_COLUMN = "MinFeret"
FERET_X_COLUMN = "FeretX"
FERET_Y_COLUMN = "FeretY"
FERET_ANGLE_COLUMN = "FeretAngle"
EQUIVALENT_DIAMETER_COLUMN = "Equivalent circular diameter"
FERET_X_SPAN_COLUMN = "Feret X span"
FERET_Y_SPAN_COLUMN = "Feret Y span"
DEFAULT_INPUT_PATH = Path("SEM_analysis_outputs") / "D1.25" / "Results.csv"
DEFAULT_OUTPUT_DIR = Path("SEM_analysis_outputs") / "D1.25"
HISTOGRAM_COLOR = "#2f6f73"
EDGE_COLOR = "#243438"
LENGTH_UNIT = "\u00b5m"
DISPLAY_MAX_PORE_LENGTH_UM = 600.0
DEFAULT_NUMBER_OF_BINS = 30


def _validate_length_scale(length_scale_to_um: float) -> float:
    """Return a valid multiplier that converts input lengths to micrometres."""
    scale = float(length_scale_to_um)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("The input-length scale must be finite and positive.")
    return scale


def load_feret_values(
    csv_path: str | Path,
    *,
    length_scale_to_um: float = 1.0,
) -> tuple[pd.Series, dict[str, int]]:
    """Read a Fiji CSV and return finite, positive Feret values in micrometres."""
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Input CSV file not found: {csv_path.resolve()}")

    measurements = pd.read_csv(csv_path)
    if FERET_COLUMN not in measurements.columns:
        available = ", ".join(str(column) for column in measurements.columns)
        raise ValueError(
            f"Required column '{FERET_COLUMN}' was not found in {csv_path.name}. "
            f"Available columns: {available}"
        )

    length_scale_to_um = _validate_length_scale(length_scale_to_um)
    numeric_feret = pd.to_numeric(measurements[FERET_COLUMN], errors="coerce")
    valid_mask = np.isfinite(numeric_feret) & (numeric_feret > 0)
    valid_values = (
        numeric_feret.loc[valid_mask].astype(float).copy() * length_scale_to_um
    )

    if valid_values.empty:
        raise ValueError(
            f"Column '{FERET_COLUMN}' contains no finite, positive measurements."
        )

    cleaning_summary = {
        "rows_read": int(len(measurements)),
        "valid_feret_values": int(valid_mask.sum()),
        "invalid_or_missing_values": int((~valid_mask).sum()),
    }
    return valid_values, cleaning_summary


def load_pore_measurements(
    csv_path: str | Path,
    *,
    length_scale_to_um: float = 1.0,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Return valid Feret/area pairs converted to micrometres and µm²."""
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Input CSV file not found: {csv_path.resolve()}")

    measurements = pd.read_csv(csv_path)
    missing_columns = [
        column
        for column in (FERET_COLUMN, AREA_COLUMN)
        if column not in measurements.columns
    ]
    if missing_columns:
        available = ", ".join(str(column) for column in measurements.columns)
        raise ValueError(
            f"Required column(s) {', '.join(missing_columns)} were not found in "
            f"{csv_path.name}. Available columns: {available}"
        )

    length_scale_to_um = _validate_length_scale(length_scale_to_um)
    optional_columns = [
        column
        for column in (
            MIN_FERET_COLUMN,
            FERET_X_COLUMN,
            FERET_Y_COLUMN,
            FERET_ANGLE_COLUMN,
        )
        if column in measurements.columns
    ]
    numeric = measurements[
        [FERET_COLUMN, AREA_COLUMN, *optional_columns]
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )
    required_numeric = numeric[[FERET_COLUMN, AREA_COLUMN]]
    valid_mask = (
        np.isfinite(required_numeric).all(axis=1)
        & (required_numeric > 0).all(axis=1)
    )
    valid_measurements = numeric.loc[valid_mask].astype(float).copy()
    valid_measurements[FERET_COLUMN] *= length_scale_to_um
    valid_measurements[AREA_COLUMN] *= length_scale_to_um**2
    if MIN_FERET_COLUMN in valid_measurements:
        valid_measurements[MIN_FERET_COLUMN] *= length_scale_to_um

    valid_measurements[EQUIVALENT_DIAMETER_COLUMN] = 2.0 * np.sqrt(
        valid_measurements[AREA_COLUMN] / np.pi
    )
    if FERET_ANGLE_COLUMN in valid_measurements:
        angle_radians = np.deg2rad(valid_measurements[FERET_ANGLE_COLUMN])
        valid_measurements[FERET_X_SPAN_COLUMN] = np.abs(
            valid_measurements[FERET_COLUMN] * np.cos(angle_radians)
        )
        valid_measurements[FERET_Y_SPAN_COLUMN] = np.abs(
            valid_measurements[FERET_COLUMN] * np.sin(angle_radians)
        )

    if valid_measurements.empty:
        raise ValueError(
            f"Columns '{FERET_COLUMN}' and '{AREA_COLUMN}' contain no rows with "
            "finite, positive paired measurements."
        )

    cleaning_summary = {
        "rows_read": int(len(measurements)),
        "valid_pore_measurements": int(valid_mask.sum()),
        "invalid_or_missing_pairs": int((~valid_mask).sum()),
    }
    return valid_measurements, cleaning_summary


def calculate_pore_geometry_statistics(
    pore_measurements: pd.DataFrame,
) -> pd.DataFrame:
    """Summarise pore area, caliper dimensions, and Feret locations."""
    rows = [
        ("Mean projected pore area", pore_measurements[AREA_COLUMN].mean(), "µm²"),
        ("Median projected pore area", pore_measurements[AREA_COLUMN].median(), "µm²"),
        (
            "Equivalent circular diameter of mean pore area",
            2.0 * np.sqrt(pore_measurements[AREA_COLUMN].mean() / np.pi),
            LENGTH_UNIT,
        ),
        (
            "Mean equivalent circular diameter",
            pore_measurements[EQUIVALENT_DIAMETER_COLUMN].mean(),
            LENGTH_UNIT,
        ),
        (
            "Median equivalent circular diameter",
            pore_measurements[EQUIVALENT_DIAMETER_COLUMN].median(),
            LENGTH_UNIT,
        ),
    ]

    optional_statistics = (
        (MIN_FERET_COLUMN, "Mean minimum Feret diameter", LENGTH_UNIT, "mean"),
        (MIN_FERET_COLUMN, "Median minimum Feret diameter", LENGTH_UNIT, "median"),
        (FERET_X_SPAN_COLUMN, "Mean Feret X span", LENGTH_UNIT, "mean"),
        (FERET_Y_SPAN_COLUMN, "Mean Feret Y span", LENGTH_UNIT, "mean"),
        (FERET_X_COLUMN, "Mean Feret X coordinate", "px", "mean"),
        (FERET_Y_COLUMN, "Mean Feret Y coordinate", "px", "mean"),
    )
    for column, label, unit, aggregation in optional_statistics:
        if column in pore_measurements:
            rows.append(
                (label, getattr(pore_measurements[column], aggregation)(), unit)
            )

    return pd.DataFrame(rows, columns=["Statistic", "Value", "Unit"])


def calculate_descriptive_statistics(feret_values: pd.Series) -> pd.DataFrame:
    """Calculate the requested pore-length summary without outlier removal."""
    values = feret_values.to_numpy(dtype=float)
    statistics = [
        ("Number of pores", float(values.size), "count"),
        ("Mean maximum Feret diameter", float(np.mean(values)), LENGTH_UNIT),
        ("Median maximum Feret diameter", float(np.median(values)), LENGTH_UNIT),
        (
            "Standard deviation of maximum Feret diameter",
            float(np.std(values, ddof=1)),
            LENGTH_UNIT,
        ),
        ("Smallest maximum Feret diameter", float(np.min(values)), LENGTH_UNIT),
        ("Largest maximum Feret diameter", float(np.max(values)), LENGTH_UNIT),
        (
            "25th percentile of maximum Feret diameter",
            float(np.percentile(values, 25)),
            LENGTH_UNIT,
        ),
        (
            "75th percentile of maximum Feret diameter",
            float(np.percentile(values, 75)),
            LENGTH_UNIT,
        ),
    ]
    return pd.DataFrame(statistics, columns=["Statistic", "Value", "Unit"])


def _style_histogram_axis(axis: plt.Axes) -> None:
    axis.set_xlabel("Maximum Feret diameter (\u00b5m)", fontsize=12)
    axis.set_ylabel("Number of pores", fontsize=12)
    axis.tick_params(axis="both", labelsize=10)
    axis.set_axisbelow(True)
    axis.grid(axis="y", color="#d9d9d9", linewidth=0.7, alpha=0.75)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def _add_statistics_and_display_limit(
    axis: plt.Axes,
    values: np.ndarray,
    *,
    display_max_pore_length_um: float | None = DISPLAY_MAX_PORE_LENGTH_UM,
    logarithmic_axis: bool = False,
    bin_description: str,
) -> None:
    """Add the median and document values beyond the visible range."""
    median_feret = float(np.median(values))

    axis.axvline(
        median_feret,
        color="#b1822b",
        linestyle=":",
        linewidth=2.0,
        label=f"Median: {median_feret:.2f} µm",
        zorder=4,
    )
    axis.lines[-1].set_label(f"Median: {median_feret:.2f} \u00b5m")

    if display_max_pore_length_um is None:
        if logarithmic_axis:
            axis.set_xlim(float(np.min(values)), float(np.max(values)))
        else:
            axis.set_xlim(0, float(np.max(values)))
    else:
        if logarithmic_axis:
            axis.set_xlim(float(np.min(values)), display_max_pore_length_um)
        else:
            axis.set_xlim(0, display_max_pore_length_um)

    axis.legend(loc="upper right", frameon=True, fontsize=9)

    hidden_count = (
        int(np.sum(values > display_max_pore_length_um))
        if display_max_pore_length_um is not None
        else 0
    )
    if hidden_count and display_max_pore_length_um is not None:
        pore_word = "pore" if hidden_count == 1 else "pores"
        axis.text(
            0.99,
            0.88,
            f"{hidden_count} {pore_word} above "
            f"{display_max_pore_length_um:g} \u00b5m not displayed",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            color="#555555",
        )


def create_linear_histogram(
    feret_values: pd.Series,
    sample_name: str = "D1.25",
    display_max_pore_length_um: float = DISPLAY_MAX_PORE_LENGTH_UM,
    number_of_bins: int = DEFAULT_NUMBER_OF_BINS,
) -> plt.Figure:
    """Create equally spaced bins across the visible pore-length range."""
    figure, axis = plt.subplots(figsize=(7.2, 5.0), facecolor="white")
    values = feret_values.to_numpy(dtype=float)
    linear_bins = np.linspace(0, display_max_pore_length_um, number_of_bins + 1)
    axis.hist(
        values,
        bins=linear_bins,
        color=HISTOGRAM_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.55,
        alpha=0.92,
    )
    bin_width = display_max_pore_length_um / number_of_bins
    axis.set_title("(a) Number-weighted pore size distribution", fontsize=14, pad=12)
    _style_histogram_axis(axis)
    _add_statistics_and_display_limit(
        axis,
        values,
        display_max_pore_length_um=display_max_pore_length_um,
        bin_description=f"Bin width: {bin_width:g} \u00b5m",
    )
    figure.tight_layout()
    return figure


def create_logarithmic_histogram(
    feret_values: pd.Series,
    sample_name: str = "D1.25",
    display_max_pore_length_um: float = DISPLAY_MAX_PORE_LENGTH_UM,
    number_of_bins: int = DEFAULT_NUMBER_OF_BINS,
) -> plt.Figure:
    """Create logarithmic bins across the complete measured diameter range."""
    values = feret_values.to_numpy(dtype=float)
    minimum = float(np.min(values))
    maximum = float(np.max(values))

    logarithmic_bins = np.geomspace(
        minimum,
        maximum,
        number_of_bins + 1,
    )
    figure, axis = plt.subplots(figsize=(7.2, 5.0), facecolor="white")
    axis.hist(
        values,
        bins=logarithmic_bins,
        color=HISTOGRAM_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.55,
        alpha=0.92,
    )
    axis.set_xscale("log")
    axis.set_title(
        f"{sample_name}: Pore size distribution (logarithmic x-axis)",
        fontsize=14,
        pad=42,
    )
    _style_histogram_axis(axis)
    _add_statistics_and_display_limit(
        axis,
        values,
        display_max_pore_length_um=None,
        logarithmic_axis=True,
        bin_description=f"{number_of_bins} logarithmic bins",
    )
    figure.tight_layout()
    return figure


def calculate_area_weighted_distribution(
    pore_measurements: pd.DataFrame,
    bin_width_um: float = 5.0,
) -> pd.DataFrame:
    """Sum projected pore area in fixed-width maximum-Feret-diameter bins."""
    if bin_width_um <= 0:
        raise ValueError("The area-weighted histogram bin width must be positive.")

    feret = pore_measurements[FERET_COLUMN].to_numpy(dtype=float)
    areas = pore_measurements[AREA_COLUMN].to_numpy(dtype=float)
    upper_edge = np.ceil(np.max(feret) / bin_width_um) * bin_width_um
    bin_edges = np.arange(0, upper_edge + bin_width_um, bin_width_um)
    area_by_bin, _ = np.histogram(feret, bins=bin_edges, weights=areas)
    count_by_bin, _ = np.histogram(feret, bins=bin_edges)
    total_area = float(np.sum(areas))

    lower_edges = bin_edges[:-1]
    upper_edges = bin_edges[1:]
    return pd.DataFrame(
        {
            "Lower maximum Feret diameter (um)": lower_edges,
            "Upper maximum Feret diameter (um)": upper_edges,
            "Bin label": [
                f"{lower:g}-{upper:g}"
                for lower, upper in zip(lower_edges, upper_edges)
            ],
            "Number of pores": count_by_bin.astype(int),
            "Projected pore area (um^2)": area_by_bin,
            "Share of total projected pore area (%)": area_by_bin / total_area * 100,
        }
    )


def create_area_weighted_histogram(
    pore_measurements: pd.DataFrame,
    sample_name: str = "D1.25",
    display_max_pore_length_um: float = DISPLAY_MAX_PORE_LENGTH_UM,
    bin_width_um: float = 5.0,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Create a fixed-bin distribution weighted by projected two-dimensional area."""
    distribution = calculate_area_weighted_distribution(
        pore_measurements,
        bin_width_um,
    )
    feret = pore_measurements[FERET_COLUMN].to_numpy(dtype=float)
    areas = pore_measurements[AREA_COLUMN].to_numpy(dtype=float)
    plot_max_um = np.ceil(np.max(feret) / bin_width_um) * bin_width_um

    figure, axis = plt.subplots(figsize=(7.2, 5.0), facecolor="white")
    axis.bar(
        distribution["Lower maximum Feret diameter (um)"],
        distribution["Share of total projected pore area (%)"],
        width=bin_width_um,
        align="edge",
        color=HISTOGRAM_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.55,
        alpha=0.92,
    )
    axis.set_xlim(0, plot_max_um)
    axis.set_title(
        "(b) Area-weighted pore size distribution",
        fontsize=14,
        pad=12,
    )
    axis.set_xlabel("Maximum Feret diameter (\u00b5m)", fontsize=12)
    axis.set_ylabel("Share of total projected pore area (%)", fontsize=12)
    axis.locator_params(axis="x", nbins=6)
    axis.tick_params(axis="both", labelsize=10)
    axis.set_axisbelow(True)
    axis.grid(axis="y", color="#d9d9d9", linewidth=0.7, alpha=0.75)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.tight_layout()
    return figure, distribution


def save_figure(figure: plt.Figure, output_stem: str | Path) -> None:
    """Save a figure as a 300 dpi PNG and a vector PDF."""
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output_stem.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    figure.savefig(
        output_stem.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
    )


def run_analysis(
    csv_path: str | Path = DEFAULT_INPUT_PATH,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    *,
    sample_name: str = "D1.25",
    length_scale_to_um: float = 1.0,
    display_max_pore_length_um: float = DISPLAY_MAX_PORE_LENGTH_UM,
    number_of_bins: int = DEFAULT_NUMBER_OF_BINS,
    close_figures: bool = True,
) -> tuple[
    pd.Series,
    pd.DataFrame,
    dict[str, int],
    tuple[plt.Figure, plt.Figure, plt.Figure],
]:
    """Run the complete analysis and write all requested output files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    feret_values, cleaning_summary = load_feret_values(
        csv_path,
        length_scale_to_um=length_scale_to_um,
    )
    statistics = calculate_descriptive_statistics(feret_values)
    statistics.to_csv(output_dir / "pore_size_descriptive_statistics.csv", index=False)

    linear_figure = create_linear_histogram(
        feret_values,
        sample_name,
        display_max_pore_length_um,
        number_of_bins,
    )
    save_figure(linear_figure, output_dir / "pore_size_distribution_linear_bins")

    logarithmic_figure = create_logarithmic_histogram(
        feret_values,
        sample_name,
        display_max_pore_length_um,
        number_of_bins,
    )
    save_figure(
        logarithmic_figure,
        output_dir / "pore_size_distribution_logarithmic_bins",
    )

    pore_measurements, _ = load_pore_measurements(
        csv_path,
        length_scale_to_um=length_scale_to_um,
    )
    pore_measurements.to_csv(
        output_dir / "pore_geometry_measurements.csv",
        index=False,
    )
    geometry_statistics = calculate_pore_geometry_statistics(pore_measurements)
    geometry_statistics.to_csv(
        output_dir / "pore_geometry_descriptive_statistics.csv",
        index=False,
    )
    bin_width_um = display_max_pore_length_um / number_of_bins
    area_weighted_figure, area_weighted_distribution = (
        create_area_weighted_histogram(
            pore_measurements,
            sample_name,
            display_max_pore_length_um,
            bin_width_um,
        )
    )
    area_weighted_distribution.to_csv(
        output_dir / "pore_area_weighted_distribution.csv",
        index=False,
    )
    save_figure(
        area_weighted_figure,
        output_dir / "pore_area_weighted_distribution",
    )

    if close_figures:
        plt.close(linear_figure)
        plt.close(logarithmic_figure)
        plt.close(area_weighted_figure)

    return (
        feret_values,
        statistics,
        cleaning_summary,
        (linear_figure, logarithmic_figure, area_weighted_figure),
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyse maximum Feret diameters measured in an SEM image."
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help=f"Fiji CSV file (default: {DEFAULT_INPUT_PATH})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--sample-name",
        default="D1.25",
        help="Sample identifier shown in plot titles (default: D1.25)",
    )
    parser.add_argument(
        "--input-unit",
        choices=("um", "mm"),
        default="um",
        help="Length unit used by the Fiji CSV (default: um)",
    )
    parser.add_argument(
        "--display-max",
        type=float,
        default=DISPLAY_MAX_PORE_LENGTH_UM,
        help="Maximum Feret diameter shown in the linear plot (default: 600)",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=DEFAULT_NUMBER_OF_BINS,
        help="Number of histogram bins (default: 30)",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    try:
        _, statistics, cleaning_summary, _ = run_analysis(
            arguments.csv_path,
            arguments.output_dir,
            sample_name=arguments.sample_name,
            length_scale_to_um=1000.0 if arguments.input_unit == "mm" else 1.0,
            display_max_pore_length_um=arguments.display_max,
            number_of_bins=arguments.bins,
        )
    except (FileNotFoundError, ValueError, pd.errors.ParserError) as error:
        raise SystemExit(f"Error: {error}") from error

    print(
        f"Analysed {cleaning_summary['valid_feret_values']} pores; "
        f"removed {cleaning_summary['invalid_or_missing_values']} invalid or missing entries."
    )
    print(statistics.to_string(index=False))
    print(f"Outputs saved to: {arguments.output_dir.resolve()}")


if __name__ == "__main__":
    main()
