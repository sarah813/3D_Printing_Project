"""Create reference-calibrated EDS weight-percent estimates from PAD spectra.

The T3 composition supplied by the user is used as a transfer reference. This
is a semi-quantitative comparison, not a replacement for a vendor ZAF/Phi-Rho-Z
area quantification export.
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline
from scipy.optimize import nnls

from sem_eds_analysis import analyze_pad


OUTPUT_ROOT = Path("SEM_analysis_outputs")
REFERENCE_SPECIMEN = "T3"
SPECIMEN_DISPLAY_NAMES = {
    "DF1.25": "VF3_DF",
    "SF1.25": "VF3_SF",
}
REFERENCE_SITE = "T3x1000-20kV"
REFERENCE_WT_PERCENT = {
    "Cu": 81.36,
    "C": 14.91,
    "N": 1.32,
    "O": 1.18,
    "Si": 0.63,
    "Al": 0.31,
    "Ti": 0.08,
    "Cl": 0.07,
    "S": 0.06,
    "Ca": 0.05,
    "P": 0.02,
}

ATOMIC_NUMBER = {
    "C": 6,
    "N": 7,
    "O": 8,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ca": 20,
    "Ti": 22,
    "Fe": 26,
    "Cu": 29,
}

# Line energies in keV and approximate relative K-beta intensity. Low-Z K-beta
# lines are omitted where they are negligible at the spectrum resolution.
PEAK_LINES = {
    "C": [(0.277, 1.0)],
    "N": [(0.392, 1.0)],
    "O": [(0.525, 1.0)],
    "Al": [(1.486, 1.0), (1.557, 0.06)],
    "Si": [(1.740, 1.0), (1.836, 0.05)],
    "P": [(2.014, 1.0), (2.139, 0.05)],
    "S": [(2.307, 1.0), (2.464, 0.06)],
    "Cl": [(2.622, 1.0), (2.816, 0.08)],
    "Ca": [(3.692, 1.0), (4.013, 0.10)],
    "Ti": [(4.511, 1.0), (4.932, 0.10)],
    "Fe": [(6.404, 1.0), (7.058, 0.13)],
    "Cu": [(8.048, 1.0), (8.905, 0.14)],
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _gaussian_template(energy: np.ndarray, lines: list[tuple[float, float]]) -> np.ndarray:
    template = np.zeros_like(energy, dtype=float)
    for line_energy, relative_intensity in lines:
        energy_ev = line_energy * 1000.0
        fwhm_ev = math.sqrt(65.0**2 + 5.55 * 0.12 * 3.8 * energy_ev)
        sigma_kev = fwhm_ev / 2355.0
        template += relative_intensity * np.exp(
            -0.5 * ((energy - line_energy) / sigma_kev) ** 2
        )
    template_sum = template.sum()
    return template / template_sum if template_sum else template


def _spline_background_basis(energy: np.ndarray, spacing_kev: float = 0.80) -> np.ndarray:
    degree = 3
    interior = np.arange(energy.min() + spacing_kev, energy.max(), spacing_kev)
    knots = np.r_[
        np.repeat(energy.min(), degree + 1),
        interior,
        np.repeat(energy.max(), degree + 1),
    ]
    basis_count = len(knots) - degree - 1
    columns = []
    for index in range(basis_count):
        coefficients = np.zeros(basis_count)
        coefficients[index] = 1.0
        columns.append(BSpline(knots, coefficients, degree)(energy))
    return np.column_stack(columns)


def fit_peak_areas(
    spectrum_path: Path,
    selected_elements: list[str],
) -> tuple[pd.Series, pd.DataFrame]:
    spectrum = pd.read_csv(spectrum_path)
    fit_mask = spectrum["Energy (keV)"].between(0.16, 9.25)
    energy = spectrum.loc[fit_mask, "Energy (keV)"].to_numpy(dtype=float)
    counts_per_second = spectrum.loc[fit_mask, "Counts per second"].to_numpy(dtype=float)

    fitted_elements = [element for element in selected_elements if element in PEAK_LINES]
    background = _spline_background_basis(energy)
    peak_templates = np.column_stack(
        [_gaussian_template(energy, PEAK_LINES[element]) for element in fitted_elements]
    )
    design = np.column_stack([background, peak_templates])
    weights = 1.0 / np.sqrt(np.maximum(counts_per_second, 1.0))
    coefficients, _ = nnls(design * weights[:, None], counts_per_second * weights)
    fitted = design @ coefficients
    peak_coefficients = coefficients[background.shape[1] :]

    diagnostic = pd.DataFrame(
        {
            "Energy (keV)": energy,
            "Observed counts per second": counts_per_second,
            "Fitted counts per second": fitted,
            "Fitted background": background @ coefficients[: background.shape[1]],
        }
    )
    return pd.Series(peak_coefficients, index=fitted_elements, dtype=float), diagnostic


def _interpolated_response_factor(
    element: str,
    response_factors: pd.Series,
) -> float:
    known = response_factors.dropna()
    known_atomic_numbers = np.asarray([ATOMIC_NUMBER[item] for item in known.index])
    order = np.argsort(known_atomic_numbers)
    return float(
        np.exp(
            np.interp(
                ATOMIC_NUMBER[element],
                known_atomic_numbers[order],
                np.log(known.to_numpy()[order]),
            )
        )
    )


def estimate_composition(
    peak_areas: pd.Series,
    response_factors: pd.Series,
) -> pd.Series:
    unnormalized = {}
    for element, peak_area in peak_areas.items():
        if element in response_factors.index:
            factor = float(response_factors[element])
        else:
            factor = _interpolated_response_factor(element, response_factors)
        unnormalized[element] = max(float(peak_area) * factor, 0.0)
    estimates = pd.Series(unnormalized, dtype=float)
    return 100.0 * estimates / estimates.sum()


def _save_table_image(table: pd.DataFrame, output_path: Path, title: str) -> None:
    ordered = table.sort_values("Estimated wt%", ascending=False).reset_index(drop=True)
    figure_height = 1.1 + 0.43 * (len(ordered) + 1)
    figure, axis = plt.subplots(figsize=(4.4, figure_height), facecolor="white")
    axis.axis("off")
    axis.set_title(title, fontsize=14, fontweight="bold", pad=10)
    display_values = [
        [
            row.Element,
            (
                f"{row['Estimated wt%']:.2f}"
                if pd.notna(row["Estimated wt%"])
                else "n.d."
            ),
        ]
        for _, row in ordered.iterrows()
    ]
    display_values.append(["Total", f"{ordered['Estimated wt%'].sum():.2f}"])
    rendered = axis.table(
        cellText=display_values,
        colLabels=["Element", "Estimated wt%"],
        cellLoc="right",
        colLoc="right",
        loc="center",
        edges="horizontal",
    )
    rendered.auto_set_font_size(False)
    rendered.set_fontsize(11)
    rendered.scale(1.0, 1.35)
    for (row, column), cell in rendered.get_celld().items():
        if row == 0 or row == len(display_values):
            cell.set_text_props(fontweight="bold")
        if column == 0:
            cell.set_text_props(ha="left")
    figure.tight_layout()
    figure.savefig(output_path, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _selected_elements(channel_table: pd.DataFrame, site_name: str) -> list[str]:
    site_channels = channel_table.loc[channel_table["Site name"] == site_name, "Element"]
    return list(dict.fromkeys(site_channels.tolist()))


def _prepare_pad(pad_path: Path, specimen: str) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    output_dir = pad_path.parent / "EDS_analysis"
    acquisitions, channels = analyze_pad(pad_path, output_dir)
    acquisitions = acquisitions.sort_values("Beam energy (kV)", ascending=False)
    site_name = str(acquisitions.iloc[0]["Site name"])
    spectrum_path = output_dir / re.sub(r"[^A-Za-z0-9._-]+", "_", site_name) / "summed_spectrum.csv"
    return spectrum_path, acquisitions, channels


def main() -> None:
    specimen_folders = [path for path in OUTPUT_ROOT.iterdir() if path.is_dir()]
    specimen_folders.sort(
        key=lambda folder: (
            not any(
                pad_path.stem.lower().startswith(folder.name.lower())
                for pad_path in folder.glob("*.pad")
            ),
            folder.name,
        )
    )
    inventory_rows = []
    unique_pad_by_specimen: dict[str, Path] = {}
    seen_hashes: dict[str, tuple[str, Path]] = {}

    for folder in specimen_folders:
        pad_files = sorted(folder.glob("*.pad"))
        specimen_hashes = []
        for pad_path in pad_files:
            digest = _sha256(pad_path)
            specimen_hashes.append(digest)
            if digest not in seen_hashes:
                seen_hashes[digest] = (folder.name, pad_path)
                unique_pad_by_specimen.setdefault(folder.name, pad_path)
        duplicate_of = ""
        for digest in specimen_hashes:
            owner, owner_path = seen_hashes[digest]
            if owner != folder.name:
                duplicate_of = f"{owner}/{owner_path.name}"
                break
        inventory_rows.append(
            {
                "Specimen": SPECIMEN_DISPLAY_NAMES.get(folder.name, folder.name),
                "PAD files found": len(pad_files),
                "Unique PAD datasets assigned": int(folder.name in unique_pad_by_specimen),
                "Duplicate of": duplicate_of,
                "Status": (
                    "available"
                    if folder.name in unique_pad_by_specimen
                    else "duplicate"
                    if duplicate_of
                    else "missing"
                ),
            }
        )

    if REFERENCE_SPECIMEN not in unique_pad_by_specimen:
        raise FileNotFoundError("The T3 reference PAD is required for transfer calibration.")

    prepared = {}
    for specimen, pad_path in unique_pad_by_specimen.items():
        spectrum_path, acquisitions, channels = _prepare_pad(pad_path, specimen)
        site_name = str(acquisitions.sort_values("Beam energy (kV)", ascending=False).iloc[0]["Site name"])
        elements = _selected_elements(channels, site_name)
        peak_areas, diagnostic = fit_peak_areas(spectrum_path, elements)
        prepared[specimen] = {
            "pad_path": pad_path,
            "site_name": site_name,
            "acquisitions": acquisitions,
            "elements": elements,
            "peak_areas": peak_areas,
            "diagnostic": diagnostic,
        }

    reference_peak_areas = prepared[REFERENCE_SPECIMEN]["peak_areas"]
    reference_weights = pd.Series(REFERENCE_WT_PERCENT, dtype=float)
    common_reference = reference_weights.index.intersection(reference_peak_areas.index)
    measurable_reference = common_reference[reference_peak_areas.loc[common_reference] > 0]
    response_factors = (
        reference_weights.loc[measurable_reference]
        / reference_peak_areas.loc[measurable_reference]
    )

    combined_rows = []
    for specimen, result in prepared.items():
        display_specimen = SPECIMEN_DISPLAY_NAMES.get(specimen, specimen)
        if specimen == REFERENCE_SPECIMEN:
            estimate = reference_weights.copy()
            method = "User-supplied T3 reference"
        else:
            estimate = estimate_composition(result["peak_areas"], response_factors)
            method = "20 kV peak fit calibrated to T3 reference"
        estimate = estimate.sort_values(ascending=False)
        reported_estimate = estimate.mask(estimate <= 0.0)
        table = pd.DataFrame(
            {
                "Element": reported_estimate.index,
                "Estimated wt%": reported_estimate.values,
            }
        )
        output_dir = result["pad_path"].parent / "EDS_analysis"
        table.to_csv(output_dir / "estimated_wt_percent.csv", index=False)
        result["diagnostic"].to_csv(output_dir / "spectrum_fit_diagnostic.csv", index=False)
        _save_table_image(
            table,
            output_dir / "estimated_wt_percent.png",
            f"{display_specimen}: estimated EDS composition",
        )
        for element, value in reported_estimate.items():
            combined_rows.append(
                {
                    "Specimen": display_specimen,
                    "Acquisition": result["site_name"],
                    "Element": element,
                    "Estimated wt%": value,
                    "Quantification status": (
                        "estimated"
                        if pd.notna(value)
                        else "not reliably determined"
                    ),
                    "Method": method,
                }
            )

    combined = pd.DataFrame(combined_rows)
    combined.to_csv(OUTPUT_ROOT / "eds_estimated_wt_percent_all_specimens.csv", index=False)
    inventory = pd.DataFrame(inventory_rows)
    inventory.to_csv(OUTPUT_ROOT / "eds_pad_inventory.csv", index=False)
    print(inventory.to_string(index=False))
    print("\nEstimated compositions:")
    print(combined.pivot(index="Specimen", columns="Element", values="Estimated wt%").round(2).to_string())


if __name__ == "__main__":
    main()
