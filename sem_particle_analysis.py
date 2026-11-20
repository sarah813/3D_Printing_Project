"""Summarise ImageJ copper-particle measurements from green filaments."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ("Area", "Feret", "Circ.", "Solidity")
DEFAULT_IMAGE_WIDTH_PX = 1024
DEFAULT_IMAGE_HEIGHT_PX = 768


def load_particle_measurements(csv_path: Path | str) -> tuple[pd.DataFrame, int]:
    """Load finite, physically meaningful particle measurements."""
    csv_path = Path(csv_path)
    measurements = pd.read_csv(csv_path)
    missing_columns = [
        column for column in REQUIRED_COLUMNS if column not in measurements.columns
    ]
    if missing_columns:
        raise ValueError(
            f"{csv_path.name} is missing required columns: "
            f"{', '.join(missing_columns)}"
        )

    numeric = measurements.loc[:, REQUIRED_COLUMNS].apply(
        pd.to_numeric, errors="coerce"
    )
    finite_mask = np.isfinite(numeric).all(axis=1)
    valid_mask = (
        finite_mask
        & (numeric["Area"] > 0)
        & (numeric["Feret"] > 0)
        & numeric["Circ."].between(0, 1)
        & numeric["Solidity"].between(0, 1)
    )
    valid_measurements = numeric.loc[valid_mask].copy()
    if valid_measurements.empty:
        raise ValueError(f"{csv_path.name} contains no valid particle measurements.")

    return valid_measurements, int((~valid_mask).sum())


def create_green_filament_particle_summary(
    specimen_csvs: Mapping[str, Path | str],
    *,
    image_width_px: int = DEFAULT_IMAGE_WIDTH_PX,
    image_height_px: int = DEFAULT_IMAGE_HEIGHT_PX,
) -> pd.DataFrame:
    """Create one thesis-table row per green-filament specimen.

    D10, D50, and D90 are number-weighted percentiles of the maximum Feret
    diameter. The area fraction is the summed detected-particle area divided by
    the complete analyzed raster area.
    """
    if image_width_px <= 0 or image_height_px <= 0:
        raise ValueError("Image width and height must be positive integers.")

    image_area_px2 = int(image_width_px) * int(image_height_px)
    rows = []
    for specimen, csv_path in specimen_csvs.items():
        measurements, _ = load_particle_measurements(csv_path)
        feret_quantiles = measurements["Feret"].quantile([0.10, 0.50, 0.90])
        rows.append(
            {
                "Specimen": specimen,
                "Detected copper particles (n)": len(measurements),
                "Copper particle area fraction (%)": (
                    100.0 * measurements["Area"].sum() / image_area_px2
                ),
                "Mean particle area (px^2)": measurements["Area"].mean(),
                "D10 maximum Feret diameter (px)": feret_quantiles.loc[0.10],
                "D50 maximum Feret diameter (px)": feret_quantiles.loc[0.50],
                "D90 maximum Feret diameter (px)": feret_quantiles.loc[0.90],
                "Mean circularity (-)": measurements["Circ."].mean(),
                "Mean solidity (-)": measurements["Solidity"].mean(),
            }
        )

    return pd.DataFrame(rows)
