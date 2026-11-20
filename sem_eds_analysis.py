"""Extract and summarize Thermo Fisher Axia ChemiSEM EDS PAD projects."""

from __future__ import annotations

import argparse
import io
import math
import re
import shutil
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_PAD = Path("SEM_analysis_outputs") / "T3" / "T3_2.pad"
DEFAULT_OUTPUT = Path("SEM_analysis_outputs") / "T3" / "EDS_analysis"

ELEMENT_ORDER = [
    "Cu",
    "O",
    "C",
    "N",
    "Si",
    "P",
    "S",
    "Cl",
    "Al",
    "Ca",
    "Ti",
    "Cr",
    "Nd",
    "La",
]

PEAKS_5_KV = [
    (0.277, "C K"),
    (0.392, "N K"),
    (0.525, "O K"),
    (0.929, "Cu L"),
    (1.740, "Si K"),
    (2.014, "P K"),
    (2.307, "S K"),
    (2.622, "Cl K"),
]

PEAKS_20_KV = [
    (0.277, "C K"),
    (0.392, "N K"),
    (0.525, "O K"),
    (0.929, "Cu L"),
    (1.486, "Al K"),
    (1.740, "Si K"),
    (2.014, "P K"),
    (2.307, "S K"),
    (2.622, "Cl K"),
    (3.690, "Ca K"),
    (4.511, "Ti K"),
    (8.048, "Cu K-alpha"),
    (8.905, "Cu K-beta"),
]


def _read_text(archive: zipfile.ZipFile, name: str) -> str:
    with archive.open(name) as stream:
        return stream.read().decode("utf-8", errors="replace")


def _first_float(text: str, pattern: str) -> float:
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        return math.nan
    return float(match.group(1))


def _first_int(text: str, pattern: str) -> int:
    value = _first_float(text, pattern)
    return int(value) if np.isfinite(value) else 0


def _site_ids(archive: zipfile.ZipFile) -> list[str]:
    pattern = re.compile(r"^sessions/[^/]+/sites/([^/]+)/metadata\.yml$")
    return sorted(
        match.group(1)
        for name in archive.namelist()
        if (match := pattern.match(name))
    )


def _site_prefix(archive: zipfile.ZipFile, site_id: str) -> str:
    suffix = f"/sites/{site_id}/metadata.yml"
    metadata_name = next(name for name in archive.namelist() if name.endswith(suffix))
    return metadata_name.removesuffix("metadata.yml")


def _site_name(site_metadata: str, site_id: str) -> str:
    match = re.search(r"(?m)^name:\s*['\"]?([^'\"\r\n]+)", site_metadata)
    return match.group(1).strip() if match else site_id


def _spectrum(detector_metadata: str) -> np.ndarray:
    match = re.search(
        r"(?ms)^spectrum:\s*\r?\n  data:\s*\r?\n"
        r"(?P<data>(?:    - [-+0-9.eE]+\s*\r?\n)+)",
        detector_metadata,
    )
    if not match:
        raise ValueError("The summed spectrum was not found in detector metadata.")
    return np.asarray(
        [float(value) for value in re.findall(r"(?m)^    - ([-+0-9.eE]+)", match["data"])],
        dtype=float,
    )


def _quant_channels(quant_view_metadata: str) -> list[dict[str, object]]:
    pattern = re.compile(
        r"(?ms)^  - line:\s*\r?\n"
        r"      atom: ([A-Za-z]+)\s*\r?\n"
        r"      family: ([A-Za-z]+)\s*\r?\n"
        r"    intensity_scale: ([-+0-9.eE]+)"
    )
    return [
        {"Element": atom, "Line family": family, "Preview intensity scale": float(scale)}
        for atom, family, scale in pattern.findall(quant_view_metadata)
    ]


def _element_from_map_name(name: str) -> str:
    match = re.search(r"map-([A-Za-z]+)-", Path(name).name)
    return match.group(1) if match else Path(name).stem


def _map_sort_key(name: str) -> tuple[int, str]:
    element = _element_from_map_name(name)
    try:
        return ELEMENT_ORDER.index(element), name
    except ValueError:
        return len(ELEMENT_ORDER), name


def _save_figure(figure: plt.Figure, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")


def _plot_spectrum(
    spectrum: pd.DataFrame,
    site_name: str,
    beam_energy_kv: float,
    output_stem: Path,
) -> None:
    upper_energy = min(beam_energy_kv, float(spectrum["Energy (keV)"].max()))
    visible = spectrum["Energy (keV)"] <= upper_energy
    figure, axis = plt.subplots(figsize=(12, 6.2), facecolor="white")
    axis.plot(
        spectrum.loc[visible, "Energy (keV)"],
        spectrum.loc[visible, "Counts per second"],
        color="#2f6f73",
        linewidth=1.1,
    )
    axis.set_yscale("symlog", linthresh=1.0)
    axis.set_xlim(0, upper_energy)
    axis.set_xlabel("Energy (keV)", fontsize=14)
    axis.set_ylabel("Counts per second", fontsize=14)
    axis.set_title(f"{site_name}: summed EDS spectrum ({beam_energy_kv:g} kV)", fontsize=17)
    axis.tick_params(axis="both", labelsize=12)
    axis.grid(color="#d9d9d9", linewidth=0.7, alpha=0.7)
    axis.spines[["top", "right"]].set_visible(False)

    peaks = PEAKS_5_KV if beam_energy_kv <= 5.5 else PEAKS_20_KV
    visible_y = spectrum.loc[visible, "Counts per second"].to_numpy()
    ymax = float(np.nanmax(visible_y))
    for index, (energy, label) in enumerate(peaks):
        if energy > upper_energy:
            continue
        channel = int(np.abs(spectrum["Energy (keV)"].to_numpy() - energy).argmin())
        y = max(float(spectrum.iloc[channel]["Counts per second"]), 1.0)
        axis.annotate(
            label,
            xy=(energy, y),
            xytext=(0, 18 + 13 * (index % 2)),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            rotation=90,
            arrowprops={"arrowstyle": "-", "color": "#666666", "linewidth": 0.7},
        )
    axis.set_ylim(bottom=0, top=max(ymax * 1.8, 10))
    figure.tight_layout()
    _save_figure(figure, output_stem)
    plt.close(figure)


def _save_maps_and_montage(
    archive: zipfile.ZipFile,
    prefix: str,
    site_name: str,
    output_dir: Path,
) -> None:
    map_names = sorted(
        [
            name
            for name in archive.namelist()
            if name.startswith(prefix + "views/quant/preview/map-") and name.endswith(".png")
        ],
        key=_map_sort_key,
    )
    image_name = next(
        name
        for name in archive.namelist()
        if name.startswith(prefix + "images/") and name.endswith("/image.png")
    )
    blended_name = prefix + "views/quant/preview/blended.png"
    export_dir = output_dir / "maps"
    export_dir.mkdir(parents=True, exist_ok=True)

    entries = [("SEM", image_name), ("Quant overlay", blended_name)] + [
        (_element_from_map_name(name), name) for name in map_names
    ]
    for label, name in entries:
        destination = export_dir / f"{label.replace(' ', '_')}.png"
        with archive.open(name) as source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)

    columns = 4
    rows = math.ceil(len(entries) / columns)
    figure, axes = plt.subplots(rows, columns, figsize=(16, 3.7 * rows), facecolor="white")
    axes = np.asarray(axes).reshape(-1)
    for axis, (label, name) in zip(axes, entries):
        image = plt.imread(io.BytesIO(archive.read(name)), format="png")
        axis.imshow(image, cmap="gray" if image.ndim == 2 else None)
        axis.set_title(label, fontsize=13)
        axis.axis("off")
    for axis in axes[len(entries) :]:
        axis.axis("off")
    figure.suptitle(f"{site_name}: EDS elemental maps", fontsize=20)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    _save_figure(figure, output_dir / "elemental_map_montage")
    plt.close(figure)


def analyze_pad(pad_path: Path, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    acquisition_rows: list[dict[str, object]] = []
    channel_rows: list[dict[str, object]] = []

    with zipfile.ZipFile(pad_path) as archive:
        for site_id in _site_ids(archive):
            prefix = _site_prefix(archive, site_id)
            site_metadata = _read_text(archive, prefix + "metadata.yml")
            detector_metadata = _read_text(
                archive, prefix + "eds/detectors/eds1/metadata.yml"
            )
            quant_view_metadata = _read_text(
                archive, prefix + "views/quant/metadata.yml"
            )
            name = _site_name(site_metadata, site_id)
            safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
            site_output = output_dir / safe_name
            site_output.mkdir(parents=True, exist_ok=True)

            energy_per_channel_ev = _first_float(
                detector_metadata, r"(?m)^    energy_per_channel:\s*([-+0-9.eE]+)"
            )
            beam_energy_ev = _first_float(
                detector_metadata, r"(?m)^    beam_energy:\s*([-+0-9.eE]+)"
            )
            live_time_s = _first_float(
                detector_metadata, r"(?m)^    live_time:\s*([-+0-9.eE]+)"
            )
            dead_time_percent = _first_float(
                detector_metadata, r"(?m)^    dead_time:\s*([-+0-9.eE]+)"
            )
            counts = _spectrum(detector_metadata)
            energy_kev = np.arange(counts.size) * energy_per_channel_ev / 1000.0
            spectrum = pd.DataFrame(
                {
                    "Energy (keV)": energy_kev,
                    "Counts": counts,
                    "Counts per second": counts / live_time_s,
                }
            )
            spectrum.to_csv(site_output / "summed_spectrum.csv", index=False)
            _plot_spectrum(
                spectrum,
                name,
                beam_energy_ev / 1000.0,
                site_output / "summed_spectrum",
            )
            _save_maps_and_montage(archive, prefix, name, site_output)

            acquisition_rows.append(
                {
                    "Site ID": site_id,
                    "Site name": name,
                    "Beam energy (kV)": beam_energy_ev / 1000.0,
                    "Beam current (nA)": _first_float(
                        detector_metadata,
                        r"(?m)^    beam_current:\s*([-+0-9.eE]+)",
                    )
                    * 1e9,
                    "Working distance (mm)": _first_float(
                        site_metadata, r"(?m)^  working_distance:\s*([-+0-9.eE]+)"
                    )
                    * 1000.0,
                    "Horizontal field width (um)": _first_float(
                        site_metadata, r"(?m)^  horizontal_field_width:\s*([-+0-9.eE]+)"
                    )
                    * 1e6,
                    "Live time (s)": live_time_s,
                    "Dead time (%)": dead_time_percent,
                    "Energy per channel (eV)": energy_per_channel_ev,
                    "Spectrum channels": counts.size,
                    "Input count rate (counts/s)": _first_int(
                        detector_metadata,
                        r"(?m)^  input_count_rate:\s*([-+0-9.eE]+)",
                    ),
                }
            )
            for channel in _quant_channels(quant_view_metadata):
                channel_rows.append({"Site name": name, **channel})

    acquisitions = pd.DataFrame(acquisition_rows).sort_values("Beam energy (kV)")
    channels = pd.DataFrame(channel_rows)
    acquisitions.to_csv(output_dir / "acquisition_summary.csv", index=False)
    channels.to_csv(output_dir / "quant_map_channels.csv", index=False)

    report = f"""# T3 EDS analysis

The PAD project contains two Axia ChemiSEM EDS map acquisitions at nominally the
same field of view: **T3-x1000-5kV** and **T3x1000-20kV**.

## Main interpretation

- Copper is the dominant metallic signal. The 20 kV spectrum provides the most
  defensible confirmation because it excites both Cu-L and the separated Cu-K
  lines near 8.05 and 8.90 keV.
- Carbon and oxygen are substantial low-energy signals. Carbon can include the
  specimen binder, mounting medium, carbon coating, or chamber contamination;
  oxygen is consistent with oxidized copper and/or oxygen-containing residue.
- In the maps, Cu is broadly distributed across the matrix, whereas O is
  concentrated around dark curved features and inclusions. This spatial
  anticorrelation is consistent with oxygen-rich non-matrix regions or surface
  oxidation, but the map alone cannot distinguish those mechanisms.
- The project also contains maps for N, Al, Si, P, S, Cl, Ca, and Ti. These
  should be treated as minor/trace candidates until their peaks and spatial
  distributions are checked against the spectrum and background.
- Cr, Nd, and La were selected automatically only in the 5 kV acquisition.
  Their low-energy L/M lines overlap strongly with O-K and Cu-L. They are not
  considered confirmed from this PAD alone and should be verified at higher
  voltage or with standards before reporting them as composition.
- The 20 kV acquisition has 44.4% detector dead time. It confirms the Cu-K
  peaks well, but precise quantitative comparison should use a lower count rate
  or vendor dead-time-corrected quantification.

## Quantification limitation

The PAD stores summed spectra and quantitative-map preview images, but it does
not contain a ready-made bulk weight-percent table. Preview intensity scales
are exported in `quant_map_channels.csv`; they are rendering scales, not bulk
composition values. Reliable wt%/at% should be exported from the vendor
software or recalculated with standards/background corrections.
"""
    (output_dir / "README.md").write_text(report, encoding="utf-8")
    return acquisitions, channels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pad_path", nargs="?", type=Path, default=DEFAULT_PAD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    acquisitions, channels = analyze_pad(arguments.pad_path, arguments.output_dir)
    print(acquisitions.to_string(index=False))
    print(f"\nExported {len(channels)} map channels to {arguments.output_dir.resolve()}")


if __name__ == "__main__":
    main()
