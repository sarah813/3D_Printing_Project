"""Load PLA DoE measurements from the editable JSON file."""

from __future__ import annotations

import json
from pathlib import Path


VACUUM_RESPONSE_FIELDS = (
    "min_pressure_mbar",
    "pump_down_speed",
    "outgassing_rate",
)
VISUAL_DEFECT_FIELDS = (
    "blob_count",
    "has_stringing",
)


def load_pla_doe_data(path=None):
    """Read the PLA DoE JSON file."""
    data_path = Path(path) if path is not None else Path(__file__).with_name(
        "pla_doe_measurements.json"
    )
    return json.loads(data_path.read_text(encoding="utf-8"))


def levels_from_json(levels):
    """Convert JSON string level keys back to numeric coded-level keys."""
    converted = {}
    for factor, mapping in levels.items():
        if not isinstance(mapping, dict):
            converted[factor] = mapping
            continue
        try:
            converted[factor] = {
                int(float(coded_level)): value
                for coded_level, value in mapping.items()
            }
        except ValueError:
            converted[factor] = mapping
    return converted


def design_by_name(data):
    """Return all JSON designs keyed by their design name."""
    return {
        design["name"]: {
            **design,
            "levels": levels_from_json(design.get("levels", {})),
        }
        for design in data.get("designs", [])
    }


def run_labels(runs):
    """Return run labels in plotting order."""
    return [run["run"] for run in runs]


def response_defaults():
    """Create empty vacuum-response fields for a run."""
    return {field: None for field in VACUUM_RESPONSE_FIELDS}


def mass_values(run):
    """Return the direct repeated weight measurements for one run."""
    measurements = run.get("measurements", {})
    return list(measurements.get("weight_g", []))


def repeat_runs(design):
    """Collect nested repeat records from all runs in one design."""
    repeats = []
    for run in design.get("runs", []):
        for repeat in run.get("measurements", {}).get("repeats", []):
            repeat_run = {
                "run": repeat["run"],
                "mass_g": mass_values(repeat),
                "responses": {**response_defaults(), **repeat.get("responses", {})},
            }
            visual_result = repeat.get("visual_result", {})
            if visual_result:
                repeat_run["visual_defects"] = visual_result
                repeat_run.update(
                    {
                        key: value
                        for key, value in visual_result.items()
                        if value is not None
                    }
                )
            repeats.append(repeat_run)
    return repeats


def run_record(run, design):
    """Convert one JSON run to the notebook-compatible run dictionary."""
    parameters = dict(run.get("real_parameters", {}))
    parameters.update(design.get("fixed_settings", {}))
    record = {
        "run": run["run"],
        "mass_g": mass_values(run),
        "responses": {**response_defaults(), **run.get("responses", {})},
        "coded_parameters": dict(run.get("matrix_code") or {}),
        "parameters": parameters,
    }
    record.update(parameters)

    visual_result = run.get("visual_result", {})
    if visual_result:
        record["visual_defects"] = visual_result
        record.update(
            {
                key: value
                for key, value in visual_result.items()
                if value is not None
            }
        )

    if "weight_mean_g" in run.get("measurements", {}):
        record["mass_mean_g"] = run["measurements"]["weight_mean_g"]

    return record


def runs_from_design(design, include_confirmation=True):
    """Convert JSON design runs to notebook-compatible run dictionaries."""
    runs = []
    base_ironing_runs = {f"I{i}" for i in range(1, 12)}
    for run in design.get("runs", []):
        is_ironing_confirmation = (
            design.get("name") == "Ironing"
            and str(run.get("run", "")) not in base_ironing_runs
        )
        if not include_confirmation and is_ironing_confirmation:
            continue
        runs.append(run_record(run, design))
    return runs


def confirmation_runs_from_design(design):
    """Return non-matrix ironing confirmation runs."""
    base_ironing_runs = {f"I{i}" for i in range(1, 12)}
    return [
        run_record(run, design)
        for run in design.get("runs", [])
        if str(run.get("run", "")) not in base_ironing_runs
    ]


def validate_unique_runs(name, runs):
    """Fail fast when a design accidentally contains duplicate run labels."""
    labels = run_labels(runs)
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"{name} contains duplicate run labels: {duplicates}")


def mass_measurements_from_runs(runs):
    """Build the column-wise mass measurement dictionary used by the notebook."""
    max_measurements = max((len(run.get("mass_g", [])) for run in runs), default=0)
    return {
        f"mass_{measurement_no}_g": [
            run.get("mass_g", [])[measurement_no - 1]
            if len(run.get("mass_g", [])) >= measurement_no
            else None
            for run in runs
        ]
        for measurement_no in range(1, max_measurements + 1)
    }


def response_measurements_from_runs(runs, response_key):
    """Build a column-wise measurement dictionary for repeated numeric responses."""
    def values_for_run(run):
        return list(run.get("measurements", {}).get(response_key, []))

    max_measurements = max((len(values_for_run(run)) for run in runs), default=0)
    return {
        f"{response_key}_{measurement_no}": [
            values_for_run(run)[measurement_no - 1]
            if len(values_for_run(run)) >= measurement_no
            else None
            for run in runs
        ]
        for measurement_no in range(1, max_measurements + 1)
    }


def repeat_mass_measurements_from_runs(runs):
    """Build the run-label keyed repeat measurement dictionary."""
    max_measurements = max((len(run.get("mass_g", [])) for run in runs), default=0)
    return {
        f"mass_{measurement_no}_g": {
            run["run"]: (
                run.get("mass_g", [])[measurement_no - 1]
                if len(run.get("mass_g", [])) >= measurement_no
                else None
            )
            for run in runs
        }
        for measurement_no in range(1, max_measurements + 1)
    }


def visual_defect_value(run, key):
    """Read one visual-defect value from a run."""
    visual_defects = run.get("visual_defects", {})
    if isinstance(visual_defects, dict):
        return visual_defects.get(key)
    return None
