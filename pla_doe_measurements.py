"""PLA DoE measurement exports loaded from ``pla_doe_measurements.json``."""

from pla_doe_measurement_loader import (
    VACUUM_RESPONSE_FIELDS,
    VISUAL_DEFECT_FIELDS,
    confirmation_runs_from_design,
    design_by_name,
    load_pla_doe_data,
    mass_measurements_from_runs,
    repeat_mass_measurements_from_runs,
    repeat_runs,
    run_labels,
    runs_from_design,
    validate_unique_runs,
    visual_defect_value,
)


PLA_DOE_DATA = load_pla_doe_data()
PLA_DOE_DESIGNS = design_by_name(PLA_DOE_DATA)

DENSITY_PLA_MATTE_G_CM3 = PLA_DOE_DATA["constants"]["density_pla_matte_g_cm3"]
V_RI_MM3 = PLA_DOE_DATA["constants"]["v_ri_mm3"]
V_RII_MM3 = PLA_DOE_DATA["constants"]["v_rii_mm3"]

FD_I = PLA_DOE_DESIGNS["FD_I"]
FD_II = PLA_DOE_DESIGNS["FD_II"]
FD_II_CONCENTRIC = PLA_DOE_DESIGNS["FD_II_concentric"]
EXTRA = PLA_DOE_DESIGNS["Extra_runs"]
BOX_BEHNKEN_NEW_SCALE = PLA_DOE_DESIGNS["Box_Behnken_new_scale"]
IRONING = PLA_DOE_DESIGNS["Ironing"]

FD_I_LEVELS = FD_I["levels"]
FD_II_LEVELS = FD_II["levels"]
EXTRA_RUN_NOZZLE_LINE_WIDTH_LEVELS = EXTRA["levels"]["nozzle_line_width_grid"]
EXTRA_RUN_BOX_BEHNKEN_LEVELS = BOX_BEHNKEN_NEW_SCALE["levels"]
IRONING_FIXED_SETTINGS = IRONING["fixed_settings"]
IRONING_LEVELS = IRONING["levels"]

FD_I_RUNS = runs_from_design(FD_I)
FD_I_REPEAT_RUNS = repeat_runs(FD_I)
FD_II_RUNS = runs_from_design(FD_II)
FD_II_CONCENTRIC_RUNS = runs_from_design(FD_II_CONCENTRIC)
EXTRA_RUNS = runs_from_design(EXTRA)
BOX_BEHNKEN_NEW_SCALE_RUNS = runs_from_design(BOX_BEHNKEN_NEW_SCALE)
IRONING_RUNS = runs_from_design(IRONING, include_confirmation=False)
IRONING_REPEAT_RUNS = repeat_runs(IRONING)
IRONING_CONFIRMATION_RUNS = confirmation_runs_from_design(IRONING)

for _name, _run_records in {
    "FD_I_RUNS": FD_I_RUNS,
    "FD_I_REPEAT_RUNS": FD_I_REPEAT_RUNS,
    "FD_II_RUNS": FD_II_RUNS,
    "FD_II_CONCENTRIC_RUNS": FD_II_CONCENTRIC_RUNS,
    "EXTRA_RUNS": EXTRA_RUNS,
    "BOX_BEHNKEN_NEW_SCALE_RUNS": BOX_BEHNKEN_NEW_SCALE_RUNS,
    "IRONING_RUNS": IRONING_RUNS,
    "IRONING_REPEAT_RUNS": IRONING_REPEAT_RUNS,
    "IRONING_CONFIRMATION_RUNS": IRONING_CONFIRMATION_RUNS,
}.items():
    validate_unique_runs(_name, _run_records)


# Compatibility exports used by the existing notebook.
FD_I_MASS_1_G = [run["mass_g"][0] for run in FD_I_RUNS]
FD_I_MASS_MEASUREMENTS_G = mass_measurements_from_runs(FD_I_RUNS)
FD_I_REPEAT_MASS_1_G = {
    run["run"]: run["mass_g"][0]
    for run in FD_I_REPEAT_RUNS
}
FD_I_REPEAT_MASS_MEASUREMENTS_G = repeat_mass_measurements_from_runs(
    FD_I_REPEAT_RUNS
)

FD_II_MASS_MEASUREMENTS_G = mass_measurements_from_runs(FD_II_RUNS)

FD_II_CONCENTRIC_RUN_LABELS = run_labels(FD_II_CONCENTRIC_RUNS)
FD_II_CONCENTRIC_MASS_MEASUREMENTS_G = mass_measurements_from_runs(
    FD_II_CONCENTRIC_RUNS
)

EXTRA_RUN_LABELS = run_labels(EXTRA_RUNS)
EXTRA_RUN_MASS_MEASUREMENTS_G = mass_measurements_from_runs(EXTRA_RUNS)

BOX_BEHNKEN_NEW_SCALE_RUN_LABELS = run_labels(BOX_BEHNKEN_NEW_SCALE_RUNS)
BOX_BEHNKEN_NEW_SCALE_MASS_MEASUREMENTS_G = mass_measurements_from_runs(
    BOX_BEHNKEN_NEW_SCALE_RUNS
)

IRONING_RUN_LABELS = run_labels(IRONING_RUNS)
IRONING_MASS_MEASUREMENTS_G = mass_measurements_from_runs(IRONING_RUNS)
IRONING_REPEAT_MASS_MEASUREMENTS_G = repeat_mass_measurements_from_runs(
    IRONING_REPEAT_RUNS
)

IRONING_BLOB_COUNTS = {
    run["run"]: visual_defect_value(run, "blob_count")
    for run in IRONING_RUNS
    if visual_defect_value(run, "blob_count") is not None
}
IRONING_STRINGING_RUNS = [
    run["run"]
    for run in IRONING_RUNS
    if visual_defect_value(run, "has_stringing") is True
]
IRONING_REPEAT_BLOB_COUNTS = {
    run["run"]: visual_defect_value(run, "blob_count")
    for run in IRONING_REPEAT_RUNS
    if visual_defect_value(run, "blob_count") is not None
}
IRONING_REPEAT_STRINGING_RUNS = [
    run["run"]
    for run in IRONING_REPEAT_RUNS
    if visual_defect_value(run, "has_stringing") is True
]
