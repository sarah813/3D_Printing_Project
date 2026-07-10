"""Input constants and measurements for the FFD notebook.

Edit the ``*_RUNS`` lists below. The older notebook constants are generated
from these records at the bottom of the file so the notebook can keep using
the existing imports while the data stays easier to edit.
"""

DENSITY_PLA_MATTE_G_CM3 = 1.31

V_RI_MM3 = 80 * 30 * 4.36  # Bambu Studio
V_RII_MM3 = 64465.765  # Fusion properties


def _run_labels(runs):
    return [run["run"] for run in runs]


def _validate_unique_runs(name, runs):
    labels = _run_labels(runs)
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"{name} contains duplicate run labels: {duplicates}")


def _mass_measurements_from_runs(runs):
    max_measurements = max(len(run.get("mass_g", [])) for run in runs)
    return {
        f"mass_{measurement_no}_g": [
            run.get("mass_g", [])[measurement_no - 1]
            if len(run.get("mass_g", [])) >= measurement_no
            else None
            for run in runs
        ]
        for measurement_no in range(1, max_measurements + 1)
    }


def _repeat_mass_measurements_from_runs(runs):
    max_measurements = max(len(run.get("mass_g", [])) for run in runs)
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


FD_I_RUNS = [
    {"run": "R1", "mass_g": [13.681]},
    {"run": "R2", "mass_g": [13.657]},
    {"run": "R3", "mass_g": [13.622]},
    {"run": "R4", "mass_g": [13.594]},
    {"run": "R5", "mass_g": [13.796]},
    {"run": "R6", "mass_g": [13.843]},
    {"run": "R7", "mass_g": [13.972]},
    {"run": "R8", "mass_g": [13.973]},
    {"run": "R9", "mass_g": [14.480]},
    {"run": "R10", "mass_g": [14.541]},
    {"run": "R11", "mass_g": [14.554]},
    {"run": "R12", "mass_g": [14.545]},
    {"run": "R13", "mass_g": [14.390]},
    {"run": "R14", "mass_g": [14.445]},
    {"run": "R15", "mass_g": [14.403]},
    {"run": "R16", "mass_g": [14.362]},
]

FD_I_REPEAT_RUNS = [
    {"run": "R2.2", "mass_g": [13.795]},
    {"run": "R9.2", "mass_g": [14.555]},
    {"run": "R12.2", "mass_g": [14.628]},
]

FD_II_RUNS = [
    {"run": "RR1", "mass_g": [84.536]},
    {"run": "RR2", "mass_g": [84.550]},
    {"run": "RR3", "mass_g": [84.349]},
    {"run": "RR4", "mass_g": [84.496]},
    {"run": "RR5", "mass_g": [85.314]},
    {"run": "RR6", "mass_g": [84.483]},
    {"run": "RR7", "mass_g": [85.337]},
    {"run": "RR8", "mass_g": [84.320]},
    {"run": "RR9", "mass_g": [88.351]},
    {"run": "RR10", "mass_g": [87.760]},
    {"run": "RR11", "mass_g": [89.017]},
    {"run": "RR12", "mass_g": [88.166]},
    {"run": "RR13", "mass_g": [87.925]},
    {"run": "RR14", "mass_g": [88.548]},
    {"run": "RR15", "mass_g": [87.809]},
    {"run": "RR16", "mass_g": [88.671]},
]

FD_II_CONCENTRIC_RUNS = [
    {"run": "RC1", "mass_g": [83.862]},
    {"run": "RC2", "mass_g": [None]},
    {"run": "RC3", "mass_g": [None]},
    {"run": "RC4", "mass_g": [None]},
    {"run": "RC5", "mass_g": [None]},
    {"run": "RC6", "mass_g": [None]},
    {"run": "RC7", "mass_g": [None]},
    {"run": "RC8", "mass_g": [None]},
    {"run": "RC9", "mass_g": [None]},
    {"run": "RC10", "mass_g": [None]},
    {"run": "RC11", "mass_g": [87.859]},
    {"run": "RC12", "mass_g": [None]},
    {"run": "RC13", "mass_g": [None]},
    {"run": "RC14", "mass_g": [None]},
    {"run": "RC15", "mass_g": [None]},
    {"run": "RC16", "mass_g": [None]},
]

EXTRA_RUNS = [
    {"run": "LD", "mass_g": [35.585, 35.710, None, None]},
    {"run": "HD", "mass_g": [None, None, None, None]},
    {"run": "E1", "mass_g": [86.642, 87.034, 86.975, None]},
    {"run": "E2", "mass_g": [None, None, None, None]},
    {"run": "E3", "mass_g": [None, None, None, None]},
    {"run": "E4", "mass_g": [None, None, None, None]},
    {"run": "E5", "mass_g": [None, None, None, None]},
    {"run": "E6", "mass_g": [None, None, None, None]},
]

BOX_BEHNKEN_NEW_SCALE_RUNS = [
    {"run": "B1", "mass_g": [89.090, 89.100]},
    {"run": "B2", "mass_g": [89.015, 89.022]},
    {"run": "B3", "mass_g": [88.986, 88.996]},
    {"run": "B4", "mass_g": [89.171, 89.183]},
    {"run": "B5", "mass_g": [89.313, 89.316]},
    {"run": "B6", "mass_g": [89.789, 89.796]},
    {"run": "B7", "mass_g": [89.720, 89.730]},
    {"run": "B8", "mass_g": [89.306, 89.312]},
    {"run": "B9", "mass_g": [88.868, 88.863]},
    {"run": "B10", "mass_g": [89.093, 89.105]},
    {"run": "B11", "mass_g": [88.492, 88.502]},
    {"run": "B12", "mass_g": [89.101, 89.111]},
    {"run": "B13", "mass_g": [88.821, 88.821]},
    {"run": "B14", "mass_g": [88.855, 88.864]},
    {"run": "B15", "mass_g": [88.890, 88.903]},
]

# Ironing screening design: planned I1-I8 half-fraction plus I9-I11 center runs.
IRONING_RUNS = [
    {"run": "I1", "mass_g": [9.452, None], "blob_count": 7, "has_stringing": False},
    {"run": "I2", "mass_g": [9.623, None], "blob_count": 11, "has_stringing": False},
    {"run": "I3", "mass_g": [9.477, None], "blob_count": 2, "has_stringing": False},
    {"run": "I4", "mass_g": [9.537, None], "blob_count": 3, "has_stringing": False},
    {"run": "I5", "mass_g": [9.551, None], "blob_count": 6, "has_stringing": False},
    {"run": "I6", "mass_g": [9.678, None], "blob_count": 10, "has_stringing": True},
    {"run": "I7", "mass_g": [9.558, None], "blob_count": 1, "has_stringing": False},
    {"run": "I8", "mass_g": [9.718, None], "blob_count": 1, "has_stringing": False},
    {"run": "I9", "mass_g": [9.500, None], "blob_count": 2, "has_stringing": False},
    {"run": "I10", "mass_g": [9.477, None], "blob_count": 2, "has_stringing": False},
    {"run": "I11", "mass_g": [9.444, None], "blob_count": 3, "has_stringing": False},
]

# Repeat prints of original ironing screening runs.
# Use the decimal suffix for the print number of the same original run.
IRONING_REPEAT_RUNS = [
    {"run": "I1.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I2.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I3.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I4.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I5.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I6.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I7.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I8.2", "mass_g": [9.734, None], "blob_count": 2, "has_stringing": False},
    {"run": "I8.3", "mass_g": [9.729, None], "blob_count": None, "has_stringing": False},
    {"run": "I8.4", "mass_g": [9.690, None], "blob_count": None, "has_stringing": False},
    {"run": "I9.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I10.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
    {"run": "I11.2", "mass_g": [None, None], "blob_count": None, "has_stringing": False},
]

# Ironing confirmation/optimization runs live here, separate from the planned
# screening matrix. Add future mixed-level checks or optimized settings below.
IRONING_CONFIRMATION_RUNS = [
    {
        "run": "IA",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 220,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 60,
        "mass_mean_g": 9.650,
        "quality_score": 2.5,
        "has_stringing": False,
    },
    {
        "run": "IB",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 220,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 45,
        "mass_mean_g": 9.627,
        "quality_score": 8.0,
        "has_stringing": True,
    },
    {
        "run": "I8-lowT",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 215,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 60,
        "mass_mean_g": 9.670,
        "quality_score": None,
        "has_stringing": None,
    },
    {
        "run": "I8-medT",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 220,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 60,
        "mass_mean_g": 9.656,
        "quality_score": None,
        "has_stringing": None,
    },
    {
        "run": "I8-T1",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 225,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 60,
        "mass_mean_g": 9.753,
        "quality_score": 8,
        "has_stringing": True,
        "extra_extrusion": True,
    },
    {
        "run": "I8-T2",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 215,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 60,
        "mass_mean_g": 9.735,
        "quality_score": 15,
        "has_stringing": True,
        "extra_extrusion": True,
    },
    {
        "run": "I8-T3",
        "A_flow_ratio": 1.01,
        "B_nozzle_temp_C": 225,
        "C_ironing_flow_percent": 2,
        "D_ironing_speed_mm_s": 45,
        "mass_mean_g": 9.688,
        "quality_score": 8,
        "has_stringing": False,
        "extra_extrusion": True,
    },
]


for _name, _runs in {
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
    _validate_unique_runs(_name, _runs)


# Compatibility exports used by the existing notebook.
FD_I_MASS_1_G = [run["mass_g"][0] for run in FD_I_RUNS]
FD_I_MASS_MEASUREMENTS_G = _mass_measurements_from_runs(FD_I_RUNS)
FD_I_REPEAT_MASS_1_G = {
    run["run"]: run["mass_g"][0]
    for run in FD_I_REPEAT_RUNS
}
FD_I_REPEAT_MASS_MEASUREMENTS_G = _repeat_mass_measurements_from_runs(
    FD_I_REPEAT_RUNS
)

FD_II_MASS_MEASUREMENTS_G = _mass_measurements_from_runs(FD_II_RUNS)

FD_II_CONCENTRIC_RUN_LABELS = _run_labels(FD_II_CONCENTRIC_RUNS)
FD_II_CONCENTRIC_MASS_MEASUREMENTS_G = _mass_measurements_from_runs(
    FD_II_CONCENTRIC_RUNS
)

EXTRA_RUN_LABELS = _run_labels(EXTRA_RUNS)
EXTRA_RUN_MASS_MEASUREMENTS_G = _mass_measurements_from_runs(EXTRA_RUNS)

BOX_BEHNKEN_NEW_SCALE_RUN_LABELS = _run_labels(BOX_BEHNKEN_NEW_SCALE_RUNS)
BOX_BEHNKEN_NEW_SCALE_MASS_MEASUREMENTS_G = _mass_measurements_from_runs(
    BOX_BEHNKEN_NEW_SCALE_RUNS
)

IRONING_RUN_LABELS = _run_labels(IRONING_RUNS)
IRONING_MASS_MEASUREMENTS_G = _mass_measurements_from_runs(IRONING_RUNS)
IRONING_REPEAT_MASS_MEASUREMENTS_G = _repeat_mass_measurements_from_runs(
    IRONING_REPEAT_RUNS
)

IRONING_BLOB_COUNTS = {
    run["run"]: run["blob_count"]
    for run in IRONING_RUNS
    if run.get("blob_count") is not None
}
IRONING_STRINGING_RUNS = [
    run["run"]
    for run in IRONING_RUNS
    if run.get("has_stringing") is True
]
IRONING_REPEAT_BLOB_COUNTS = {
    run["run"]: run["blob_count"]
    for run in IRONING_REPEAT_RUNS
    if run.get("blob_count") is not None
}
IRONING_REPEAT_STRINGING_RUNS = [
    run["run"]
    for run in IRONING_REPEAT_RUNS
    if run.get("has_stringing") is True
]
