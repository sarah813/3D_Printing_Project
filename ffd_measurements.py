"""Input constants and compact experiment data for the FFD notebook.

Edit the experiment dictionaries below. Each experiment can contain:
- ``factors`` and ``levels`` for the setup
- ``runs`` for planned design runs
- ``repeats`` for repeat prints/measurements
- ``additional_runs`` for confirmation or custom runs

Run tuple format with coded settings:
    ("R1", [-1, -1, -1, -1, 1, 1], [13.681])

Run tuple format without coded settings, useful for repeats:
    ("R2.2", [13.795])

The measurement list can contain several measurements, for example
``[9.718, 9.706, None]``. Additional/custom runs may use coded values outside
``-1`` and ``1`` or specify real parameter values directly.
"""

DENSITY_PLA_MATTE_G_CM3 = 1.31

V_RI_MM3 = 80 * 30 * 4.36  # Bambu Studio
V_RII_MM3 = 64465.765  # Fusion properties

VACUUM_RESPONSE_FIELDS = (
    "min_pressure_mbar",
    "pump_down_speed",
    "outgassing_rate",
)
VISUAL_DEFECT_FIELDS = (
    "blob_count",
    "has_stringing",
)


FD_FULL_CODED = [
    [-1, -1, -1, -1, 1, 1],
    [-1, -1, -1, 1, 1, 1],
    [-1, -1, 1, -1, 1, -1],
    [-1, -1, 1, 1, 1, -1],
    [-1, 1, -1, -1, -1, 1],
    [-1, 1, -1, 1, -1, 1],
    [-1, 1, 1, -1, -1, -1],
    [-1, 1, 1, 1, -1, -1],
    [1, -1, -1, -1, -1, -1],
    [1, -1, -1, 1, -1, -1],
    [1, -1, 1, -1, -1, 1],
    [1, -1, 1, 1, -1, 1],
    [1, 1, -1, -1, 1, -1],
    [1, 1, -1, 1, 1, -1],
    [1, 1, 1, -1, 1, 1],
    [1, 1, 1, 1, 1, 1],
]
FD_NO_WALLS_CODED = [row[:-1] for row in FD_FULL_CODED]

FD_I = {
    "description": "First fractional factorial design.",
    "factors": [
        "A_flow_ratio",
        "B_nozzle_temp",
        "C_speed",
        "D_line_width",
        "E_layer_height",
        "F_walls",
    ],
    "levels": {
        "A_flow_ratio": {-1: 1.00, 1: 1.05},
        "B_nozzle_temp": {-1: 220, 1: 230},
        "C_speed": {-1: 50, 1: 200},
        "D_line_width": {-1: 0.62, 1: 0.75},
        "E_layer_height": {-1: 0.12, 1: 0.18},
        "F_walls": {-1: 4, 1: 10},
    },
    "runs": [
        ("R1", FD_FULL_CODED[0], [13.681]),
        ("R2", FD_FULL_CODED[1], [13.657]),
        ("R3", FD_FULL_CODED[2], [13.622]),
        ("R4", FD_FULL_CODED[3], [13.594]),
        ("R5", FD_FULL_CODED[4], [13.796]),
        ("R6", FD_FULL_CODED[5], [13.843]),
        ("R7", FD_FULL_CODED[6], [13.972]),
        ("R8", FD_FULL_CODED[7], [13.973]),
        ("R9", FD_FULL_CODED[8], [14.480]),
        ("R10", FD_FULL_CODED[9], [14.541]),
        ("R11", FD_FULL_CODED[10], [14.554]),
        ("R12", FD_FULL_CODED[11], [14.545]),
        ("R13", FD_FULL_CODED[12], [14.390]),
        ("R14", FD_FULL_CODED[13], [14.445]),
        ("R15", FD_FULL_CODED[14], [14.403]),
        ("R16", FD_FULL_CODED[15], [14.362]),
    ],
    "repeats": [
        ("R2.2", [13.795]),
        ("R9.2", [14.555]),
        ("R12.2", [14.628]),
    ],
}

FD_II = {
    "description": "Second fractional factorial design without wall-count factor.",
    "factors": [
        "A_flow_ratio",
        "B_nozzle_temp",
        "C_speed",
        "D_line_width",
        "E_layer_height",
    ],
    "levels": {
        "A_flow_ratio": {-1: 1.00, 1: 1.05},
        "B_nozzle_temp": {-1: 220, 1: 230},
        "C_speed": {-1: 50, 1: 100},
        "D_line_width": {-1: 0.62, 1: 0.75},
        "E_layer_height": {-1: 0.12, 1: 0.18},
    },
    "runs": [
        ("RR1", FD_NO_WALLS_CODED[0], [84.536]),
        ("RR2", FD_NO_WALLS_CODED[1], [84.550]),
        ("RR3", FD_NO_WALLS_CODED[2], [84.349]),
        ("RR4", FD_NO_WALLS_CODED[3], [84.496]),
        ("RR5", FD_NO_WALLS_CODED[4], [85.314]),
        ("RR6", FD_NO_WALLS_CODED[5], [84.483]),
        ("RR7", FD_NO_WALLS_CODED[6], [85.337]),
        ("RR8", FD_NO_WALLS_CODED[7], [84.320]),
        ("RR9", FD_NO_WALLS_CODED[8], [88.351]),
        ("RR10", FD_NO_WALLS_CODED[9], [87.760]),
        ("RR11", FD_NO_WALLS_CODED[10], [89.017]),
        ("RR12", FD_NO_WALLS_CODED[11], [88.166]),
        ("RR13", FD_NO_WALLS_CODED[12], [87.925]),
        ("RR14", FD_NO_WALLS_CODED[13], [88.548]),
        ("RR15", FD_NO_WALLS_CODED[14], [87.809]),
        ("RR16", FD_NO_WALLS_CODED[15], [88.671]),
    ],
}

FD_II_CONCENTRIC = {
    "description": "Concentric toolpath follow-up using FD II settings.",
    "factors": FD_II["factors"],
    "levels": FD_II["levels"],
    "runs": [
        ("RC1", FD_NO_WALLS_CODED[0], [83.862]),
        ("RC2", FD_NO_WALLS_CODED[1], [None]),
        ("RC3", FD_NO_WALLS_CODED[2], [None]),
        ("RC4", FD_NO_WALLS_CODED[3], [None]),
        ("RC5", FD_NO_WALLS_CODED[4], [None]),
        ("RC6", FD_NO_WALLS_CODED[5], [None]),
        ("RC7", FD_NO_WALLS_CODED[6], [None]),
        ("RC8", FD_NO_WALLS_CODED[7], [None]),
        ("RC9", FD_NO_WALLS_CODED[8], [None]),
        ("RC10", FD_NO_WALLS_CODED[9], [None]),
        ("RC11", FD_NO_WALLS_CODED[10], [87.859]),
        ("RC12", FD_NO_WALLS_CODED[11], [None]),
        ("RC13", FD_NO_WALLS_CODED[12], [None]),
        ("RC14", FD_NO_WALLS_CODED[13], [None]),
        ("RC15", FD_NO_WALLS_CODED[14], [None]),
        ("RC16", FD_NO_WALLS_CODED[15], [None]),
    ],
}

EXTRA = {
    "description": "Additional exploratory runs.",
    "levels": {
        "nozzle_line_width_grid": {
            "B_nozzle_temp": [225, 230, 235],
            "D_line_width": [0.62, 0.70, 0.80],
        },
    },
    "runs": [
        ("LD", [None, None, None, None]),
        ("HD", [None, None, None, None]),
        ("E1", [88.964, None, None, None]),
        ("E2", [None, None, None, None]),
        ("E3", [None, None, None, None]),
        ("E4", [None, None, None, None]),
        ("E5", [None, None, None, None]),
        ("E6", [None, None, None, None]),
    ],
}

BOX_BEHNKEN_CODED = [
    [0, -1, -1],
    [0, -1, 1],
    [0, 1, -1],
    [0, 1, 1],
    [-1, 0, -1],
    [-1, 0, 1],
    [1, 0, -1],
    [1, 0, 1],
    [-1, -1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [1, 1, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
]

BOX_BEHNKEN_NEW_SCALE = {
    "description": "Box-Behnken follow-up measured on the new scale.",
    "factors": [
        "B_nozzle_temp",
        "D_line_width",
        "E_layer_height",
    ],
    "fixed_settings": {
        "flow_ratio": 1.05,
        "volumetric_speed_mm3_s": 8,
    },
    "levels": {
        "B_nozzle_temp": {-1: 225, 0: 230, 1: 235},
        "D_line_width": {-1: 0.62, 0: 0.70, 1: 0.80},
        "E_layer_height": {-1: 0.10, 0: 0.12, 1: 0.14},
    },
    "runs": [
        ("B1", BOX_BEHNKEN_CODED[0], [89.090, 89.100]),
        ("B2", BOX_BEHNKEN_CODED[1], [89.015, 89.022]),
        ("B3", BOX_BEHNKEN_CODED[2], [88.986, 88.996]),
        ("B4", BOX_BEHNKEN_CODED[3], [89.171, 89.183]),
        ("B5", BOX_BEHNKEN_CODED[4], [89.313, 89.316]),
        ("B6", BOX_BEHNKEN_CODED[5], [89.789, 89.796]),
        ("B7", BOX_BEHNKEN_CODED[6], [89.720, 89.730]),
        ("B8", BOX_BEHNKEN_CODED[7], [89.306, 89.312]),
        ("B9", BOX_BEHNKEN_CODED[8], [88.868, 88.863]),
        ("B10", BOX_BEHNKEN_CODED[9], [89.093, 89.105]),
        ("B11", BOX_BEHNKEN_CODED[10], [88.492, 88.502]),
        ("B12", BOX_BEHNKEN_CODED[11], [89.101, 89.111]),
        ("B13", BOX_BEHNKEN_CODED[12], [88.821, 88.821]),
        ("B14", BOX_BEHNKEN_CODED[13], [88.855, 88.864]),
        ("B15", BOX_BEHNKEN_CODED[14], [88.890, 88.903]),
    ],
}

IRONING = {
    "description": "Ironing half-fraction screening design plus center runs.",
    "factors": [
        "A_flow_ratio",
        "B_nozzle_temp_C",
        "C_ironing_flow_percent",
        "D_ironing_speed_mm_s",
    ],
    "fixed_settings": {
        "C_ironing_line_spacing_mm": 0.15,
        "E_layer_height": 0.12,
        "F_line_width_mm": 0.62,
    },
    "levels": {
        "A_flow_ratio": {-1: 0.99, 0: 1.00, 1: 1.01},
        "B_nozzle_temp_C": {-1: 215, 0: 220, 1: 225},
        "C_ironing_flow_percent": {-1: 0, 0: 1, 1: 2},
        "D_ironing_speed_mm_s": {-1: 30, 0: 45, 1: 60},
    },
    "runs": [
        ("I1", [-1, -1, -1, -1], [9.452, None], {"blob_count": 7, "has_stringing": False}),
        ("I2", [-1, -1, 1, 1], [9.623, None], {"blob_count": 11, "has_stringing": False}),
        ("I3", [-1, 1, -1, 1], [9.477, None], {"blob_count": 2, "has_stringing": False}),
        ("I4", [-1, 1, 1, -1], [9.537, None], {"blob_count": 3, "has_stringing": False}),
        ("I5", [1, -1, -1, 1], [9.551, None], {"blob_count": 6, "has_stringing": False}),
        ("I6", [1, -1, 1, -1], [9.678, None], {"blob_count": 10, "has_stringing": True}),
        ("I7", [1, 1, -1, -1], [9.558, None], {"blob_count": 1, "has_stringing": False}),
        ("I8", [1, 1, 1, 1], [9.718, None], {"blob_count": 1, "has_stringing": False}),
        ("I9", [0, 0, 0, 0], [9.500, None], {"blob_count": 2, "has_stringing": False}),
        ("I10", [0, 0, 0, 0], [9.477, None], {"blob_count": 2, "has_stringing": False}),
        ("I11", [0, 0, 0, 0], [9.444, None], {"blob_count": 3, "has_stringing": False}),
    ],
    "repeats": [
        ("I1.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I2.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I3.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I4.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I5.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I6.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I7.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I8.2", [9.734, None], {"blob_count": 2, "has_stringing": False}),
        ("I8.3", [9.729, None], {"blob_count": None, "has_stringing": False}),
        ("I8.4", [9.690, None], {"blob_count": None, "has_stringing": False}),
        ("I9.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I10.2", [None, None], {"blob_count": None, "has_stringing": False}),
        ("I11.2", [None, None], {"blob_count": None, "has_stringing": False}),
    ],
    "additional_runs": [
        {
            "run": "IA",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 220,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 60,
            },
            "mass_mean_g": 9.650,
            "quality_score": 2.5,
            "has_stringing": False,
        },
        {
            "run": "IB",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 220,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 45,
            },
            "mass_mean_g": 9.627,
            "quality_score": 8.0,
            "has_stringing": True,
        },
        {
            "run": "I8-lowT",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 215,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 60,
            },
            "mass_mean_g": 9.670,
            "quality_score": None,
            "has_stringing": None,
        },
        {
            "run": "I8-medT",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 220,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 60,
            },
            "mass_mean_g": 9.656,
            "quality_score": None,
            "has_stringing": None,
        },
        {
            "run": "I8-T1",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 225,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 60,
            },
            "mass_mean_g": 9.753,
            "quality_score": 8,
            "has_stringing": True,
            "extra_extrusion": True,
        },
        {
            "run": "I8-T2",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 215,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 60,
            },
            "mass_mean_g": 9.735,
            "quality_score": 15,
            "has_stringing": True,
            "extra_extrusion": True,
        },
        {
            "run": "I8-T3",
            "parameters": {
                "A_flow_ratio": 1.01,
                "B_nozzle_temp_C": 225,
                "C_ironing_flow_percent": 2,
                "D_ironing_speed_mm_s": 45,
            },
            "mass_mean_g": 9.688,
            "quality_score": 8,
            "has_stringing": False,
            "extra_extrusion": True,
        },
    ],
}


def _run_labels(runs):
    return [run["run"] for run in runs]


def _validate_unique_runs(name, runs):
    labels = _run_labels(runs)
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"{name} contains duplicate run labels: {duplicates}")


def _response_defaults():
    return {field: None for field in VACUUM_RESPONSE_FIELDS}


def _split_meta(meta):
    responses = _response_defaults()
    visual_defects = {}
    other = {}

    for key, value in meta.items():
        if key in VACUUM_RESPONSE_FIELDS:
            responses[key] = value
        elif key in VISUAL_DEFECT_FIELDS:
            visual_defects[key] = value
        else:
            other[key] = value

    return responses, visual_defects, other


def _real_value_from_coded(factor, coded_value, levels):
    mapping = levels[factor]
    if coded_value in mapping:
        return mapping[coded_value]

    coded_levels = sorted(mapping)
    if len(coded_levels) < 2:
        raise ValueError(f"Cannot interpolate {factor} with fewer than two levels.")

    for lower, upper in zip(coded_levels, coded_levels[1:]):
        if lower <= coded_value <= upper:
            low_real = mapping[lower]
            high_real = mapping[upper]
            fraction = (coded_value - lower) / (upper - lower)
            return low_real + fraction * (high_real - low_real)

    raise ValueError(f"Coded value {coded_value} is outside levels for {factor}.")


def _parameters_from_coded(coded_values, experiment):
    factors = experiment.get("factors", [])
    levels = experiment.get("levels", {})
    parameters = {
        factor: _real_value_from_coded(factor, coded_value, levels)
        for factor, coded_value in zip(factors, coded_values)
    }
    parameters.update(experiment.get("fixed_settings", {}))
    return parameters


def _record_to_run(record, experiment, has_coded=True):
    if isinstance(record, dict):
        run = dict(record)
        parameters = dict(run.pop("parameters", {}))
        parameters.update(experiment.get("fixed_settings", {}))
        if parameters:
            run["parameters"] = parameters
            run.update(parameters)
        run["responses"] = {**_response_defaults(), **run.get("responses", {})}
        visual_defects = {
            key: run[key]
            for key in VISUAL_DEFECT_FIELDS
            if key in run
        }
        if visual_defects:
            run["visual_defects"] = {
                **run.get("visual_defects", {}),
                **visual_defects,
            }
        return run

    if has_coded:
        run_name, coded_values, mass_g, *rest = record
        meta = rest[0] if rest else {}
        responses, visual_defects, other = _split_meta(meta)
        run = {
            "run": run_name,
            "mass_g": mass_g,
            "responses": responses,
            "coded_parameters": dict(zip(experiment["factors"], coded_values)),
            "parameters": _parameters_from_coded(coded_values, experiment),
            **other,
        }
        run.update(run["parameters"])
    else:
        run_name, mass_g, *rest = record
        meta = rest[0] if rest else {}
        responses, visual_defects, other = _split_meta(meta)
        run = {
            "run": run_name,
            "mass_g": mass_g,
            "responses": responses,
            **other,
        }

    if visual_defects:
        run["visual_defects"] = visual_defects
        run.update(visual_defects)

    return run


def _runs(experiment, key="runs", has_coded=True):
    return [
        _record_to_run(record, experiment, has_coded=has_coded)
        for record in experiment.get(key, [])
    ]


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


def _visual_defect_value(run, key):
    visual_defects = run.get("visual_defects", {})
    if isinstance(visual_defects, dict):
        return visual_defects.get(key)
    return None


FD_I_LEVELS = FD_I["levels"]
FD_II_LEVELS = FD_II["levels"]
EXTRA_RUN_NOZZLE_LINE_WIDTH_LEVELS = EXTRA["levels"]["nozzle_line_width_grid"]
EXTRA_RUN_BOX_BEHNKEN_LEVELS = BOX_BEHNKEN_NEW_SCALE["levels"]
IRONING_FIXED_SETTINGS = IRONING["fixed_settings"]
IRONING_LEVELS = IRONING["levels"]

FD_I_RUNS = _runs(FD_I)
FD_I_REPEAT_RUNS = _runs(FD_I, "repeats", has_coded=False)
FD_II_RUNS = _runs(FD_II)
FD_II_CONCENTRIC_RUNS = _runs(FD_II_CONCENTRIC)
EXTRA_RUNS = _runs(EXTRA, has_coded=False)
BOX_BEHNKEN_NEW_SCALE_RUNS = _runs(BOX_BEHNKEN_NEW_SCALE)
IRONING_RUNS = _runs(IRONING)
IRONING_REPEAT_RUNS = _runs(IRONING, "repeats", has_coded=False)
IRONING_CONFIRMATION_RUNS = _runs(IRONING, "additional_runs", has_coded=False)

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
    _validate_unique_runs(_name, _run_records)


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
    run["run"]: _visual_defect_value(run, "blob_count")
    for run in IRONING_RUNS
    if _visual_defect_value(run, "blob_count") is not None
}
IRONING_STRINGING_RUNS = [
    run["run"]
    for run in IRONING_RUNS
    if _visual_defect_value(run, "has_stringing") is True
]
IRONING_REPEAT_BLOB_COUNTS = {
    run["run"]: _visual_defect_value(run, "blob_count")
    for run in IRONING_REPEAT_RUNS
    if _visual_defect_value(run, "blob_count") is not None
}
IRONING_REPEAT_STRINGING_RUNS = [
    run["run"]
    for run in IRONING_REPEAT_RUNS
    if _visual_defect_value(run, "has_stringing") is True
]
