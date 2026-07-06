"""Input constants and measurements for the FFD notebook."""

DENSITY_PLA_MATTE_G_CM3 = 1.31

V_RI_MM3 = 80 * 30 * 4.36  # Bambu Studio
V_RII_MM3 = 64465.765  # Fusion properties


FD_I_MASS_1_G = [
    13.681,  # R1
    13.657,  # R2
    13.622,  # R3
    13.594,  # R4
    13.796,  # R5
    13.843,  # R6
    13.972,  # R7
    13.973,  # R8
    14.480,  # R9
    14.541,  # R10
    14.554,  # R11
    14.545,  # R12
    14.390,  # R13
    14.445,  # R14
    14.403,  # R15
    14.362,  # R16
]

FD_I_MASS_MEASUREMENTS_G = {
    "mass_1_g": FD_I_MASS_1_G,
}

FD_I_REPEAT_MASS_1_G = {
    "R2.2": 13.795,
    "R9.2": 14.555,
    "R12.2": 14.628,
}

FD_I_REPEAT_MASS_MEASUREMENTS_G = {
    "mass_1_g": FD_I_REPEAT_MASS_1_G,
}

FD_II_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        84.536,  # R1
        84.550,  # R2
        84.349,  # R3
        84.496,  # R4
        85.314,  # R5
        84.483,  # R6
        85.337,  # R7
        84.320,  # R8
        88.351,  # R9
        87.760,  # R10
        89.017,  # R11
        88.166,  # R12
        87.925,  # R13
        88.548,  # R14
        87.809,  # R15
        88.671,  # R16
    ],
}

FD_II_CONCENTRIC_RUN_LABELS = [
    "RC1",
    "RC2",
    "RC3",
    "RC4",
    "RC5",
    "RC6",
    "RC7",
    "RC8",
    "RC9",
    "RC10",
    "RC11",
    "RC12",
    "RC13",
    "RC14",
    "RC15",
    "RC16",
]

FD_II_CONCENTRIC_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        83.862,  # RC1
        None,  # RC2
        None,  # RC3
        None,  # RC4
        None,  # RC5
        None,  # RC6
        None,  # RC7
        None,  # RC8
        None,  # RC9
        None,  # RC10
        87.859,  # RC11
        None,  # RC12
        None,  # RC13
        None,  # RC14
        None,  # RC15
        None,  # RC16
    ],
}


EXTRA_RUN_LABELS = [
    "LD",
    "HD",
    "E1",
    "E2",
    "E3",
    "E4",
    "E5",
    "E6",
]

EXTRA_RUN_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        35.585,  # LD
        None,  # HD
        86.642,  # E1
        None,  # E2
        None,  # E3
        None,  # E4
        None,  # E5
        None,  # E6
    ],
    "mass_2_g": [
        35.710,  # LD
        None,  # HD
        87.034,  # E1
        None,  # E2
        None,  # E3
        None,  # E4
        None,  # E5
        None,  # E6
    ],
    "mass_3_g": [
        None,  # LD
        None,  # HD
        86.975,  # E1
        None,  # E2
        None,  # E3
        None,  # E4
        None,  # E5
        None,  # E6
    ],
    "mass_4_g": [
        None,  # LD
        None,  # HD
        None,  # E1
        None,  # E2
        None,  # E3
        None,  # E4
        None,  # E5
        None,  # E6
    ],
}


BOX_BEHNKEN_NEW_SCALE_RUN_LABELS = [
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8",
    "B9",
    "B10",
    "B11",
    "B12",
    "B13",
    "B14",
    "B15",
]

BOX_BEHNKEN_NEW_SCALE_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        89.090,  # B1
        89.015,  # B2
        88.986,  # B3
        89.171,  # B4
        89.313,  # B5
        89.789,  # B6
        89.720,  # B7
        89.306,  # B8
        88.868,  # B9
        89.093,  # B10
        88.492,  # B11
        89.101,  # B12
        88.821,  # B13
        88.855,  # B14
        88.890,  # B15
    ],
    "mass_2_g": [
        89.100,  # B1
        89.022,  # B2
        88.996,  # B3
        89.183,  # B4
        89.316,  # B5
        89.796,  # B6
        89.730,  # B7
        89.312,  # B8
        88.863,  # B9
        89.105,  # B10
        88.502,  # B11
        89.111,  # B12
        88.821,  # B13
        88.864,  # B14
        88.903,  # B15
    ],
}


# Ironing screening design: planned I1-I8 half-fraction plus I9-I11 center runs.
IRONING_RUN_LABELS = [
    "I1",
    "I2",
    "I3",
    "I4",
    "I5",
    "I6",
    "I7",
    "I8",
    "I9",
    "I10",
    "I11",
]

IRONING_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        9.452,  # I1
        9.623,  # I2
        9.477,  # I3
        9.537,  # I4
        9.551,  # I5
        9.678,  # I6
        9.558,  # I7
        9.718,  # I8
        9.500,  # I9
        9.477,  # I10
        9.444,  # I11
    ],
    "mass_2_g": [
        None,  # I1
        None,  # I2
        None,  # I3
        None,  # I4
        None,  # I5
        None,  # I6
        None,  # I7
        None,  # I8
        None,  # I9
        None,  # I10
        None,  # I11
    ],
}

# Repeat prints of original ironing screening runs.
# Use the decimal suffix for the print number of the same original run.
IRONING_REPEAT_MASS_MEASUREMENTS_G = {
    "mass_1_g": {
        "I1.2": None,
        "I2.2": None,
        "I3.2": None,
        "I4.2": None,
        "I5.2": None,
        "I6.2": None,
        "I7.2": None,
        "I8.2": 9.712,
        "I8.3": 9.690,
        "I9.2": None,
        "I10.2": None,
        "I11.2": None,
    },
    "mass_2_g": {
        "I1.2": None,
        "I2.2": None,
        "I3.2": None,
        "I4.2": None,
        "I5.2": None,
        "I6.2": None,
        "I7.2": None,
        "I8.2": None,
        "I8.3": None,
        "I9.2": None,
        "I10.2": None,
        "I11.2": None,
    },
}

IRONING_BLOB_COUNTS = {
    "I1": 7,
    "I2": 11,
    "I3": 2,
    "I4": 3,
    "I5": 6,
    "I6": 10,
    "I7": 1,
    "I8": 1,
    "I9": 2,
    "I10": 2,
    "I11": 3,
}

IRONING_STRINGING_RUNS = [
    "I6",
]

IRONING_REPEAT_BLOB_COUNTS = {
    "I8.2": 2,
}

IRONING_REPEAT_STRINGING_RUNS = []

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
]
