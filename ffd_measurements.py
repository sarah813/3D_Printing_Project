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
    "mass_3_g": [
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
    "mass_4_g": [
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
