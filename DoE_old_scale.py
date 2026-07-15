"""Archived old-scale Box-Behnken mass measurements.

These values are kept separate from the current DoE analysis because the
notebook now uses only the new-scale measurements.
"""

FD_I_OLD_SCALE_MASS_1_G = [
    12.907,  # R1
    13.334,  # R2
    13.294,  # R3
    13.259,  # R4
    13.454,  # R5
    13.490,  # R6
    13.622,  # R7
    13.637,  # R8
    14.136,  # R9
    14.188,  # R10
    14.206,  # R11
    14.200,  # R12
    14.049,  # R13
    14.101,  # R14
    14.061,  # R15
    14.025,  # R16
]

FD_I_OLD_SCALE_MASS_2_G = [
    13.272,  # R1
    13.296,  # R2
    13.286,  # R3
    13.240,  # R4
    13.432,  # R5
    13.448,  # R6
    13.492,  # R7
    13.595,  # R8
    14.126,  # R9
    14.173,  # R10
    14.197,  # R11
    14.180,  # R12
    14.017,  # R13
    14.066,  # R14
    14.025,  # R15
    13.991,  # R16
]

FD_I_OLD_SCALE_MASS_MEASUREMENTS_G = {
    "mass_1_g": FD_I_OLD_SCALE_MASS_1_G,
    "mass_2_g": FD_I_OLD_SCALE_MASS_2_G,
}

FD_I_OLD_SCALE_REPEAT_MASS_1_G = {
    "R2.2": 13.475,
    "R9.2": 14.214,
    "R12.2": 14.285,
}

FD_I_OLD_SCALE_REPEAT_MASS_2_G = {
    "R2.2": 13.470,
    "R9.2": 14.217,
    "R12.2": 14.286,
}

FD_I_OLD_SCALE_REPEAT_MASS_MEASUREMENTS_G = {
    "mass_1_g": FD_I_OLD_SCALE_REPEAT_MASS_1_G,
    "mass_2_g": FD_I_OLD_SCALE_REPEAT_MASS_2_G,
}


FD_II_OLD_SCALE_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        82.556,  # RR1
        83.023,  # RR2
        82.284,  # RR3
        82.416,  # RR4
        83.200,  # RR5
        83.047,  # RR6
        83.093,  # RR7
        82.062,  # RR8
        86.030,  # RR9
        85.403,  # RR10
        86.736,  # RR11
        85.744,  # RR12
        85.872,  # RR13
        86.470,  # RR14
        85.715,  # RR15
        86.624,  # RR16
    ],
    "mass_2_g": [
        82.670,  # RR1
        82.844,  # RR2
        81.861,  # RR3
        82.028,  # RR4
        82.788,  # RR5
        82.046,  # RR6
        83.138,  # RR7
        81.796,  # RR8
        85.997,  # RR9
        85.419,  # RR10
        86.623,  # RR11
        85.788,  # RR12
        84.972,  # RR13
        85.857,  # RR14
        85.693,  # RR15
        86.490,  # RR16
    ],
    "mass_3_g": [
        82.777,  # RR1
        82.683,  # RR2
        81.894,  # RR3
        82.255,  # RR4
        82.790,  # RR5
        82.062,  # RR6
        82.875,  # RR7
        81.847,  # RR8
        86.050,  # RR9
        85.470,  # RR10
        86.662,  # RR11
        85.717,  # RR12
        85.321,  # RR13
        86.421,  # RR14
        85.685,  # RR15
        86.169,  # RR16
    ],
    "mass_4_g": [
        82.398,  # RR1
        82.818,  # RR2
        81.901,  # RR3
        81.960,  # RR4
        82.762,  # RR5
        82.078,  # RR6
        82.846,  # RR7
        82.304,  # RR8
        85.682,  # RR9
        85.200,  # RR10
        86.255,  # RR11
        85.462,  # RR12
        85.284,  # RR13
        86.492,  # RR14
        85.681,  # RR15
        86.314,  # RR16
    ],
}


EXTRA_RUN_OLD_SCALE_RUN_LABELS = [
    "LD",
    "HD",
    "E1",
    "E2",
    "E3",
    "E4",
    "E5",
    "E6",
]

EXTRA_RUN_OLD_SCALE_MASS_MEASUREMENTS_G = {
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


BOX_BEHNKEN_OLD_SCALE_RUN_LABELS = [
    "O1",
    "O2",
    "O3",
    "O4",
    "O5",
    "O6",
    "O7",
    "O8",
    "O9",
    "O10",
    "O11",
    "O12",
    "O13",
    "O14",
    "O15",
]

BOX_BEHNKEN_OLD_SCALE_MASS_MEASUREMENTS_G = {
    "mass_1_g": [
        87.239,  # O1
        86.898,  # O2
        86.471,  # O3
        87.345,  # O4
        87.315,  # O5
        87.334,  # O6
        87.519,  # O7
        87.117,  # O8
        86.742,  # O9
        86.994,  # O10
        86.435,  # O11
        87.038,  # O12
        86.702,  # O13
        86.165,  # O14
        86.786,  # O15
    ],
    "mass_2_g": [
        86.848,  # O1
        86.790,  # O2
        86.840,  # O3
        87.186,  # O4
        87.102,  # O5
        87.592,  # O6
        87.402,  # O7
        87.405,  # O8
        86.683,  # O9
        86.862,  # O10
        86.634,  # O11
        86.936,  # O12
        86.705,  # O13
        86.237,  # O14
        86.784,  # O15
    ],
    "mass_3_g": [
        87.020,  # O1
        86.914,  # O2
        86.523,  # O3
        87.099,  # O4
        87.020,  # O5
        87.806,  # O6
        87.767,  # O7
        87.310,  # O8
        86.851,  # O9
        87.150,  # O10
        86.399,  # O11
        86.493,  # O12
        86.175,  # O13
        86.475,  # O14
        None,  # O15
    ],
    "mass_4_g": [
        86.525,  # O1
        86.575,  # O2
        87.019,  # O3
        86.701,  # O4
        86.949,  # O5
        87.819,  # O6
        87.366,  # O7
        86.846,  # O8
        86.394,  # O9
        87.071,  # O10
        85.785,  # O11
        86.603,  # O12
        86.602,  # O13
        86.435,  # O14
        None,  # O15
    ],
}
