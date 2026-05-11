"""Reusable FFD calculations and plots."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FULL_FACTOR_COLUMNS = [
    "A_flow_ratio",
    "B_nozzle_temp",
    "C_speed",
    "D_line_width",
    "E_layer_height",
    "F_walls",
]


def create_fractional_factorial_design(include_walls=True):
    """Create the 16-run coded fractional factorial design."""
    base_design = np.array(
        [
            [a, b, c, d]
            for a in [-1, 1]
            for b in [-1, 1]
            for c in [-1, 1]
            for d in [-1, 1]
        ],
        dtype=int,
    )
    a = base_design[:, 0]
    b = base_design[:, 1]
    c = base_design[:, 2]
    e = a * b
    f = a * c
    design = np.column_stack([base_design, e, f])
    coded_df = pd.DataFrame(design, columns=FULL_FACTOR_COLUMNS)
    coded_df.index += 1
    coded_df.index.name = "Run"

    if not include_walls:
        coded_df = coded_df.drop(columns=["F_walls"])

    return coded_df.astype(int)


def map_coded_to_real(coded_df, levels):
    """Map -1/+1 coded factor values to real parameter values."""
    real_df = coded_df.copy()
    for col, mapping in levels.items():
        real_df[col] = real_df[col].map(mapping)
    return real_df


def style_low_high(real_df, levels):
    """Style low values blue and high values red."""
    def color_low_high(val, col):
        low = levels[col][-1]
        high = levels[col][1]

        if val == high:
            return "background-color: #f8d7da; color: #8b0000"
        if val == low:
            return "background-color: #dbeafe; color: #0b3d91"
        return ""

    styled = real_df.style
    for col in real_df.columns:
        styled = styled.map(lambda v, c=col: color_low_high(v, c), subset=[col])
    return styled


def randomized_second_run_order(random_state=123):
    """Return randomized R*.2 run order while keeping R2.2 first."""
    runs = pd.DataFrame({"Original_Run": [f"R{i}.2" for i in range(1, 17)]})
    first = runs[runs["Original_Run"] == "R2.2"]
    rest = runs[runs["Original_Run"] != "R2.2"].sample(
        frac=1,
        random_state=random_state,
    )
    rep2 = pd.concat([first, rest]).reset_index(drop=True)
    rep2.index += 1
    return rep2


def volume_cm3(volume_mm3):
    return volume_mm3 / 1000


def true_mass_g(volume_mm3, density_g_cm3):
    return volume_cm3(volume_mm3) * density_g_cm3


def calculated_masses_table(density_g_cm3, v_ri_mm3, v_rii_mm3):
    v_ri_cm3 = volume_cm3(v_ri_mm3)
    v_rii_cm3 = volume_cm3(v_rii_mm3)

    return pd.DataFrame({
        "plate": ["RI", "RII"],
        "volume_cm3": [v_ri_cm3, v_rii_cm3],
        "mass_g": [
            v_ri_cm3 * density_g_cm3,
            v_rii_cm3 * density_g_cm3,
        ],
    })


def build_single_mass_table(masses_g, run_prefix="R"):
    mass_measurements_g = pd.DataFrame({
        "run": [f"{run_prefix}{i}" for i in range(1, len(masses_g) + 1)],
        "mass_1_g": masses_g,
    })
    mass_measurements_g["mass_mean_g"] = mass_measurements_g["mass_1_g"]
    mass_measurements_g["mass_std_g"] = np.nan
    return mass_measurements_g


def build_repeat_mass_table(repeat_masses_g):
    return pd.DataFrame({
        "run": list(repeat_masses_g.keys()),
        "mass_1_g": list(repeat_masses_g.values()),
    })


def build_mass_measurements(measurements_g, run_prefix="RR"):
    """Build mass table and calculate trimmed mean from the middle two values."""
    mass_columns = list(measurements_g.keys())
    row_count = max(len(values) for values in measurements_g.values())

    mass_measurements_g = pd.DataFrame({
        "run": [f"{run_prefix}{i}" for i in range(1, row_count + 1)],
        **measurements_g,
    })

    sorted_mass = np.sort(mass_measurements_g[mass_columns].values, axis=1)
    middle_mass = sorted_mass[:, 1:3]

    mass_measurements_g["mass_middle_1_g"] = middle_mass[:, 0]
    mass_measurements_g["mass_middle_2_g"] = middle_mass[:, 1]
    mass_measurements_g["mass_mean_g"] = middle_mass.mean(axis=1).round(1)
    mass_measurements_g["mass_std_g"] = middle_mass.std(axis=1, ddof=1)

    return mass_measurements_g, mass_columns


def calculate_bulk_density(mass_g, volume_cm3_value):
    return mass_g / volume_cm3_value


def calculate_relative_density_percent(mass_g, true_mass_g_value):
    return (mass_g / true_mass_g_value) * 100


def calculate_signed_porosity_percent(mass_g, true_mass_g_value):
    return (1 - (mass_g / true_mass_g_value)) * 100


def calculate_mass_porosity_results(mass_measurements_g, true_mass_g_value, volume_cm3_value):
    results = mass_measurements_g[["run", "mass_mean_g", "mass_std_g"]].copy()
    results["true_mass_g"] = true_mass_g_value
    results["bulk_density_g_cm3"] = results["mass_mean_g"].apply(
        lambda mass: calculate_bulk_density(mass, volume_cm3_value)
    )
    results["relative_density_percent"] = results["mass_mean_g"].apply(
        lambda mass: calculate_relative_density_percent(mass, true_mass_g_value)
    )
    results["signed_porosity_percent"] = results["mass_mean_g"].apply(
        lambda mass: calculate_signed_porosity_percent(mass, true_mass_g_value)
    )

    return results.round({
        "mass_mean_g": 1,
        "mass_std_g": 3,
        "true_mass_g": 3,
        "bulk_density_g_cm3": 3,
        "relative_density_percent": 1,
        "signed_porosity_percent": 1,
    })


def build_result_matrix(coded_df, results_df, response="mass_mean_g"):
    response_df = results_df[[response]].reset_index(drop=True)
    response_df.index = coded_df.index
    return pd.concat([coded_df.copy(), response_df], axis=1)


def compute_main_effect(matrix, effect, response="mass_mean_g"):
    mean_high = matrix.loc[matrix[effect] == 1, response].mean()
    mean_low = matrix.loc[matrix[effect] == -1, response].mean()
    return mean_high - mean_low


def compute_main_effects(matrix, response="mass_mean_g"):
    factor_cols = [col for col in matrix.columns if col not in ["Run", response]]
    return pd.DataFrame({
        "Effect": factor_cols,
        "Main effect": [
            compute_main_effect(matrix, factor, response=response)
            for factor in factor_cols
        ],
    })


def plot_mass_measurements(
    mass_measurements_g,
    mass_columns,
    y_min=None,
    mean_label="Mean",
):
    plot_df = mass_measurements_g.set_index("run")[mass_columns + ["mass_mean_g"]]
    run_labels = plot_df.index.tolist()
    x = np.arange(len(plot_df))
    bar_width = 0.22

    fig, ax = plt.subplots(figsize=(12, 4))
    offsets = np.linspace(
        -bar_width * (len(mass_columns) - 1) / 2,
        bar_width * (len(mass_columns) - 1) / 2,
        len(mass_columns),
    )
    colors = ["#9ecae1", "#4292c6", "#08519c", "#08306b"]

    for offset, column, color in zip(offsets, mass_columns, colors):
        ax.bar(x + offset, plot_df[column], width=bar_width, label=column, color=color)

    ax.plot(
        x,
        plot_df["mass_mean_g"],
        color="#000000",
        marker="o",
        linewidth=2,
        label=mean_label,
    )

    for xi, mass in zip(x, plot_df["mass_mean_g"]):
        if pd.notna(mass):
            ax.annotate(
                f"{mass:.1f}",
                xy=(xi, mass),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_xlabel("Test plate")
    ax.set_ylabel("Mass (g)")
    ax.set_title("Mass measurements")
    if y_min is not None:
        ax.set_ylim(bottom=y_min)
    ax.set_xticks(x)
    ax.set_xticklabels(run_labels, rotation=45)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_porosity_and_relative_density(results_df):
    plot_results = results_df.dropna(
        subset=["relative_density_percent", "signed_porosity_percent"]
    )
    x = np.arange(len(plot_results))

    fig, ax = plt.subplots(figsize=(8, 4))
    porosity_colors = [
        "#9ecae1" if value >= 0 else "#de2d26"
        for value in plot_results["signed_porosity_percent"]
    ]
    ax.bar(x, plot_results["signed_porosity_percent"], color=porosity_colors)

    for xi, porosity in zip(x, plot_results["signed_porosity_percent"]):
        ax.annotate(
            f"{porosity:.1f}%",
            xy=(xi, porosity),
            xytext=(0, -6 if porosity >= 0 else 6),
            textcoords="offset points",
            ha="center",
            va="top" if porosity >= 0 else "bottom",
            fontsize=9,
            color="#333333",
        )

    ax.axhline(0, color="#333333", linewidth=1)
    ax.set_xlabel("Test plate")
    ax.set_ylabel("Signed porosity (%)")
    ax.set_title("Signed porosity from average mass")
    ax.set_xticks(x)
    ax.set_xticklabels(plot_results["run"], rotation=45)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(
        x,
        plot_results["relative_density_percent"],
        color="#08519c",
        marker="o",
        linewidth=2,
    )
    ax.axhline(100, color="#333333", linestyle="-", linewidth=1.3, alpha=0.75, label="True mass (100%)")
    ax.axhline(105, color="#de2d26", linestyle="--", linewidth=1.5, alpha=0.8, label="Flow rate 1.05 (105%)")
    ax.set_xlabel("Test plate")
    ax.set_ylabel("Relative density (% of true mass)")
    ax.set_title("Relative density from average mass")
    ax.set_xticks(x)
    ax.set_xticklabels(plot_results["run"], rotation=45)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_main_effects_bar(main_effect_df, response="mass_mean_g"):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.bar(main_effect_df["Effect"], main_effect_df["Main effect"])
    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylabel(f"Main effect on {response}")
    ax.set_title("Main effects")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.show()


def plot_main_effects_line(matrix, response="mass_mean_g"):
    factor_cols = [col for col in matrix.columns if col not in ["Run", response]]
    overall_mean = matrix[response].mean()

    fig, ax = plt.subplots(figsize=(10, 5))
    xticks = []
    xtick_labels = []

    for i, factor in enumerate(factor_cols):
        mean_low = matrix.loc[matrix[factor] == -1, response].mean()
        mean_high = matrix.loc[matrix[factor] == 1, response].mean()

        x1 = i * 3 + 1
        x2 = i * 3 + 2

        ax.plot([x1, x2], [mean_low, mean_high], marker="o")
        xticks.extend([x1, x2])
        xtick_labels.extend(["-", "+"])

        ax.text(
            (x1 + x2) / 2,
            -0.10,
            factor,
            ha="center",
            va="top",
            fontsize=9,
            fontweight="bold",
            transform=ax.get_xaxis_transform(),
        )

    ax.axhline(overall_mean, linestyle="--", linewidth=1)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.tick_params(axis="x", pad=4)
    ax.set_ylabel("Average Response")
    ax.set_xlabel("Factors", labelpad=22)
    ax.set_title("Main Effects Plot")
    fig.subplots_adjust(bottom=0.18)
    plt.show()

import itertools

def add_interactions(matrix):
    matrix_interactions = matrix.copy()
    
    factor_cols = [col for col in matrix.columns if col not in ["Run", "mass_mean_g"]]
    
    for factor1, factor2 in itertools.combinations(factor_cols, 2):
        interaction_name = f"{factor1} × {factor2}"
        matrix_interactions[interaction_name] = matrix[factor1] * matrix[factor2]
    
    return matrix_interactions

def compute_main_effect(matrix, effect):
    
    if effect not in matrix.columns:
        print(f"Effect '{effect}' not found in matrix columns.")
        return   
    
    mean_high = matrix.loc[matrix[effect] == 1, "mass_mean_g"].mean()
    mean_low = matrix.loc[matrix[effect] == -1, "mass_mean_g"].mean()
    
    main_effect = mean_high - mean_low
    
    print(effect, ":", main_effect)
    
    return main_effect

import pandas as pd
import matplotlib.pyplot as plt

def plot_interaction_effects(interaction_effects):
    interaction_effects_df = pd.DataFrame({
        "Interaction": interaction_effects.keys(),
        "Effect": interaction_effects.values()
    })

    interaction_effects_df = interaction_effects_df.sort_values(
        by="Effect",
        key=abs,
        ascending=False
    )

    plt.figure(figsize=(12, 6))
    plt.bar(interaction_effects_df["Interaction"], interaction_effects_df["Effect"])
    plt.axhline(0, color="black", linewidth=1)

    plt.ylabel("Interaction effect on mass_mean_g")
    plt.title("Interaction Effects Plot")
    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    plt.show()


