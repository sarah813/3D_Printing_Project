"""Reusable FFD calculations and plots."""

import itertools

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

FACTOR_COLORS = {
    "A_flow_ratio": "#1f77b4",
    "B_nozzle_temp": "#ff7f0e",
    "B_nozzle_temp_C": "#ff7f0e",
    "C_speed": "#2ca02c",
    "C_ironing_flow_percent": "#2ca02c",
    "D_line_width": "#d62728",
    "D_ironing_speed_mm_s": "#d62728",
    "E_layer_height": "#9467bd",
    "F_walls": "#8c564b",
}


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


def create_box_behnken_design(factor_levels, center_repeats=3):
    """Create a 3-factor Box-Behnken design with optional center repeats.

    factor_levels must map each factor to {-1: low, 0: middle, 1: high}.
    The design has edge-center points plus repeated all-middle center points.
    """
    factor_cols = list(factor_levels)
    if len(factor_cols) != 3:
        raise ValueError("Box-Behnken design here expects exactly 3 factors.")

    rows = []
    for fixed_factor in factor_cols:
        active_factors = [factor for factor in factor_cols if factor != fixed_factor]
        for low_high in itertools.product([-1, 1], repeat=2):
            row = {factor: 0 for factor in factor_cols}
            row[active_factors[0]] = low_high[0]
            row[active_factors[1]] = low_high[1]
            rows.append(row)

    rows.extend({factor: 0 for factor in factor_cols} for _ in range(center_repeats))
    coded_df = pd.DataFrame(rows, columns=factor_cols)
    coded_df.index += 1
    coded_df.index.name = "Run"

    real_df = map_coded_to_real(coded_df, factor_levels)
    return coded_df, real_df


def create_two_level_factorial_design(factor_levels, center_repeats=1):
    """Create a compact 2-level factorial design plus optional center repeats.

    factor_levels must map each factor to {-1: low, 0: middle, 1: high}.
    The factorial part uses only -1/+1 combinations; center repeats add the
    default all-middle setting for a quick curvature/repeatability check.
    """
    factor_cols = list(factor_levels)
    rows = [
        dict(zip(factor_cols, values))
        for values in itertools.product([-1, 1], repeat=len(factor_cols))
    ]
    rows.extend({factor: 0 for factor in factor_cols} for _ in range(center_repeats))

    coded_df = pd.DataFrame(rows, columns=factor_cols)
    coded_df.index += 1
    coded_df.index.name = "Run"

    real_df = map_coded_to_real(coded_df, factor_levels)
    return coded_df, real_df


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


def ideal_mass_g(volume_mm3, density_g_cm3):
    """Mass expected from specimen volume and material density."""
    return true_mass_g(volume_mm3, density_g_cm3)


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


def calculate_mean_mass_from_measurements(mass_values_g):
    """Calculate the representative mean mass for one run.

    Current rule for four measurements:
    sort the values, discard the lowest and highest value, and average the
    middle two. This reduces the influence of one unusually low/high reading.
    """
    valid_masses = pd.Series(mass_values_g).dropna().to_numpy(dtype=float)
    measurement_count = len(valid_masses)

    if measurement_count == 0:
        return np.nan

    if measurement_count == 4:
        sorted_masses = np.sort(valid_masses)
        middle_masses = sorted_masses[1:3]
        return middle_masses.mean()

    if measurement_count == 2 or measurement_count == 3:
        return valid_masses.mean()
    
    if measurement_count == 1:
        return valid_masses[0]

    raise ValueError(
        f"Expected 1 to 4 mass measurements, got {measurement_count}."
    )


def calculate_mass_std_from_measurements(mass_values_g):
    """Calculate the spread of the measurements used for the mean mass.

    For four measurements this matches the trimmed-mean rule and uses only the
    middle two values. For two or three measurements it uses all available
    values because there is not enough data to safely trim high/low outliers.
    """
    valid_masses = pd.Series(mass_values_g).dropna().to_numpy(dtype=float)
    measurement_count = len(valid_masses)

    if measurement_count <= 1:
        return np.nan
    if measurement_count == 4:
        values_for_std = np.sort(valid_masses)[1:3]
    else:
        values_for_std = valid_masses

    return values_for_std.std(ddof=1)


def middle_two_masses_or_nan(mass_values_g):
    """Return the middle two masses only when four measurements are available."""
    valid_masses = pd.Series(mass_values_g).dropna().to_numpy(dtype=float)
    if len(valid_masses) != 4:
        return pd.Series([np.nan, np.nan])

    middle_masses = np.sort(valid_masses)[1:3]
    return pd.Series(middle_masses)


def build_mass_measurements(
    measurements_g,
    run_prefix="RR",
    run_labels=None,
    mean_round_digits=1,
):
    """Build mass table and calculate the representative mean mass."""
    mass_columns = list(measurements_g.keys())
    row_count = max(len(values) for values in measurements_g.values())
    if run_labels is None:
        run_labels = [f"{run_prefix}{i}" for i in range(1, row_count + 1)]
    elif len(run_labels) != row_count:
        raise ValueError("run_labels must have the same length as measurements.")

    mass_measurements_g = pd.DataFrame({
        "run": run_labels,
        **measurements_g,
    })

    # These columns are filled only for four-measurement rows, where the mean
    # uses the middle two values. They stay NaN for one/two/three measurements.
    mass_measurements_g[["mass_middle_1_g", "mass_middle_2_g"]] = (
        mass_measurements_g[mass_columns].apply(middle_two_masses_or_nan, axis=1)
    )
    mass_measurements_g["mass_mean_g"] = mass_measurements_g[mass_columns].apply(
        calculate_mean_mass_from_measurements,
        axis=1,
    ).round(mean_round_digits)
    mass_measurements_g["mass_std_g"] = mass_measurements_g[mass_columns].apply(
        calculate_mass_std_from_measurements,
        axis=1,
    )

    return mass_measurements_g, mass_columns


def build_repeat_mass_measurements(repeat_measurements_g, mean_round_digits=1):
    """Build repeat-run mass table and calculate the repeat mean mass."""
    mass_columns = list(repeat_measurements_g.keys())
    run_names = list(next(iter(repeat_measurements_g.values())).keys())

    repeat_measurements_g = pd.DataFrame({
        "run": run_names,
        **{
            column: [measurements.get(run) for run in run_names]
            for column, measurements in repeat_measurements_g.items()
        },
    })
    repeat_measurements_g["mass_mean_g"] = repeat_measurements_g[mass_columns].apply(
        calculate_mean_mass_from_measurements,
        axis=1,
    ).round(mean_round_digits)
    repeat_measurements_g["mass_std_g"] = repeat_measurements_g[mass_columns].apply(
        calculate_mass_std_from_measurements,
        axis=1,
    )

    return repeat_measurements_g, mass_columns


def calculate_repeatability(original_measurements_g, repeat_measurements_g):
    """Compare repeated runs against their original run mean masses.

    Repeat run names are expected to include the original run before the dot,
    for example R2.2 is compared with original run R2.
    """
    repeat_df = repeat_measurements_g.copy()
    if "mass_mean_g" not in repeat_df.columns:
        mass_columns = [
            col for col in repeat_df.columns
            if col.startswith("mass_") and col.endswith("_g")
        ]
        repeat_df["mass_mean_g"] = repeat_df[mass_columns].apply(
            calculate_mean_mass_from_measurements,
            axis=1,
        )
    repeat_df = repeat_df.dropna(subset=["mass_mean_g"])

    original_df = original_measurements_g[["run", "mass_mean_g"]].rename(
        columns={"mass_mean_g": "original_mass_mean_g"}
    )
    repeat_df = repeat_df[["run", "mass_mean_g"]].rename(
        columns={"mass_mean_g": "repeat_mass_mean_g"}
    )
    repeat_df["original_run"] = repeat_df["run"].str.split(".").str[0]
    repeat_df = (
        repeat_df.groupby("original_run", as_index=False)["repeat_mass_mean_g"]
        .mean()
    )
    repeat_df["run_repeat"] = repeat_df["original_run"] + " repeat mean"

    repeatability_df = repeat_df.merge(
        original_df,
        left_on="original_run",
        right_on="run",
        suffixes=("_repeat", "_original"),
    )
    repeatability_df["difference_g"] = (
        repeatability_df["repeat_mass_mean_g"]
        - repeatability_df["original_mass_mean_g"]
    )
    repeatability_df["abs_difference_g"] = repeatability_df["difference_g"].abs()
    repeatability_df["difference_percent"] = (
        repeatability_df["difference_g"]
        / repeatability_df["original_mass_mean_g"]
        * 100
    )

    return repeatability_df[[
        "run_repeat",
        "original_run",
        "original_mass_mean_g",
        "repeat_mass_mean_g",
        "difference_g",
        "abs_difference_g",
        "difference_percent",
    ]].round({
        "original_mass_mean_g": 3,
        "repeat_mass_mean_g": 3,
        "difference_g": 3,
        "abs_difference_g": 3,
        "difference_percent": 2,
    })


def calculate_center_point_repeatability(
    coded_df,
    mass_measurements_g,
    response="mass_mean_g",
    center_level=0,
):
    """Estimate process repeatability from repeated center-point runs.

    Box-Behnken designs usually include repeated center points. These have the
    same coded factor settings, so variation between their mean masses is an
    estimate of print-to-print repeatability at the center condition.
    """
    factor_cols = list(coded_df.columns)
    repeat_df = coded_df.copy()
    repeat_df["run"] = mass_measurements_g["run"].to_numpy()
    repeat_df[response] = mass_measurements_g[response].to_numpy()
    if "mass_std_g" in mass_measurements_g.columns:
        repeat_df["mass_std_g"] = mass_measurements_g["mass_std_g"].to_numpy()

    center_mask = repeat_df[factor_cols].eq(center_level).all(axis=1)
    center_df = repeat_df.loc[center_mask].dropna(subset=[response]).copy()

    if center_df.empty:
        summary_df = pd.DataFrame([{
            "center_repeats_with_data": 0,
            "center_mean_g": np.nan,
            "between_print_std_g": np.nan,
            "center_range_g": np.nan,
            "max_pair_difference_g": np.nan,
            "max_pair_difference_percent": np.nan,
            "mean_within_print_std_g": np.nan,
        }])
        return center_df, summary_df

    center_values = center_df[response].to_numpy(dtype=float)
    center_mean = center_values.mean()
    center_range = center_values.max() - center_values.min()
    between_std = np.nan if len(center_values) < 2 else center_values.std(ddof=1)
    mean_within_std = (
        center_df["mass_std_g"].mean()
        if "mass_std_g" in center_df.columns
        else np.nan
    )

    center_df["difference_from_center_mean_g"] = center_df[response] - center_mean
    center_df["abs_difference_from_center_mean_g"] = (
        center_df["difference_from_center_mean_g"].abs()
    )

    summary_df = pd.DataFrame([{
        "center_repeats_with_data": len(center_df),
        "center_mean_g": center_mean,
        "between_print_std_g": between_std,
        "center_range_g": center_range,
        "max_pair_difference_g": center_range,
        "max_pair_difference_percent": center_range / center_mean * 100,
        "mean_within_print_std_g": mean_within_std,
    }])

    return (
        center_df.round({
            response: 3,
            "mass_std_g": 3,
            "difference_from_center_mean_g": 3,
            "abs_difference_from_center_mean_g": 3,
        }),
        summary_df.round({
            "center_mean_g": 3,
            "between_print_std_g": 3,
            "center_range_g": 3,
            "max_pair_difference_g": 3,
            "max_pair_difference_percent": 2,
            "mean_within_print_std_g": 3,
        }),
    )


def calculate_center_curvature(
    coded_df,
    mass_measurements_g,
    response="mass_mean_g",
    center_level=0,
):
    """Estimate overall curvature from factorial and center-point means.

    A two-level design cannot separate the quadratic effect of each factor.
    This calculation therefore reports only one combined curvature estimate.
    """
    factor_cols = list(coded_df.columns)
    analysis_df = coded_df.copy()
    analysis_df[response] = mass_measurements_g[response].to_numpy()

    center_mask = analysis_df[factor_cols].eq(center_level).all(axis=1)
    factorial_mask = analysis_df[factor_cols].abs().eq(1).all(axis=1)
    center_values = analysis_df.loc[center_mask, response].dropna()
    factorial_values = analysis_df.loc[factorial_mask, response].dropna()

    center_mean = center_values.mean() if len(center_values) else np.nan
    factorial_mean = factorial_values.mean() if len(factorial_values) else np.nan
    curvature = (
        center_mean - factorial_mean
        if pd.notna(center_mean) and pd.notna(factorial_mean)
        else np.nan
    )

    return pd.DataFrame([{
        "factorial_runs_with_data": len(factorial_values),
        "center_runs_with_data": len(center_values),
        "factorial_mean_g": factorial_mean,
        "center_mean_g": center_mean,
        "overall_curvature_g": curvature,
        "center_std_g": (
            center_values.std(ddof=1) if len(center_values) > 1 else np.nan
        ),
    }])


def plot_center_curvature(curvature_summary, experiment=None):
    """Plot factorial and center means for an overall curvature check."""
    row = curvature_summary.iloc[0]
    if pd.isna(row["factorial_mean_g"]) or pd.isna(row["center_mean_g"]):
        return

    values = [row["factorial_mean_g"], row["center_mean_g"]]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(
        ["Factorial mean", "Center mean"],
        values,
        color=["#4c78a8", "#f28e2b"],
        width=0.55,
    )
    for index, value in enumerate(values):
        ax.annotate(
            f"{value:.3f} g",
            (index, value),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
        )
    ax.set_ylabel("Mass (g)")
    ax.set_title(experiment_title("Overall center-point curvature", experiment))
    ax.set_ylim(0, max(values) * 1.08)
    ax.text(
        0.98,
        0.05,
        f"Center - factorial mean = {row['overall_curvature_g']:+.3f} g",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "#bbbbbb",
            "alpha": 0.9,
        },
    )
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_center_point_repeatability(
    center_repeatability_df,
    summary_df,
    response="mass_mean_g",
    experiment=None,
):
    """Plot repeated center-point masses with mean and repeatability band."""
    if center_repeatability_df.empty:
        return

    fig, ax = plt.subplots(figsize=(6, 3.8))
    x = np.arange(len(center_repeatability_df))
    center_mean = summary_df.loc[0, "center_mean_g"]
    between_std = summary_df.loc[0, "between_print_std_g"]

    ax.plot(
        x,
        center_repeatability_df[response],
        marker="o",
        color="#08519c",
        linewidth=1.8,
    )
    ax.axhline(
        center_mean,
        color="#333333",
        linestyle="--",
        linewidth=1.2,
        label=f"center mean = {center_mean:.3f} g",
    )

    if pd.notna(between_std):
        ax.fill_between(
            [-0.4, len(center_repeatability_df) - 0.6],
            center_mean - between_std,
            center_mean + between_std,
            color="#9ecae1",
            alpha=0.25,
            label=f"+/-1 repeat std = {between_std:.3f} g",
        )

    for xi, (_, row) in zip(x, center_repeatability_df.iterrows()):
        ax.annotate(
            f"{row[response]:.3f}",
            xy=(xi, row[response]),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(center_repeatability_df["run"])
    ax.set_ylabel("Mass (g)")
    ax.set_title(experiment_title("Center-point repeatability", experiment))
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    plt.tight_layout()
    plt.show(block=False)


def plot_repeatability_summary_comparison(summary_df, experiment=None):
    """Compare repeatability summary metrics across measurement sets."""
    metric_cols = [
        "between_print_std_g",
        "center_range_g",
        "mean_within_print_std_g",
    ]
    available_metrics = [
        metric for metric in metric_cols
        if metric in summary_df.columns and summary_df[metric].notna().any()
    ]
    if not available_metrics:
        return

    plot_df = summary_df.set_index("measurement_set")[available_metrics]
    ax = plot_df.plot(kind="bar", figsize=(7, 4), color=["#08519c", "#f28e2b", "#6baed6"])
    ax.set_ylabel("Mass variation (g)")
    ax.set_title(experiment_title("Center-point repeatability comparison", experiment))
    ax.tick_params(axis="x", rotation=0)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def calculate_bulk_density(mass_g, volume_cm3_value):
    return mass_g / volume_cm3_value


def calculate_relative_density_percent(mass_g, true_mass_g_value):
    return (mass_g / true_mass_g_value) * 100


def calculate_signed_porosity_percent(mass_g, true_mass_g_value):
    return (1 - (mass_g / true_mass_g_value)) * 100


def calculate_mass_deviation_from_ideal_percent(mass_g, ideal_mass_g_value):
    """Positive means measured mass is below ideal; negative means above ideal."""
    return (1 - (mass_g / ideal_mass_g_value)) * 100


def calculate_mass_deficit_from_ideal_percent(mass_g, ideal_mass_g_value):
    """Deprecated name; use calculate_mass_deviation_from_ideal_percent."""
    return calculate_mass_deviation_from_ideal_percent(mass_g, ideal_mass_g_value)


def calculate_measured_to_ideal_mass_percent(mass_g, ideal_mass_g_value):
    return (mass_g / ideal_mass_g_value) * 100


def calculate_ideal_mass_deviation_results(
    mass_measurements_g,
    ideal_mass_g_value,
    volume_cm3_value,
):
    results = mass_measurements_g[["run", "mass_mean_g", "mass_std_g"]].copy()
    results["ideal_mass_g"] = ideal_mass_g_value
    results["bulk_density_g_cm3"] = results["mass_mean_g"].apply(
        lambda mass: calculate_bulk_density(mass, volume_cm3_value)
    )
    results["measured_to_ideal_mass_percent"] = results["mass_mean_g"].apply(
        lambda mass: calculate_measured_to_ideal_mass_percent(mass, ideal_mass_g_value)
    )
    results["mass_deviation_from_ideal_percent"] = results["mass_mean_g"].apply(
        lambda mass: calculate_mass_deviation_from_ideal_percent(mass, ideal_mass_g_value)
    )
    # Backward-compatible column names used by earlier notebook cells.
    results["true_mass_g"] = ideal_mass_g_value
    results["relative_density_percent"] = results["measured_to_ideal_mass_percent"]
    results["mass_deficit_from_ideal_percent"] = results["mass_deviation_from_ideal_percent"]
    results["signed_porosity_percent"] = results["mass_deviation_from_ideal_percent"]

    return results.round({
        "mass_mean_g": 1,
        "mass_std_g": 3,
        "ideal_mass_g": 3,
        "true_mass_g": 3,
        "bulk_density_g_cm3": 3,
        "measured_to_ideal_mass_percent": 1,
        "mass_deviation_from_ideal_percent": 1,
        "mass_deficit_from_ideal_percent": 1,
        "relative_density_percent": 1,
        "signed_porosity_percent": 1,
    })


def calculate_mass_porosity_results(mass_measurements_g, true_mass_g_value, volume_cm3_value):
    """Deprecated name; use calculate_ideal_mass_deviation_results."""
    return calculate_ideal_mass_deviation_results(
        mass_measurements_g,
        true_mass_g_value,
        volume_cm3_value,
    )


def build_result_matrix(coded_df, results_df, response="mass_mean_g"):
    response_df = results_df[[response]].reset_index(drop=True)
    response_df.index = coded_df.index
    return pd.concat([coded_df.copy(), response_df], axis=1)


def build_ironing_quality_table(
    blob_counts,
    run_labels=None,
    stringing_runs=None,
    stringing_penalty=5,
):
    """Build a visual-defect score table for ironing runs.

    The quality score is intentionally simple and transparent:
    blob count plus a configurable penalty for visible stringing.
    Lower scores are better.
    """
    stringing_runs = set(stringing_runs or [])
    if isinstance(blob_counts, dict):
        labels = list(run_labels) if run_labels is not None else list(blob_counts)
        counts = [blob_counts.get(label) for label in labels]
    else:
        labels = list(run_labels) if run_labels is not None else [
            f"I{i}" for i in range(1, len(blob_counts) + 1)
        ]
        counts = list(blob_counts)

    rows = []
    for label, count in zip(labels, counts):
        has_stringing = label in stringing_runs
        rows.append({
            "run": label,
            "blob_count": count,
            "has_stringing": has_stringing,
            "stringing_penalty": stringing_penalty if has_stringing else 0,
            "quality_score": (
                np.nan
                if pd.isna(count)
                else count + (stringing_penalty if has_stringing else 0)
            ),
        })
    return pd.DataFrame(rows)


def build_mass_quality_tradeoff(
    mass_measurements_g,
    quality_table,
    target_mass_g=None,
    mass_weight=0.5,
    quality_weight=0.5,
):
    """Combine ironing mass and visual-quality results into a tradeoff table.

    If target_mass_g is supplied, the mass part rewards closeness to target.
    Otherwise it rewards higher measured mass. The combined score is scaled
    from the measured data only; lower combined scores are better.
    """
    tradeoff_df = mass_measurements_g[["run", "mass_mean_g"]].merge(
        quality_table[["run", "blob_count", "has_stringing", "quality_score"]],
        on="run",
        how="outer",
    )

    if target_mass_g is None:
        tradeoff_df["mass_metric"] = -tradeoff_df["mass_mean_g"]
        tradeoff_df["mass_metric_label"] = "higher mass is better"
    else:
        tradeoff_df["mass_error_abs_g"] = (
            tradeoff_df["mass_mean_g"] - target_mass_g
        ).abs()
        tradeoff_df["mass_metric"] = tradeoff_df["mass_error_abs_g"]
        tradeoff_df["mass_metric_label"] = "closer to target mass is better"

    def normalize_lower_is_better(series):
        valid = series.dropna()
        if valid.empty:
            return pd.Series(np.nan, index=series.index)
        spread = valid.max() - valid.min()
        if spread == 0:
            return pd.Series(0.0, index=series.index)
        return (series - valid.min()) / spread

    tradeoff_df["mass_score_0_best"] = normalize_lower_is_better(
        tradeoff_df["mass_metric"]
    )
    tradeoff_df["quality_score_0_best"] = normalize_lower_is_better(
        tradeoff_df["quality_score"]
    )
    tradeoff_df["combined_score_0_best"] = (
        mass_weight * tradeoff_df["mass_score_0_best"]
        + quality_weight * tradeoff_df["quality_score_0_best"]
    )
    tradeoff_df["combined_rank"] = tradeoff_df["combined_score_0_best"].rank(
        method="min",
    )
    return tradeoff_df.sort_values("combined_score_0_best")


def build_desirability_table(
    response_df,
    mass_col="mass_mean_g",
    quality_col="quality_score",
    mass_low=None,
    mass_high=None,
    quality_best=None,
    quality_worst=None,
    mass_weight=1.0,
    quality_weight=1.0,
):
    """Calculate mass, quality, and combined response desirabilities.

    Mass uses a higher-is-better desirability and visual quality uses a
    lower-is-better desirability. The overall desirability is their weighted
    geometric mean. All desirabilities range from 0 to 1.
    """
    result_df = response_df.copy()
    valid_mass = result_df[mass_col].dropna()
    valid_quality = result_df[quality_col].dropna()
    if valid_mass.empty or valid_quality.empty:
        raise ValueError("Mass and quality responses are required.")

    mass_low = valid_mass.min() if mass_low is None else mass_low
    mass_high = valid_mass.max() if mass_high is None else mass_high
    quality_best = valid_quality.min() if quality_best is None else quality_best
    quality_worst = valid_quality.max() if quality_worst is None else quality_worst

    if mass_high <= mass_low:
        raise ValueError("mass_high must be greater than mass_low.")
    if quality_worst <= quality_best:
        raise ValueError("quality_worst must be greater than quality_best.")
    if mass_weight <= 0 or quality_weight <= 0:
        raise ValueError("Desirability weights must be positive.")

    result_df["mass_desirability"] = (
        (result_df[mass_col] - mass_low) / (mass_high - mass_low)
    ).clip(0, 1)
    result_df["quality_desirability"] = (
        (quality_worst - result_df[quality_col])
        / (quality_worst - quality_best)
    ).clip(0, 1)

    total_weight = mass_weight + quality_weight
    result_df["overall_desirability"] = (
        result_df["mass_desirability"].pow(mass_weight)
        * result_df["quality_desirability"].pow(quality_weight)
    ).pow(1 / total_weight)
    result_df["desirability_rank"] = result_df["overall_desirability"].rank(
        method="min",
        ascending=False,
    )
    result_df.attrs["desirability_limits"] = {
        "mass_low": mass_low,
        "mass_high": mass_high,
        "quality_best": quality_best,
        "quality_worst": quality_worst,
        "mass_weight": mass_weight,
        "quality_weight": quality_weight,
    }
    return result_df.sort_values("overall_desirability", ascending=False)


def calculate_score_repeatability(
    original_table,
    repeat_table,
    score_col,
    run_col="run",
):
    """Compare repeated non-mass scores such as visual quality scores."""
    if repeat_table.empty:
        return pd.DataFrame()

    comparison_df = repeat_table[[run_col, score_col]].dropna(
        subset=[score_col]
    ).copy()
    if comparison_df.empty:
        return pd.DataFrame()

    comparison_df["original_run"] = comparison_df[run_col].str.replace(
        r"\.\d+$",
        "",
        regex=True,
    )
    comparison_df = comparison_df.merge(
        original_table[[run_col, score_col]].rename(
            columns={
                run_col: "original_run",
                score_col: f"original_{score_col}",
            }
        ),
        on="original_run",
        how="left",
    )
    comparison_df = comparison_df.rename(columns={score_col: f"repeat_{score_col}"})
    comparison_df["difference"] = (
        comparison_df[f"repeat_{score_col}"]
        - comparison_df[f"original_{score_col}"]
    )
    comparison_df["abs_difference"] = comparison_df["difference"].abs()
    return comparison_df


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


def build_factor_level_coverage(matrix, response="mass_mean_g", factor_cols=None):
    """Count available responses at each coded level for every factor."""
    if factor_cols is None:
        factor_cols = [col for col in matrix.columns if col not in ["Run", response]]

    rows = []
    for factor in factor_cols:
        for level in [-1, 0, 1]:
            rows.append({
                "factor": factor,
                "coded_level": level,
                "available_runs": int(
                    matrix.loc[matrix[factor] == level, response].notna().sum()
                ),
            })
    return pd.DataFrame(rows)


def fit_main_effect_screening_model(matrix, response="mass_mean_g", factor_cols=None):
    """Fit an intercept and coded main effects using all available responses.

    This supports incomplete screening designs. The model is estimable only
    when its design matrix has full column rank. Interactions and curvature are
    intentionally excluded.
    """
    if factor_cols is None:
        factor_cols = [col for col in matrix.columns if col not in ["Run", response]]

    analysis_df = matrix.dropna(subset=[response]).copy()
    design_matrix = pd.DataFrame(
        {"Intercept": np.ones(len(analysis_df))},
        index=analysis_df.index,
    )
    for factor in factor_cols:
        design_matrix[factor] = analysis_df[factor].astype(float)

    y = analysis_df[response].to_numpy(dtype=float)
    x = design_matrix.to_numpy(dtype=float)
    coefficients, _, rank, singular_values = np.linalg.lstsq(x, y, rcond=None)
    predicted = x @ coefficients
    residual = y - predicted
    residual_df = len(analysis_df) - rank
    residual_ss = float(np.sum(residual**2))
    residual_ms = residual_ss / residual_df if residual_df > 0 else np.nan
    total_ss = float(np.sum((y - y.mean()) ** 2))
    r_squared = np.nan if total_ss == 0 else 1 - residual_ss / total_ss

    coefficient_series = pd.Series(coefficients, index=design_matrix.columns)
    statistics_df = coefficient_series.rename("coefficient").reset_index()
    statistics_df = statistics_df.rename(columns={"index": "term"})
    statistics_df["effect"] = statistics_df["coefficient"]
    statistics_df.loc[
        statistics_df["term"] != "Intercept",
        "effect",
    ] *= 2

    if residual_df > 0 and rank == design_matrix.shape[1]:
        xtx_inv = np.linalg.pinv(x.T @ x)
        standard_errors = np.sqrt(np.diag(xtx_inv) * residual_ms)
        statistics_df["standard_error"] = standard_errors
        statistics_df["t_value"] = (
            statistics_df["coefficient"] / statistics_df["standard_error"]
        )
        try:
            from scipy import stats

            statistics_df["p_value"] = 2 * stats.t.sf(
                statistics_df["t_value"].abs(),
                df=residual_df,
            )
        except Exception:
            statistics_df["p_value"] = np.nan
    else:
        statistics_df["standard_error"] = np.nan
        statistics_df["t_value"] = np.nan
        statistics_df["p_value"] = np.nan

    fitted_df = analysis_df.copy()
    fitted_df["predicted_response"] = predicted
    fitted_df["residual"] = residual

    info_df = pd.DataFrame([{
        "available_runs": len(analysis_df),
        "required_model_terms": design_matrix.shape[1],
        "rank": rank,
        "full_rank": rank == design_matrix.shape[1],
        "residual_df": residual_df,
        "r_squared": r_squared,
        "residual_std": np.sqrt(residual_ms) if residual_df > 0 else np.nan,
    }])

    return {
        "factor_cols": factor_cols,
        "response": response,
        "coefficients": coefficient_series,
        "statistics_df": statistics_df,
        "fitted_df": fitted_df,
        "info_df": info_df,
        "singular_values": singular_values,
    }


def plot_screening_model_effects(model, experiment=None):
    """Plot signed coded main effects from a screening model."""
    effect_df = model["statistics_df"].copy()
    effect_df = effect_df[effect_df["term"] != "Intercept"]
    effect_df = effect_df.sort_values("effect", key=abs, ascending=True)
    colors = [FACTOR_COLORS.get(term, "#4c78a8") for term in effect_df["term"]]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.barh(effect_df["term"], effect_df["effect"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel(
        f"Estimated main effect on {model['response']} "
        "(+ means higher response at high factor level)"
    )
    ax.set_title(experiment_title("Main-effect screening model", experiment))
    ax.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_combined_screening_model_effects(
    models,
    response_labels=None,
    normalize=True,
    experiment=None,
):
    """Plot comparable signed main effects for multiple screening responses."""
    if response_labels is None:
        response_labels = {
            label: model["response"]
            for label, model in models.items()
        }

    effect_frames = []
    for label, model in models.items():
        effects = model["statistics_df"].copy()
        effects = effects[effects["term"] != "Intercept"]
        effects = effects.rename(columns={"term": "factor", "effect": "main_effect"})

        scale = 1.0
        axis_label = "Main effect"
        if normalize:
            fitted = model["fitted_df"][model["response"]].dropna()
            response_range = fitted.max() - fitted.min()
            if pd.notna(response_range) and response_range != 0:
                scale = response_range
            axis_label = "Main effect / observed response range"

        effects["response"] = response_labels.get(label, label)
        effects["plot_effect"] = effects["main_effect"] / scale
        effect_frames.append(effects[["factor", "response", "plot_effect"]])

    effect_df = pd.concat(effect_frames, ignore_index=True)
    factors = list(models[next(iter(models))]["factor_cols"])
    responses = [
        response_labels.get(label, label)
        for label in models
    ]

    x = np.arange(len(factors))
    bar_width = min(0.34, 0.8 / max(len(responses), 1))
    offsets = np.linspace(
        -bar_width * (len(responses) - 1) / 2,
        bar_width * (len(responses) - 1) / 2,
        len(responses),
    )
    colors = ["#08519c", "#de2d26", "#31a354", "#756bb1"]

    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    for offset, response, color in zip(offsets, responses, colors):
        response_df = effect_df[effect_df["response"] == response].set_index("factor")
        values = [
            response_df.loc[factor, "plot_effect"]
            if factor in response_df.index
            else np.nan
            for factor in factors
        ]
        ax.bar(
            x + offset,
            values,
            width=bar_width,
            label=response,
            color=color,
            alpha=0.88,
        )

    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(factors, rotation=35, ha="right")
    ax.set_ylabel(axis_label)
    ax.set_title(experiment_title("Combined main-effect screening model", experiment))
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_screening_observed_vs_predicted(model, experiment=None):
    """Plot observed versus predicted responses for a screening model."""
    fitted_df = model["fitted_df"]
    response = model["response"]

    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(
        fitted_df[response],
        fitted_df["predicted_response"],
        color="#08519c",
    )
    minimum = min(fitted_df[response].min(), fitted_df["predicted_response"].min())
    maximum = max(fitted_df[response].max(), fitted_df["predicted_response"].max())
    ax.plot([minimum, maximum], [minimum, maximum], color="#333333", linestyle="--")

    for run, observed, predicted in zip(
        fitted_df.index,
        fitted_df[response],
        fitted_df["predicted_response"],
    ):
        ax.annotate(
            str(run),
            (observed, predicted),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(f"Observed {response}")
    ax.set_ylabel(f"Predicted {response}")
    ax.set_title(experiment_title("Observed vs predicted", experiment))
    ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_quality_counts(quality_table, experiment=None):
    """Plot visual-defect counts and stringing penalty for ironing runs."""
    plot_df = quality_table.dropna(subset=["quality_score"]).copy()
    x = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(
        x,
        plot_df["blob_count"],
        color="#74a9cf",
        label="Blob count",
    )
    ax.bar(
        x,
        plot_df["stringing_penalty"].fillna(0),
        bottom=plot_df["blob_count"],
        color="#ef8a62",
        label="Stringing penalty",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["run"], rotation=45)
    ax.set_ylabel("Visual defect score")
    ax.set_title(experiment_title("Blob and stringing score", experiment))
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_mass_quality_tradeoff(tradeoff_df, target_mass_g=None, experiment=None):
    """Plot mass versus visual-quality score for tradeoff decisions."""
    plot_df = tradeoff_df.dropna(subset=["mass_mean_g", "quality_score"]).copy()
    if "overall_desirability" in plot_df.columns:
        color_values = plot_df["overall_desirability"]
        color_label = "Overall desirability (higher is better)"
        cmap = "viridis"
    else:
        color_values = plot_df["combined_score_0_best"]
        color_label = "Combined score (lower is better)"
        cmap = "viridis_r"

    fig, ax = plt.subplots(figsize=(7.2, 5))
    scatter = ax.scatter(
        plot_df["mass_mean_g"],
        plot_df["quality_score"],
        c=color_values,
        cmap=cmap,
        s=85,
        edgecolor="#333333",
    )
    for _, row in plot_df.iterrows():
        ax.annotate(
            row["run"],
            (row["mass_mean_g"], row["quality_score"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9,
        )
    if target_mass_g is not None:
        ax.axvline(
            target_mass_g,
            color="#f28e2b",
            linewidth=1.5,
            alpha=0.75,
            label="Target mass",
        )
        ax.legend()
    ax.set_xlabel("Mass mean (g)")
    ax.set_ylabel("Visual defect score (lower is better)")
    ax.set_title(experiment_title("Mass and visual-quality tradeoff", experiment))
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(color_label)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_desirability_ranking(desirability_df, experiment=None):
    """Plot overall desirability ranked from best to worst."""
    plot_df = desirability_df.dropna(subset=["overall_desirability"]).copy()
    plot_df = plot_df.sort_values("overall_desirability", ascending=True)
    colors = plt.cm.viridis(plot_df["overall_desirability"])

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.barh(
        plot_df["run"],
        plot_df["overall_desirability"],
        color=colors,
    )
    ax.set_xlim(0, 1)
    ax.set_xlabel("Overall desirability (0 = unacceptable, 1 = best)")
    ax.set_ylabel("Ironing run")
    ax.set_title(experiment_title("Mass and quality desirability", experiment))
    ax.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def predict_response_from_main_effects(
    coded_setting,
    grand_mean,
    main_effects,
    factor_cols,
):
    """Predict a response from an additive main-effects model.

    This is the standard screening-model estimate:
    predicted response = grand mean + sum(coded level * effect / 2).
    It does not include interactions or curvature, so use it for planning or
    confirmation comparisons rather than as proof of the true optimum.
    """
    return grand_mean + sum(
        coded_setting[factor] * main_effects[factor] / 2
        for factor in factor_cols
    )


def real_settings_from_coded(coded_setting, levels, factor_cols=None):
    """Convert a coded setting dictionary into real printer settings."""
    if factor_cols is None:
        factor_cols = list(levels)

    return {
        factor: levels[factor][coded_setting[factor]]
        for factor in factor_cols
    }


def build_screening_confirmation_comparison(
    model,
    confirmation_runs,
    levels,
    response="mass_mean_g",
):
    """Compare measured confirmation runs with a fitted screening model."""
    confirmation_df = pd.DataFrame(confirmation_runs).copy()
    factor_cols = model["factor_cols"]

    def real_to_coded(factor, real_value):
        mapping = levels[factor]
        for coded_value, mapped_value in mapping.items():
            if np.isclose(real_value, mapped_value):
                return float(coded_value)
        coded_levels = np.array(sorted(mapping), dtype=float)
        real_levels = np.array(
            [mapping[coded_value] for coded_value in coded_levels],
            dtype=float,
        )
        return float(np.interp(real_value, real_levels, coded_levels))

    coded_df = pd.DataFrame(index=confirmation_df.index)
    for factor in factor_cols:
        coded_df[factor] = confirmation_df[factor].apply(
            lambda value: real_to_coded(factor, value)
        )

    coefficients = model["coefficients"]
    predicted = np.full(len(coded_df), coefficients["Intercept"], dtype=float)
    for factor in factor_cols:
        predicted += coded_df[factor].to_numpy() * coefficients[factor]

    result_df = confirmation_df.copy()
    for factor in factor_cols:
        result_df[f"{factor}_coded"] = coded_df[factor]
    result_df["predicted_response"] = predicted
    result_df["prediction_error"] = result_df[response] - predicted
    result_df["abs_prediction_error"] = result_df["prediction_error"].abs()
    return result_df


def plot_screening_confirmation_comparison(
    comparison_df,
    response="mass_mean_g",
    target_value=None,
    experiment=None,
):
    """Plot observed and screening-model-predicted confirmation responses."""
    plot_df = comparison_df.dropna(
        subset=[response, "predicted_response"]
    ).copy()
    x = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(
        x,
        plot_df[response],
        color="#08519c",
        marker="o",
        linewidth=2,
        label="Observed",
    )
    ax.plot(
        x,
        plot_df["predicted_response"],
        color="#de2d26",
        marker="s",
        linewidth=2,
        label="Main-effect prediction",
    )
    if target_value is not None:
        ax.axhline(
            target_value,
            color="#f28e2b",
            linewidth=1.5,
            alpha=0.65,
            label="Target mass",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["run"])
    ax.set_ylabel(response)
    ax.set_title(experiment_title("Confirmation runs", experiment))
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def build_main_effect_prediction_comparison(
    result_matrix,
    main_effect_df,
    levels,
    response="mass_mean_g",
    factor_cols=None,
    run_prefix="RR",
    experiment="Experiment",
    confirmation_measurements=None,
    confirmation_run="E1",
    confirmation_coded_setting=None,
    include_screening_rows=True,
):
    """Compare main-effect predictions with measured confirmation values.

    By default the table includes screening-model context rows plus optionally
    one confirmation run such as E1. Set include_screening_rows=False when the
    notebook should show only the confirmation run. Provide
    confirmation_coded_setting when that confirmation run should also get a
    model prediction.
    """
    if factor_cols is None:
        factor_cols = list(levels)

    grand_mean = result_matrix[response].mean()
    main_effects = main_effect_df.set_index("Effect")["Main effect"]

    def predict(coded_setting):
        return predict_response_from_main_effects(
            coded_setting,
            grand_mean,
            main_effects,
            factor_cols,
        )

    all_positive_coded = {factor: 1 for factor in factor_cols}
    main_effect_max_coded = {
        factor: (1 if main_effects[factor] >= 0 else -1)
        for factor in factor_cols
    }

    all_positive_observed = result_matrix.loc[
        (result_matrix[factor_cols] == pd.Series(all_positive_coded)).all(axis=1),
        response,
    ]
    all_positive_observed_value = (
        all_positive_observed.iloc[0] if not all_positive_observed.empty else pd.NA
    )

    observed_max_run = result_matrix[response].idxmax()
    observed_max_coded = (
        result_matrix.loc[observed_max_run, factor_cols]
        .astype(int)
        .to_dict()
    )
    observed_max_value = result_matrix.loc[observed_max_run, response]

    rows = []

    if include_screening_rows:
        rows.extend([
            {
                "case": f"{experiment} all +1 settings",
                "coded_setting": all_positive_coded,
                "real_setting": real_settings_from_coded(
                    all_positive_coded,
                    levels,
                    factor_cols,
                ),
                "predicted_mass_g": predict(all_positive_coded),
                "observed_mass_g": all_positive_observed_value,
                "note": "All factors high is not necessarily the response maximum.",
            },
            {
                "case": f"{experiment} main-effect predicted max",
                "coded_setting": main_effect_max_coded,
                "real_setting": real_settings_from_coded(
                    main_effect_max_coded,
                    levels,
                    factor_cols,
                ),
                "predicted_mass_g": predict(main_effect_max_coded),
                "observed_mass_g": pd.NA,
                "note": "Uses only main effects; confirm because interactions, aliasing, and curvature are not included.",
            },
            {
                "case": f"{experiment} observed max run {run_prefix}{observed_max_run}",
                "coded_setting": observed_max_coded,
                "real_setting": real_settings_from_coded(
                    observed_max_coded,
                    levels,
                    factor_cols,
                ),
                "predicted_mass_g": predict(observed_max_coded),
                "observed_mass_g": observed_max_value,
                "note": "Highest measured response inside the original matrix.",
            },
        ])

    if confirmation_measurements is not None:
        confirmation_observed = confirmation_measurements.loc[
            confirmation_measurements["run"] == confirmation_run,
            response,
        ]
        confirmation_value = (
            confirmation_observed.iloc[0]
            if not confirmation_observed.empty
            else pd.NA
        )
        rows.append({
            "case": f"{confirmation_run} confirmation run",
            "coded_setting": (
                confirmation_coded_setting
                if confirmation_coded_setting is not None
                else pd.NA
            ),
            "real_setting": (
                real_settings_from_coded(
                    confirmation_coded_setting,
                    levels,
                    factor_cols,
                )
                if confirmation_coded_setting is not None
                else pd.NA
            ),
            "predicted_mass_g": (
                predict(confirmation_coded_setting)
                if confirmation_coded_setting is not None
                else pd.NA
            ),
            "observed_mass_g": confirmation_value,
            "note": (
                "Predicted mass uses the FD_II main-effect maximum setting; "
                "compare this measured value against that screening-model prediction."
            ),
        })

    comparison = pd.DataFrame(rows)
    comparison["predicted_mass_g"] = pd.to_numeric(
        comparison["predicted_mass_g"],
        errors="coerce",
    ).round(2)
    comparison["observed_mass_g"] = pd.to_numeric(
        comparison["observed_mass_g"],
        errors="coerce",
    ).round(2)

    return comparison


def coded_setting_from_effects(effect_df, objective):
    """Choose coded factor settings that minimize or maximize a response.

    A positive main effect means the response is higher at +1 than at -1.
    For minimization we choose the opposite sign; for maximization we choose
    the same sign. This is a simple main-effects prediction, so interactions
    and curvature are not included.
    """
    if objective not in {"minimize", "maximize"}:
        raise ValueError("objective must be 'minimize' or 'maximize'")

    signs = effect_df.set_index("Effect")["Main effect"].apply(
        lambda effect: 1 if effect > 0 else -1
    )
    if objective == "minimize":
        signs = -signs

    return signs.astype(int).to_dict()


def build_quadratic_model_matrix(coded_df, factor_cols):
    """Build the quadratic response-surface model matrix for coded factors."""
    model_terms = {"Intercept": np.ones(len(coded_df))}

    for factor in factor_cols:
        model_terms[factor] = coded_df[factor]

    for factor in factor_cols:
        model_terms[f"{factor}^2"] = coded_df[factor] ** 2

    for factor1, factor2 in itertools.combinations(factor_cols, 2):
        model_terms[f"{factor1}:{factor2}"] = coded_df[factor1] * coded_df[factor2]

    return pd.DataFrame(model_terms, index=coded_df.index)


def fit_quadratic_response_surface(matrix, response="mass_mean_g", factor_cols=None):
    """Fit a quadratic response-surface model for a Box-Behnken design.

    Rows without a response are ignored, so the model can be prepared while
    measurements are still incomplete. With exactly as many measured rows as
    coefficients, the fit is interpolating and residual diagnostics are limited.
    """
    if factor_cols is None:
        factor_cols = [col for col in matrix.columns if col not in ["Run", response]]

    analysis_df = matrix.dropna(subset=[response]).copy()
    design_matrix = build_quadratic_model_matrix(analysis_df, factor_cols)
    y = analysis_df[response].to_numpy(dtype=float)

    coefficients, residuals, rank, singular_values = np.linalg.lstsq(
        design_matrix.to_numpy(dtype=float),
        y,
        rcond=None,
    )
    coefficient_series = pd.Series(coefficients, index=design_matrix.columns)
    predicted = design_matrix.to_numpy(dtype=float) @ coefficients
    residual_values = y - predicted
    residual_df = len(analysis_df) - rank

    total_ss = np.sum((y - y.mean()) ** 2)
    residual_ss = np.sum(residual_values ** 2)
    residual_ms = residual_ss / residual_df if residual_df > 0 else np.nan
    r_squared = np.nan if total_ss == 0 else 1 - residual_ss / total_ss

    try:
        xtx_inv = np.linalg.pinv(
            design_matrix.to_numpy(dtype=float).T @ design_matrix.to_numpy(dtype=float)
        )
        standard_errors = np.sqrt(np.diag(xtx_inv) * residual_ms)
    except Exception:
        standard_errors = np.full(len(coefficient_series), np.nan)

    fitted_df = analysis_df[[*factor_cols, response]].copy()
    if "run" in analysis_df.columns:
        fitted_df.insert(0, "run", analysis_df["run"])
    else:
        fitted_df.insert(0, "run", analysis_df.index)
    fitted_df["predicted_response"] = predicted
    fitted_df["residual"] = residual_values

    info_df = pd.DataFrame([{
        "measured_runs": len(analysis_df),
        "model_terms": len(coefficient_series),
        "rank": rank,
        "residual_df": residual_df,
        "r_squared": r_squared,
        "residual_std": np.sqrt(residual_ms) if residual_df > 0 else np.nan,
    }])

    coefficient_df = coefficient_series.rename("coefficient").reset_index()
    coefficient_df = coefficient_df.rename(columns={"index": "term"})
    term_statistics_df = coefficient_df.copy()
    term_statistics_df["standard_error"] = standard_errors
    term_statistics_df["t_value"] = (
        term_statistics_df["coefficient"] / term_statistics_df["standard_error"]
    )
    term_statistics_df["absolute_t_value"] = term_statistics_df["t_value"].abs()
    term_statistics_df["F_value"] = term_statistics_df["t_value"] ** 2
    try:
        from scipy import stats

        term_statistics_df["p_value"] = 2 * stats.t.sf(
            term_statistics_df["absolute_t_value"],
            df=residual_df,
        )
    except Exception:
        term_statistics_df["p_value"] = np.nan

    return {
        "coefficients": coefficient_series,
        "coefficient_df": coefficient_df,
        "term_statistics_df": term_statistics_df,
        "factor_cols": factor_cols,
        "response": response,
        "fitted_df": fitted_df,
        "info_df": info_df,
        "design_matrix": design_matrix,
        "singular_values": singular_values,
    }


def build_quadratic_anova_table(model):
    """Build an overall ANOVA table for a fitted quadratic response model."""
    response = model["response"]
    fitted_df = model["fitted_df"]
    observed = fitted_df[response].to_numpy(dtype=float)
    residual = fitted_df["residual"].to_numpy(dtype=float)

    total_ss = float(np.sum((observed - observed.mean()) ** 2))
    residual_ss = float(np.sum(residual**2))
    model_ss = total_ss - residual_ss

    rank = int(model["info_df"]["rank"].iloc[0])
    model_df = rank - 1
    residual_df = int(model["info_df"]["residual_df"].iloc[0])
    total_df = len(fitted_df) - 1

    model_ms = model_ss / model_df if model_df > 0 else np.nan
    residual_ms = residual_ss / residual_df if residual_df > 0 else np.nan
    f_value = model_ms / residual_ms if residual_ms and not np.isnan(residual_ms) else np.nan
    try:
        from scipy import stats

        p_value = stats.f.sf(f_value, model_df, residual_df)
    except Exception:
        p_value = np.nan

    return pd.DataFrame({
        "source": ["Model", "Residual", "Total"],
        "df": [model_df, residual_df, total_df],
        "sum_sq": [model_ss, residual_ss, total_ss],
        "mean_sq": [model_ms, residual_ms, np.nan],
        "F_value": [f_value, np.nan, np.nan],
        "p_value": [p_value, np.nan, np.nan],
    })


def predict_quadratic_response_surface(model, coded_df):
    """Predict a response from a fitted quadratic response-surface model."""
    design_matrix = build_quadratic_model_matrix(coded_df, model["factor_cols"])
    return design_matrix @ model["coefficients"]


def plot_quadratic_pareto(model, experiment=None, alpha=0.05):
    """Plot a Pareto chart of standardized quadratic-model effects."""
    term_df = model["term_statistics_df"].copy()
    term_df = term_df[term_df["term"] != "Intercept"]
    term_df = term_df.sort_values("absolute_t_value", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [
        FACTOR_COLORS.get(term.split(":")[0].replace("^2", ""), "#7f7f7f")
        for term in term_df["term"]
    ]
    ax.barh(term_df["term"], term_df["absolute_t_value"], color=colors)
    ax.set_xlabel("Absolute standardized effect (|t|)")
    ax.set_title(experiment_title("Pareto chart of quadratic model effects", experiment))
    ax.grid(axis="x", alpha=0.25)

    residual_df = int(model["info_df"]["residual_df"].iloc[0])
    if residual_df > 0:
        try:
            from scipy import stats

            critical_value = stats.t.ppf(1 - alpha / 2, df=residual_df)
            ax.axvline(
                critical_value,
                color="#d95f02",
                linewidth=1.5,
                label=f"alpha={alpha:g} threshold",
            )
            ax.legend(loc="lower right")
        except Exception:
            pass

    plt.tight_layout()
    plt.show(block=False)


def plot_quadratic_coefficients(model, experiment=None):
    """Plot fitted quadratic model coefficients."""
    coefficient_df = model["coefficient_df"].copy()
    coefficient_df = coefficient_df[coefficient_df["term"] != "Intercept"]

    fig, ax = plt.subplots(figsize=(10, 4))
    colors = [
        FACTOR_COLORS.get(term.split(":")[0].replace("^2", ""), "#7f7f7f")
        for term in coefficient_df["term"]
    ]
    ax.bar(coefficient_df["term"], coefficient_df["coefficient"], color=colors)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylabel(f"Coefficient on {model['response']}")
    ax.set_title(experiment_title("Quadratic response-surface coefficients", experiment))
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.show(block=False)


def plot_quadratic_observed_vs_predicted(model, experiment=None):
    """Plot observed versus fitted response values for the quadratic model."""
    fitted_df = model["fitted_df"]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        fitted_df[model["response"]],
        fitted_df["predicted_response"],
        color="#08519c",
    )

    min_value = min(fitted_df[model["response"]].min(), fitted_df["predicted_response"].min())
    max_value = max(fitted_df[model["response"]].max(), fitted_df["predicted_response"].max())
    ax.plot([min_value, max_value], [min_value, max_value], color="#333333", linestyle="--")

    for _, row in fitted_df.iterrows():
        run_label = row["run"]
        if pd.notna(run_label):
            if float(run_label).is_integer():
                run_label = int(run_label)
        ax.annotate(
            run_label,
            xy=(row[model["response"]], row["predicted_response"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(f"Observed {model['response']}")
    ax.set_ylabel(f"Predicted {model['response']}")
    ax.set_title(experiment_title("Observed vs predicted", experiment))
    ax.text(
        1.04,
        0.5,
        "Model fit check:\n"
        "Each point represents one\n"
        "measured Box-Behnken run.\n\n"
        "The dashed line indicates\n"
        "perfect agreement between\n"
        "measured and predicted mass.\n\n"
        "Larger deviations from the\n"
        "line indicate larger model\n"
        "residuals for that run.",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#bbbbbb",
            "alpha": 0.92,
        },
    )
    ax.grid(alpha=0.25)
    fig.subplots_adjust(right=0.68)
    plt.show(block=False)


def plot_quadratic_response_surface_slices(
    model,
    fixed_level=0,
    grid_points=41,
    experiment=None,
):
    """Plot pairwise response-surface contour slices with other factors fixed."""
    factor_cols = model["factor_cols"]
    factor_pairs = list(itertools.combinations(factor_cols, 2))
    fig, axes = plt.subplots(
        1,
        len(factor_pairs),
        figsize=(5 * len(factor_pairs), 4),
        constrained_layout=True,
    )
    if len(factor_pairs) == 1:
        axes = [axes]

    grid = np.linspace(-1, 1, grid_points)

    for ax, (factor_x, factor_y) in zip(axes, factor_pairs):
        xx, yy = np.meshgrid(grid, grid)
        rows = []
        for x_value, y_value in zip(xx.ravel(), yy.ravel()):
            row = {factor: fixed_level for factor in factor_cols}
            row[factor_x] = x_value
            row[factor_y] = y_value
            rows.append(row)

        prediction_df = pd.DataFrame(rows)
        zz = predict_quadratic_response_surface(model, prediction_df).to_numpy()
        zz = zz.reshape(xx.shape)

        contour = ax.contourf(xx, yy, zz, levels=15, cmap="viridis")
        ax.contour(xx, yy, zz, levels=8, colors="white", linewidths=0.5, alpha=0.7)
        ax.set_xlabel(factor_x)
        ax.set_ylabel(factor_y)
        ax.set_title(f"{factor_x} vs {factor_y}")
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])
        fig.colorbar(contour, ax=ax, label=model["response"])

    fig.suptitle(
        experiment_title(
            f"Response-surface slices (other factors fixed at {fixed_level})",
            experiment,
        ),
    )
    plt.show(block=False)


def plot_quadratic_surface_and_contour_pairs(
    model,
    levels=None,
    fixed_level=0,
    grid_points=41,
    experiment=None,
):
    """Plot paired 3D surface and 2D contour maps for each factor pair.

    The model is fitted in coded -1/0/+1 units. If levels are supplied, axes
    are labeled with the real printer settings while predictions still use the
    coded grid internally.
    """
    factor_cols = model["factor_cols"]
    factor_pairs = list(itertools.combinations(factor_cols, 2))
    grid = np.linspace(-1, 1, grid_points)

    def coded_grid_to_real(factor, coded_values):
        coded_levels = np.array([-1, 0, 1], dtype=float)
        real_levels = np.array(
            [levels[factor][coded_level] for coded_level in [-1, 0, 1]],
            dtype=float,
        )
        return np.interp(coded_values, coded_levels, real_levels)

    for row_number, (factor_x, factor_y) in enumerate(factor_pairs, start=1):
        xx_coded, yy_coded = np.meshgrid(grid, grid)
        rows = []
        for x_value, y_value in zip(xx_coded.ravel(), yy_coded.ravel()):
            row = {factor: fixed_level for factor in factor_cols}
            row[factor_x] = x_value
            row[factor_y] = y_value
            rows.append(row)

        prediction_df = pd.DataFrame(rows)
        zz = predict_quadratic_response_surface(model, prediction_df).to_numpy()
        zz = zz.reshape(xx_coded.shape)

        if levels is None:
            xx_plot = xx_coded
            yy_plot = yy_coded
            x_ticks = [-1, 0, 1]
            y_ticks = [-1, 0, 1]
            x_tick_labels = ["-1", "0", "+1"]
            y_tick_labels = ["-1", "0", "+1"]
        else:
            xx_plot = coded_grid_to_real(factor_x, xx_coded)
            yy_plot = coded_grid_to_real(factor_y, yy_coded)
            x_ticks = [levels[factor_x][level] for level in [-1, 0, 1]]
            y_ticks = [levels[factor_y][level] for level in [-1, 0, 1]]
            x_tick_labels = [f"{value:g}" for value in x_ticks]
            y_tick_labels = [f"{value:g}" for value in y_ticks]

        fig = plt.figure(figsize=(15, 5.4))
        ax_surface = fig.add_subplot(1, 2, 1, projection="3d")
        ax_contour = fig.add_subplot(1, 2, 2)

        surface = ax_surface.plot_surface(
            xx_plot,
            yy_plot,
            zz,
            cmap="turbo",
            edgecolor="black",
            linewidth=0.18,
            antialiased=True,
            alpha=0.96,
        )
        z_floor = np.nanmin(zz)
        ax_surface.contour(
            xx_plot,
            yy_plot,
            zz,
            zdir="z",
            offset=z_floor,
            levels=8,
            cmap="turbo",
            linewidths=1,
        )
        ax_surface.set_zlim(z_floor, np.nanmax(zz))
        ax_surface.set_xlabel(factor_x)
        ax_surface.set_ylabel(factor_y)
        ax_surface.set_zlabel(model["response"])
        ax_surface.xaxis.labelpad = 10
        ax_surface.yaxis.labelpad = 10
        ax_surface.zaxis.labelpad = 12
        ax_surface.set_xticks(x_ticks)
        ax_surface.set_yticks(y_ticks)
        ax_surface.set_xticklabels(x_tick_labels)
        ax_surface.set_yticklabels(y_tick_labels)
        ax_surface.view_init(elev=26, azim=-135)

        fitted_df = model["fitted_df"].copy()
        plotted_points = fitted_df[
            fitted_df[factor_cols]
            .drop(columns=[factor_x, factor_y])
            .eq(fixed_level)
            .all(axis=1)
        ]
        if not plotted_points.empty:
            if levels is None:
                point_x = plotted_points[factor_x]
                point_y = plotted_points[factor_y]
            else:
                point_x = plotted_points[factor_x].map(levels[factor_x])
                point_y = plotted_points[factor_y].map(levels[factor_y])
            ax_surface.scatter(
                point_x,
                point_y,
                plotted_points[model["response"]],
                color="#de2d26",
                edgecolor="black",
                s=45,
                depthshade=False,
            )

        filled = ax_contour.contourf(
            xx_plot,
            yy_plot,
            zz,
            levels=30,
            cmap="turbo",
        )
        contour_lines = ax_contour.contour(
            xx_plot,
            yy_plot,
            zz,
            levels=8,
            colors="black",
            linewidths=0.9,
        )
        ax_contour.clabel(contour_lines, inline=True, fontsize=8, fmt="%.1f")
        if not plotted_points.empty:
            ax_contour.scatter(
                point_x,
                point_y,
                color="#de2d26",
                edgecolor="black",
                s=35,
                zorder=3,
            )
        ax_contour.set_xlabel(factor_x)
        ax_contour.set_ylabel(factor_y)
        ax_contour.set_xticks(x_ticks)
        ax_contour.set_yticks(y_ticks)
        ax_contour.set_xticklabels(x_tick_labels)
        ax_contour.set_yticklabels(y_tick_labels)
        fig.colorbar(filled, ax=ax_contour, label=model["response"])

        fig.suptitle(
            experiment_title(
                f"Surface and contour {row_number}: {factor_x} vs {factor_y}",
                experiment,
            ),
        )
        fig.subplots_adjust(
            left=0.08,
            right=0.93,
            bottom=0.13,
            top=0.88,
            wspace=0.28,
        )
        plt.show(block=False)


def real_value_from_coded(factor, coded_value, levels):
    """Convert one coded factor value to a real printer setting.

    -1 and +1 use the low/high values from levels. A coded value of 0 means
    the midpoint between low and high. Custom values, such as flow ratio 1.10,
    are already real settings and pass through unchanged.
    """
    if coded_value in levels[factor]:
        return levels[factor][coded_value]
    if coded_value == 0:
        return sum(levels[factor].values()) / len(levels[factor])
    return coded_value


def candidate_real_table(candidate_coded_df, levels):
    """Convert a candidate run table from coded values to real settings."""
    real_df = candidate_coded_df.copy()
    for factor in levels:
        real_df[factor] = real_df[factor].apply(
            lambda coded_value: real_value_from_coded(factor, coded_value, levels)
        )
    return real_df


def build_extra_run_candidates(
    min_response_coded,
    max_response_coded,
    levels,
    extra_flow_ratios=(1.08, 1.10),
    middle_factors=None,
):
    """Build named extra-run candidates for confirmation and curvature checks."""
    if middle_factors is None:
        middle_factors = list(levels)

    candidate_rows = [
        {
            "Extra_run_name": "Density_max",
            "Comment": "Confirm predicted density maximum (minimize mass deviation)",
            **min_response_coded,
        },
        {
            "Extra_run_name": "Density_min",
            "Comment": "Confirm predicted density minimum (maximize mass deviation)",
            **max_response_coded,
        },
        {
            "Extra_run_name": "All_middle_values",
            "Comment": "All factors at middle values; single center-point check",
            **{factor: 0 for factor in levels},
        },
    ]

    for factor in middle_factors:
        row = {
            "Extra_run_name": f"Middle_{factor}_at_density_max",
            "Comment": f"Density-max settings with {factor} at middle value",
            **min_response_coded,
        }
        row[factor] = 0
        candidate_rows.append(row)

    for flow_ratio in extra_flow_ratios:
        row = {
            "Extra_run_name": f"Flow_{flow_ratio:g}",
            "Comment": "Higher-flow exploratory run using low mass deviation settings",
            **min_response_coded,
        }
        row["A_flow_ratio"] = flow_ratio
        candidate_rows.append(row)

    candidate_coded_df = pd.DataFrame(candidate_rows).set_index("Extra_run_name")
    candidate_real_df = candidate_real_table(
        candidate_coded_df.drop(columns="Comment"),
        levels,
    )
    candidate_real_df.insert(0, "Comment", candidate_coded_df["Comment"])

    return candidate_coded_df, candidate_real_df


def style_extra_run_candidates(candidate_df, comment_column="Comment"):
    """Style candidate tables with compact numbers and a wider comment column."""
    comment_style = {
        "white-space": "normal",
        "min-width": "430px",
        "max-width": "620px",
    }
    number_format = {
        column: "{:g}"
        for column in candidate_df.columns
        if column != comment_column
    }

    return candidate_df.style.format(number_format).set_properties(
        subset=[comment_column],
        **comment_style,
    )


def experiment_title(title, experiment=None):
    return f"{experiment} - {title}" if experiment else title


def factor_code(factor):
    return factor.split("_", 1)[0]

def plot_mass_measurements(
    mass_measurements_g,
    mass_columns,
    y_min=None,
    y_max=None,
    mean_label="Mean",
    experiment=None,
    repeat_measurements_g=None,
    reference_mass_g=None,
    reference_label="Reference",
    reference_color="#f28e2b",
    reference_details=None,
    value_digits=1,
    repeat_bar_offset=None,
    extra_measurements_g=None,
    extra_label="extra mean",
    extra_color="#756bb1",
):
    plot_columns = mass_columns + ["mass_mean_g"]
    if "mass_std_g" in mass_measurements_g.columns:
        plot_columns.append("mass_std_g")

    plot_df = mass_measurements_g.set_index("run")[plot_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )
    run_labels = plot_df.index.tolist()
    x = np.arange(len(plot_df))
    bar_width = 0.22
    summary_lines = []
    repeat_original_runs_with_data = set()

    fig, ax = plt.subplots(figsize=(16, 4.3))
    offsets = np.linspace(
        -bar_width * (len(mass_columns) - 1) / 2,
        bar_width * (len(mass_columns) - 1) / 2,
        len(mass_columns),
    )
    colors = ["#9ecae1", "#4292c6", "#08519c", "#08306b"]

    for offset, column, color in zip(offsets, mass_columns, colors):
        ax.bar(x + offset, plot_df[column], width=bar_width, label=column, color=color)

    tick_positions = x.tolist()
    tick_labels = run_labels.copy()
    if reference_mass_g is not None and pd.notna(reference_mass_g):
        reference_x = -1
        ax.axhline(
            reference_mass_g,
            color=reference_color,
            linestyle="-",
            linewidth=1.5,
            alpha=0.35,
            zorder=0,
            label=f"{reference_label} level",
        )
        ax.bar(
            reference_x,
            reference_mass_g,
            width=bar_width * 1.4,
            color=reference_color,
            label=reference_label,
            alpha=0.9,
        )
        ax.annotate(
            f"{reference_mass_g:.{value_digits}f}",
            xy=(reference_x, reference_mass_g),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#7a3b00",
        )
        tick_positions = [reference_x, *tick_positions]
        tick_labels = [reference_label, *tick_labels]

    if reference_details:
        ax.text(
            1.02,
            0.5,
            reference_details,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=9,
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "white",
                "edgecolor": reference_color,
                "alpha": 0.9,
            },
        )

    if repeat_measurements_g is not None and not repeat_measurements_g.empty:
        repeat_df = repeat_measurements_g.copy()
        repeat_df["original_run"] = repeat_df["run"].str.split(".").str[0]
        run_position = {run: position for position, run in enumerate(run_labels)}
        repeat_df["x"] = repeat_df["original_run"].map(run_position)
        repeat_df = repeat_df.dropna(subset=["x", "mass_mean_g"])

        if not repeat_df.empty:
            repeat_df = (
                repeat_df.groupby("original_run", as_index=False)
                .agg({
                    "x": "first",
                    "mass_mean_g": "mean",
                })
            )
            repeat_original_runs_with_data = set(repeat_df["original_run"])
            repeat_x = repeat_df["x"].to_numpy(dtype=float)
            if repeat_bar_offset is None:
                repeat_bar_offset = (len(mass_columns) + 2.4) * bar_width / 2
            repeat_x += repeat_bar_offset
            ax.bar(
                repeat_x,
                repeat_df["mass_mean_g"],
                width=bar_width,
                color="#31a354",
                label="repeat mean",
                alpha=0.85,
            )

            for xi, run, mass in zip(
                repeat_x,
                repeat_df["original_run"],
                repeat_df["mass_mean_g"],
            ):
                ax.annotate(
                    f"{mass:.{value_digits}f}",
                    xy=(xi, mass),
                    xytext=(0, 10),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="black",
                )

            repeatability_df = calculate_repeatability(
                mass_measurements_g,
                repeat_measurements_g,
            )
            if not repeatability_df.empty:
                summary_lines.append(
                    "Repeatability: "
                    f"mean abs diff {repeatability_df['abs_difference_g'].mean():.3f} g; "
                    f"max {repeatability_df['abs_difference_g'].max():.3f} g"
                )

    if extra_measurements_g is not None and not extra_measurements_g.empty:
        extra_df = extra_measurements_g[["run", "mass_mean_g"]].copy()
        extra_df["mass_mean_g"] = pd.to_numeric(
            extra_df["mass_mean_g"],
            errors="coerce",
        )
        extra_df = extra_df.dropna(subset=["mass_mean_g"])

        if not extra_df.empty:
            extra_x = np.arange(
                len(plot_df) + 1,
                len(plot_df) + 1 + len(extra_df),
                dtype=float,
            )
            ax.bar(
                extra_x,
                extra_df["mass_mean_g"],
                width=bar_width,
                color=extra_color,
                label=extra_label,
                alpha=0.88,
            )
            for xi, run, mass in zip(
                extra_x,
                extra_df["run"],
                extra_df["mass_mean_g"],
            ):
                ax.annotate(
                    f"{mass:.{value_digits}f}",
                    xy=(xi, mass),
                    xytext=(0, 10),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color=extra_color,
                )
            tick_positions.extend(extra_x.tolist())
            tick_labels.extend(extra_df["run"].tolist())

    ax.plot(
        x,
        plot_df["mass_mean_g"],
        color="#000000",
        marker="o",
        linewidth=2,
        label=mean_label,
    )

    for xi, run, mass in zip(x, plot_df.index, plot_df["mass_mean_g"]):
        if pd.notna(mass) and run not in repeat_original_runs_with_data:
            ax.annotate(
                f"{mass:.{value_digits}f}",
                xy=(xi, mass),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_xlabel("Test plate")
    ax.set_ylabel("Mass (g)")
    ax.set_title(experiment_title("Mass measurements", experiment))
    if "mass_std_g" in plot_df.columns and plot_df["mass_std_g"].notna().any():
        measurement_uncertainty_g = np.sqrt((plot_df["mass_std_g"] ** 2).mean())
        ax.text(
            0.99,
            0.96,
            f"Mean measurement uncertainty: ±{measurement_uncertainty_g:.3f} g",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#bbbbbb",
                "alpha": 0.85,
            },
        )
    if summary_lines:
        ax.text(
            0.99,
            0.04,
            "\n".join(summary_lines),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=9,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#31a354",
                "alpha": 0.85,
            },
        )
    if y_min is not None or y_max is not None:
        ax.set_ylim(bottom=y_min, top=y_max)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    ax.grid(axis="y", alpha=0.3)

    if reference_details:
        fig.subplots_adjust(right=0.70)
    else:
        fig.subplots_adjust(right=0.82, bottom=0.20)
    plt.show(block=False)


def plot_combined_mass_measurement_sets(
    measurement_sets,
    minimum_mass=80,
    experiment=None,
):
    """Compare mean masses from multiple run sets on one shared plot."""
    colors = ["#2b8cbe", "#31a354", "#756bb1", "#e6550d"]
    fig, ax = plt.subplots(figsize=(15, 5))
    tick_positions = []
    tick_labels = []
    current_x = 0
    plotted_values = []
    nonempty_sets = []

    for set_label, measurements in measurement_sets.items():
        plot_df = measurements.copy()
        plot_df["mass_mean_g"] = pd.to_numeric(
            plot_df["mass_mean_g"],
            errors="coerce",
        )
        plot_df = plot_df.loc[plot_df["mass_mean_g"] > minimum_mass].copy()
        if not plot_df.empty:
            nonempty_sets.append((set_label, plot_df))

    for set_index, (set_label, plot_df) in enumerate(nonempty_sets):
        x = np.arange(current_x, current_x + len(plot_df))
        color = colors[set_index % len(colors)]
        y_error = None
        if "mass_std_g" in plot_df.columns:
            std_values = pd.to_numeric(plot_df["mass_std_g"], errors="coerce")
            if std_values.notna().any():
                y_error = std_values.fillna(0)

        ax.errorbar(
            x,
            plot_df["mass_mean_g"],
            yerr=y_error,
            color=color,
            marker="o",
            markersize=5,
            linewidth=1.5,
            capsize=3,
            label=set_label,
        )
        tick_positions.extend(x.tolist())
        tick_labels.extend(plot_df["run"].tolist())
        plotted_values.extend(plot_df["mass_mean_g"].tolist())

        current_x += len(plot_df)
        if set_index < len(nonempty_sets) - 1:
            ax.axvline(current_x - 0.5, color="#777777", linewidth=1, alpha=0.5)
            current_x += 1

    if plotted_values:
        value_range = max(plotted_values) - min(plotted_values)
        margin = max(0.25, value_range * 0.08)
        ax.set_ylim(min(plotted_values) - margin, max(plotted_values) + margin)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax.set_xlabel("Run")
    ax.set_ylabel("Mean mass (g)")
    ax.set_title(experiment_title("Combined mass measurements", experiment))
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_ideal_mass_deviation(results_df, experiment=None):
    plot_results = results_df.dropna(
        subset=["measured_to_ideal_mass_percent", "mass_deviation_from_ideal_percent"]
    )
    x = np.arange(len(plot_results))

    fig, ax = plt.subplots(figsize=(8, 4))
    deviation_colors = [
        "#9ecae1" if value >= 0 else "#de2d26"
        for value in plot_results["mass_deviation_from_ideal_percent"]
    ]
    ax.bar(x, plot_results["mass_deviation_from_ideal_percent"], color=deviation_colors)

    for xi, deviation in zip(x, plot_results["mass_deviation_from_ideal_percent"]):
        ax.annotate(
            f"{deviation:.1f}%",
            xy=(xi, deviation),
            xytext=(0, -6 if deviation >= 0 else 6),
            textcoords="offset points",
            ha="center",
            va="top" if deviation >= 0 else "bottom",
            fontsize=9,
            color="#333333",
        )

    ax.axhline(0, color="#333333", linewidth=1)
    ax.set_xlabel("Test plate")
    ax.set_ylabel("Mass deviation from ideal (%)")
    ax.set_title(experiment_title("Mass deviation from ideal mass", experiment))
    ax.set_xticks(x)
    ax.set_xticklabels(plot_results["run"], rotation=45)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.show(block=False)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(
        x,
        plot_results["measured_to_ideal_mass_percent"],
        color="#08519c",
        marker="o",
        linewidth=2,
    )
    ax.axhline(100, color="#333333", linestyle="-", linewidth=1.3, alpha=0.75, label="Ideal mass (100%)")
    ax.axhline(105, color="#de2d26", linestyle="--", linewidth=1.5, alpha=0.8, label="Flow rate 1.05 (105%)")
    ax.set_xlabel("Test plate")
    ax.set_ylabel("Measured / ideal mass (%)")
    ax.set_title(experiment_title("Measured mass relative to ideal mass", experiment))
    ax.set_xticks(x)
    ax.set_xticklabels(plot_results["run"], rotation=45)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.show(block=False)


def plot_porosity_and_relative_density(results_df, experiment=None):
    """Deprecated name; use plot_ideal_mass_deviation."""
    plot_ideal_mass_deviation(results_df, experiment=experiment)


def plot_main_effects_bar(main_effect_df, response="mass_mean_g", experiment=None):
    fig, ax = plt.subplots(figsize=(8, 3))
    colors = [
        FACTOR_COLORS.get(effect, "#7f7f7f")
        for effect in main_effect_df["Effect"]
    ]
    ax.bar(main_effect_df["Effect"], main_effect_df["Main effect"], color=colors)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylabel(f"Main effect on {response}")
    ax.set_title(experiment_title("Main effects", experiment))
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.show(block=False)


def plot_main_effects_line(matrix, response="mass_mean_g", experiment=None):
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

        ax.plot(
            [x1, x2],
            [mean_low, mean_high],
            marker="o",
            color=FACTOR_COLORS.get(factor, "#7f7f7f"),
        )
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
    ax.set_title(experiment_title(f"Main effects plot ({response})", experiment))
    fig.subplots_adjust(bottom=0.18)
    plt.show(block=False)


def plot_main_effects_line_comparison(
    matrices,
    labels,
    response="mass_mean_g",
    experiment=None,
):
    fig, axes = plt.subplots(1, len(matrices), figsize=(6 * len(matrices), 4), sharey=True)
    if len(matrices) == 1:
        axes = [axes]

    for ax, matrix, label in zip(axes, matrices, labels):
        factor_cols = [col for col in matrix.columns if col not in ["Run", response]]
        overall_mean = matrix[response].mean()
        xticks = []
        xtick_labels = []

        for i, factor in enumerate(factor_cols):
            mean_low = matrix.loc[matrix[factor] == -1, response].mean()
            mean_high = matrix.loc[matrix[factor] == 1, response].mean()
            x1 = i * 3 + 1
            x2 = i * 3 + 2

            ax.plot(
                [x1, x2],
                [mean_low, mean_high],
                marker="o",
                color=FACTOR_COLORS.get(factor, "#7f7f7f"),
            )
            xticks.extend([x1, x2])
            xtick_labels.extend(["-", "+"])
            ax.text(
                (x1 + x2) / 2,
                -0.14,
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
        ax.set_title(label)
        ax.set_xlabel("Factors", labelpad=26)
        ax.grid(axis="y", alpha=0.2)

    axes[0].set_ylabel("Average Response")
    fig.suptitle(experiment_title(f"Conditional trends ({response})", experiment))
    fig.subplots_adjust(bottom=0.22, top=0.84, wspace=0.18)
    plt.show(block=False)


def plot_three_level_main_effects(
    matrix,
    response="mass_mean_g",
    factor_cols=None,
    levels=None,
    experiment=None,
):
    """Plot mean response at -1, 0, and +1 for each three-level factor."""
    if factor_cols is None:
        factor_cols = [col for col in matrix.columns if col not in ["Run", response]]

    plot_matrix = matrix.dropna(subset=[response])
    fig, axes = plt.subplots(
        1,
        len(factor_cols),
        figsize=(4.3 * len(factor_cols), 3.6),
        sharey=True,
    )
    if len(factor_cols) == 1:
        axes = [axes]

    for ax, factor in zip(axes, factor_cols):
        coded_levels = [-1, 0, 1]
        means = [
            plot_matrix.loc[plot_matrix[factor] == coded_level, response].mean()
            for coded_level in coded_levels
        ]
        if levels is None:
            x_values = coded_levels
            x_labels = ["-1", "0", "+1"]
        else:
            x_values = [levels[factor][coded_level] for coded_level in coded_levels]
            x_labels = [f"{value:g}" for value in x_values]

        ax.plot(
            x_values,
            means,
            marker="o",
            color=FACTOR_COLORS.get(factor, "#008b8b"),
            linewidth=1.8,
        )
        ax.set_title(factor)
        ax.set_xlabel(factor)
        ax.set_xticks(x_values)
        ax.set_xticklabels(x_labels)
        ax.grid(axis="y", alpha=0.25)

    axes[0].set_ylabel(response)
    fig.suptitle(experiment_title(f"Three-level main effects ({response})", experiment))
    fig.subplots_adjust(bottom=0.22, top=0.82, wspace=0.18)
    plt.show(block=False)


def add_interactions(matrix, response="mass_mean_g"):
    matrix_interactions = matrix.copy()
    
    factor_cols = [col for col in matrix.columns if col not in ["Run", response]]
    
    for factor1, factor2 in itertools.combinations(factor_cols, 2):
        interaction_name = f"{factor1} × {factor2}"
        matrix_interactions[interaction_name] = matrix[factor1] * matrix[factor2]
    
    return matrix_interactions

def add_interactions(matrix, response="mass_mean_g"):
    matrix_interactions = matrix.copy()
    factor_cols = [col for col in matrix.columns if col not in ["Run", response]]

    for factor1, factor2 in itertools.combinations(factor_cols, 2):
        interaction_name = f"{factor1} x {factor2}"
        matrix_interactions[interaction_name] = matrix[factor1] * matrix[factor2]

    return matrix_interactions


def interaction_alias_note(factor1, factor2, generators=None):
    if not generators:
        return ""

    pair = {factor1, factor2}
    for generated_factor, source_factors in generators.items():
        if pair == set(source_factors):
            generated_code = factor_code(generated_factor)
            source_code = "".join(factor_code(factor) for factor in source_factors)
            return f"{generated_factor}; generator {generated_code}={source_code}"

    return ""


def compute_interaction_effects(matrix, response="mass_mean_g", generators=None):
    factor_cols = [col for col in matrix.columns if col not in ["Run", response]]
    interaction_matrix = add_interactions(matrix, response=response)
    rows = []

    for factor1, factor2 in itertools.combinations(factor_cols, 2):
        interaction_name = f"{factor1} x {factor2}"
        rows.append({
            "Interaction": interaction_name,
            "Alias / generator note": interaction_alias_note(
                factor1,
                factor2,
                generators=generators,
            ),
            "Effect": compute_main_effect(
                interaction_matrix,
                interaction_name,
                response=response,
            ),
        })

    return pd.DataFrame(rows).sort_values(
        by="Effect",
        key=abs,
        ascending=False,
    )


def plot_interaction_effects(interaction_effects, experiment=None):
    if isinstance(interaction_effects, pd.DataFrame):
        interaction_effects_df = interaction_effects.copy()
    else:
        interaction_effects_df = pd.DataFrame({
            "Interaction": interaction_effects.keys(),
            "Effect": interaction_effects.values(),
        })

    interaction_effects_df = interaction_effects_df.sort_values(
        by="Effect",
        key=abs,
        ascending=False,
    )

    interaction_effects_df["Color factor"] = (
        interaction_effects_df["Interaction"].str.split(" x ").str[0]
    )
    colors = [
        FACTOR_COLORS.get(factor, "#7f7f7f")
        for factor in interaction_effects_df["Color factor"]
    ]

    plt.figure(figsize=(16, 5))
    plt.bar(
        interaction_effects_df["Interaction"],
        interaction_effects_df["Effect"],
        color=colors,
    )
    plt.axhline(0, color="black", linewidth=1)

    plt.ylabel("Effect")
    plt.title(experiment_title("Interaction effects plot", experiment))
    plt.xticks(rotation=35, ha="right")

    legend_factors = list(dict.fromkeys(interaction_effects_df["Color factor"]))
    legend_handles = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            color=FACTOR_COLORS.get(factor, "#7f7f7f"),
            label=factor,
        )
        for factor in legend_factors
    ]
    plt.legend(
        handles=legend_handles,
        title="First factor",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
    )

    plt.tight_layout(rect=[0, 0, 0.84, 1])
    plt.show(block=False)
