"""Analysis helpers for the copper definitive screening design."""

from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


FACTOR_COLORS = {
    "A_layer_height_mm": "#1f77b4",
    "B_flow_ratio": "#ff7f0e",
    "C_nozzle_temp_C": "#2ca02c",
    "D_speed_mm_s": "#d62728",
    "E_infill_percent": "#9467bd",
}


FACTOR_SHORT_NAMES = {
    "A_layer_height_mm": "A",
    "B_flow_ratio": "B",
    "C_nozzle_temp_C": "C",
    "D_speed_mm_s": "D",
    "E_infill_percent": "E",
}


FACTOR_LABELS = {
    "A_layer_height_mm": "Layer height (mm)",
    "B_flow_ratio": "Flow ratio (-)",
    "C_nozzle_temp_C": r"Nozzle temperature ($^\circ$C)",
    "D_speed_mm_s": "Print speed (mm/s)",
    "E_infill_percent": "Infill (%)",
}


INTERACTION_FACTOR_LABELS = {
    "A_layer_height_mm": "A_layer_height (mm)",
    "B_flow_ratio": "B_flow_ratio",
    "C_nozzle_temp_C": "C_nozzle_temp (\N{DEGREE SIGN}C)",
    "D_speed_mm_s": "D_speed (mm/s)",
    "E_infill_percent": "E_infill (%)",
}


RESPONSE_LABELS = {
    "weight_g": "Mass (g)",
    "mass_g": "Mass (g)",
    "weight_in_air_g": "Mass in air (g)",
    "weight_in_water_g": "Mass in water (g)",
    "volume_cm3": "Volume (cm3)",
    "density_g_cm3": "Density (g/cm3)",
    "relative_density_percent": "Reported Archimedes percentage (%)",
    "porosity_percent": "Porosity (%)",
    "geometric_density_g_cm3": "Geometric density (g/cm3)",
    "avg_density_g_cm3": "Average density (g/cm3)",
    "avg_relative_density_percent": "Average reported percentage (%)",
    "sintered_density_g_cm3": "Sintered density (g/cm3)",
}


RESPONSE_UNITS = {
    "weight_g": "g",
    "mass_g": "g",
    "weight_in_air_g": "g",
    "weight_in_water_g": "g",
    "volume_cm3": "cm3",
    "density_g_cm3": "g/cm3",
    "relative_density_percent": "%",
    "porosity_percent": "%",
    "geometric_density_g_cm3": "g/cm3",
    "avg_density_g_cm3": "g/cm3",
    "avg_relative_density_percent": "%",
    "sintered_density_g_cm3": "g/cm3",
}


def response_label(response):
    """Human-readable response label for plots."""
    return RESPONSE_LABELS.get(response, response)


def response_unit(response):
    """Unit string for a response, if known."""
    return RESPONSE_UNITS.get(response, "")


def residual_label(response):
    """Residual axis label with response unit when available."""
    unit = response_unit(response)
    return f"Residual ({unit})" if unit else "Residual"


def factor_color_for_term(term):
    """Use the owning factor color for main, quadratic, and interaction terms."""
    base_factor = term.split(":")[0].replace("^2", "")
    return FACTOR_COLORS.get(base_factor, "#777777")


def short_label_for_term(term):
    """Make compact plot labels such as A, B^2, or A:B."""
    if ":" in term:
        return ":".join(short_label_for_term(part) for part in term.split(":"))
    if term.endswith("^2"):
        return f"{short_label_for_term(term[:-2])}^2"
    return FACTOR_SHORT_NAMES.get(term, term)


def default_dsd_terms(factor_names, include_quadratic=False, interactions=None):
    """Build a reduced DSD model term list.

    Start with main effects. Add quadratic terms once all 13 runs are measured,
    and add only selected physically meaningful interactions.
    """
    terms = list(factor_names)
    if include_quadratic:
        terms.extend(f"{factor}^2" for factor in factor_names)
    if interactions:
        terms.extend(f"{factor1}:{factor2}" for factor1, factor2 in interactions)
    return terms


def build_copper_model_matrix(coded_df, terms):
    """Create the regression matrix from coded -1/0/+1 factor settings."""
    model_matrix = pd.DataFrame({"Intercept": np.ones(len(coded_df))}, index=coded_df.index)
    for term in terms:
        if term.endswith("^2"):
            factor = term[:-2]
            model_matrix[term] = coded_df[factor] ** 2
        elif ":" in term:
            factor1, factor2 = term.split(":")
            model_matrix[term] = coded_df[factor1] * coded_df[factor2]
        else:
            model_matrix[term] = coded_df[term]
    return model_matrix


def fit_copper_response_model(coded_with_response_df, response, terms):
    """Fit a reduced linear/quadratic model for one copper response."""
    analysis_df = coded_with_response_df.dropna(subset=[response]).copy()
    if analysis_df.empty:
        raise ValueError(f"No measured rows are available for {response}.")

    model_matrix = build_copper_model_matrix(analysis_df, terms)
    x = model_matrix.to_numpy(dtype=float)
    y = analysis_df[response].to_numpy(dtype=float)

    coefficients, _, rank, singular_values = np.linalg.lstsq(x, y, rcond=None)
    coefficient_series = pd.Series(coefficients, index=model_matrix.columns)
    predicted = x @ coefficients
    residual = y - predicted
    residual_df = len(analysis_df) - rank
    residual_ss = float(np.sum(residual**2))
    residual_ms = residual_ss / residual_df if residual_df > 0 else np.nan

    standard_errors = np.full(len(coefficient_series), np.nan)
    if residual_df > 0 and rank == model_matrix.shape[1]:
        xtx_inv = np.linalg.pinv(x.T @ x)
        standard_errors = np.sqrt(np.diag(xtx_inv) * residual_ms)

    factor_names = []
    for term in terms:
        if term.endswith("^2"):
            factor = term[:-2]
            if factor not in factor_names:
                factor_names.append(factor)
        elif ":" in term:
            for factor in term.split(":"):
                if factor not in factor_names:
                    factor_names.append(factor)
        elif term not in factor_names:
            factor_names.append(term)

    coefficient_df = coefficient_series.rename("coefficient").reset_index()
    coefficient_df = coefficient_df.rename(columns={"index": "term"})
    coefficient_df["standard_error"] = standard_errors
    coefficient_df["standardized_effect"] = (
        coefficient_df["coefficient"] / coefficient_df["standard_error"]
    )
    coefficient_df["t_value"] = coefficient_df["standardized_effect"]
    coefficient_df["F_value"] = coefficient_df["t_value"] ** 2
    try:
        from scipy import stats

        coefficient_df["p_value"] = (
            2 * stats.t.sf(coefficient_df["t_value"].abs(), df=residual_df)
            if residual_df > 0
            else np.nan
        )
    except Exception:
        coefficient_df["p_value"] = np.nan
    coefficient_df["absolute_standardized_effect"] = (
        coefficient_df["standardized_effect"].abs()
    )
    coefficient_df["absolute_coefficient"] = coefficient_df["coefficient"].abs()

    fitted_df = analysis_df[["Run", *factor_names, response]].copy()
    fitted_df["predicted_response"] = predicted
    fitted_df["residual"] = residual
    fitted_df["print_order"] = np.arange(1, len(fitted_df) + 1)

    return {
        "response": response,
        "terms": terms,
        "factor_names": factor_names,
        "coefficients": coefficient_series,
        "coefficient_df": coefficient_df,
        "fitted_df": fitted_df,
        "rank": rank,
        "residual_df": residual_df,
        "singular_values": singular_values,
        "full_rank": rank == model_matrix.shape[1],
        "model_matrix_columns": list(model_matrix.columns),
    }


def build_copper_term_significance_table(model, include_intercept=False):
    """Return term-wise coefficient significance statistics for a copper model."""
    term_df = model["coefficient_df"].copy()
    if not include_intercept:
        term_df = term_df[term_df["term"] != "Intercept"]
    return term_df.reindex(columns=[
        "term",
        "coefficient",
        "standard_error",
        "t_value",
        "F_value",
        "p_value",
    ])


def build_copper_term_anova_table(model):
    """Build partial term-wise ANOVA tests by dropping one term at a time.

    Each term sum of squares is the increase in residual sum of squares when
    that term is removed while all other terms remain in the model. These are
    Type III-style partial tests and therefore need not add to the model sum of
    squares for a non-orthogonal design.
    """
    response = model["response"]
    fitted_df = model["fitted_df"]
    terms = list(model["terms"])
    y = fitted_df[response].to_numpy(dtype=float)
    full_residual = fitted_df["residual"].to_numpy(dtype=float)
    full_residual_ss = float(np.sum(full_residual**2))
    full_rank = int(model["rank"])
    residual_df = int(model["residual_df"])
    residual_ms = (
        full_residual_ss / residual_df
        if residual_df > 0
        else np.nan
    )

    rows = []
    for term in terms:
        reduced_terms = [candidate for candidate in terms if candidate != term]
        reduced_matrix = build_copper_model_matrix(fitted_df, reduced_terms)
        reduced_x = reduced_matrix.to_numpy(dtype=float)
        reduced_coefficients, _, reduced_rank, _ = np.linalg.lstsq(
            reduced_x,
            y,
            rcond=None,
        )
        reduced_residual = y - reduced_x @ reduced_coefficients
        reduced_residual_ss = float(np.sum(reduced_residual**2))

        term_df = full_rank - int(reduced_rank)
        partial_ss = max(reduced_residual_ss - full_residual_ss, 0.0)
        mean_sq = partial_ss / term_df if term_df > 0 else np.nan
        f_value = (
            mean_sq / residual_ms
            if term_df > 0 and pd.notna(residual_ms) and residual_ms > 0
            else np.nan
        )
        try:
            from scipy import stats

            p_value = (
                stats.f.sf(f_value, term_df, residual_df)
                if pd.notna(f_value) and residual_df > 0
                else np.nan
            )
        except Exception:
            p_value = np.nan

        rows.append({
            "source": term,
            "df": term_df,
            "sum_sq": partial_ss,
            "mean_sq": mean_sq,
            "F_value": f_value,
            "p_value": p_value,
        })

    rows.append({
        "source": "Residual",
        "df": residual_df,
        "sum_sq": full_residual_ss,
        "mean_sq": residual_ms,
        "F_value": np.nan,
        "p_value": np.nan,
    })
    return pd.DataFrame(rows)


def build_copper_anova_table(model):
    """Build an overall ANOVA table for a fitted copper regression model."""
    response = model["response"]
    fitted_df = model["fitted_df"]
    observed = fitted_df[response].to_numpy(dtype=float)
    residual = fitted_df["residual"].to_numpy(dtype=float)

    total_ss = float(np.sum((observed - observed.mean()) ** 2))
    residual_ss = float(np.sum(residual**2))
    model_ss = total_ss - residual_ss

    model_df = int(model["rank"]) - 1
    residual_df = int(model["residual_df"])
    total_df = len(fitted_df) - 1

    model_ms = model_ss / model_df if model_df > 0 else np.nan
    residual_ms = residual_ss / residual_df if residual_df > 0 else np.nan
    f_value = (
        model_ms / residual_ms
        if residual_ms and not np.isnan(residual_ms)
        else np.nan
    )
    try:
        from scipy import stats

        p_value = (
            stats.f.sf(f_value, model_df, residual_df)
            if pd.notna(f_value) and residual_df > 0
            else np.nan
        )
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


def plot_copper_anova(
    anova_table,
    response="weight_g",
    title=None,
    variation_percent=False,
):
    """Plot explained and residual variation for a copper model."""
    plot_df = anova_table[anova_table["source"].isin(["Model", "Residual"])].copy()
    if plot_df.empty:
        return

    if variation_percent:
        total_sum_sq = plot_df["sum_sq"].sum()
        if not np.isfinite(total_sum_sq) or total_sum_sq <= 0:
            return
        plot_df["plot_value"] = 100.0 * plot_df["sum_sq"] / total_sum_sq
        plot_df["plot_label"] = plot_df["source"].map({
            "Model": "Explained by model",
            "Residual": "Residual",
        })
        colors = [
            "#4f9275" if source == "Model" else "#e07a3f"
            for source in plot_df["source"]
        ]
        value_format = lambda value: f"{value:.1f}%"
        y_label = "Share of total density variation (%)"
    else:
        plot_df["plot_value"] = plot_df["sum_sq"]
        plot_df["plot_label"] = plot_df["source"]
        colors = [
            "#4c78a8" if source == "Model" else "#bbbbbb"
            for source in plot_df["source"]
        ]
        value_format = lambda value: f"{value:.3g}"
        unit = response_unit(response)
        y_label = f"Sum of squares ({unit}^2)" if unit else "Sum of squares"

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bars = ax.bar(plot_df["plot_label"], plot_df["plot_value"], color=colors)
    for bar, value in zip(bars, plot_df["plot_value"]):
        ax.annotate(
            value_format(value),
            (bar.get_x() + bar.get_width() / 2, value),
            xytext=(0, -8),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=14,
            fontweight="bold",
            color=bar.get_facecolor(),
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.78,
                "pad": 0.25,
            },
            zorder=4,
        )

    ax.set_ylabel(y_label, fontsize=14)
    ax.set_title(
        title or f"Copper DSD - ANOVA for {response_label(response)}",
        pad=14,
        fontsize=17,
    )
    ax.tick_params(axis="both", labelsize=13)
    if variation_percent:
        ax.set_ylim(0, max(105, plot_df["plot_value"].max() * 1.15))
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_pareto(model, title=None, alpha=0.05, secondary_alpha=None):
    """Pareto chart of standardized effects, falling back to coefficients."""
    effect_df = model["coefficient_df"].query("term != 'Intercept'").copy()
    use_standardized = effect_df["absolute_standardized_effect"].notna().any()
    value_col = (
        "absolute_standardized_effect"
        if use_standardized
        else "absolute_coefficient"
    )
    x_label = (
        "Absolute standardized effect (|t|)"
        if use_standardized
        else "Absolute model coefficient"
    )

    effect_df = effect_df.sort_values(value_col, ascending=True)
    colors = [factor_color_for_term(term) for term in effect_df["term"]]

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.barh(effect_df["term"], effect_df[value_col], color=colors)
    ax.set_xlabel(x_label)
    ax.set_title(
        title or f"Copper DSD - Pareto chart for {response_label(model['response'])}"
    )
    ax.grid(axis="x", alpha=0.25)
    if use_standardized and model["residual_df"] > 0:
        try:
            from scipy import stats

            critical_value = stats.t.ppf(1 - alpha / 2, df=model["residual_df"])
            ax.axvline(
                critical_value,
                color="#d55e00",
                linewidth=1.5,
                label=f"alpha={alpha:g} threshold",
            )
            if secondary_alpha is not None:
                secondary_critical_value = stats.t.ppf(
                    1 - secondary_alpha / 2,
                    df=model["residual_df"],
                )
                ax.axvline(
                    secondary_critical_value,
                    color="#333333",
                    linewidth=1.5,
                    linestyle="--",
                    label=f"alpha={secondary_alpha:g} threshold",
                )
            ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, -0.16),
                ncol=1,
                frameon=True,
                facecolor="white",
                framealpha=0.92,
                edgecolor="#d0d0d0",
            )
        except Exception:
            pass
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.28)
    plt.show(block=False)


def plot_copper_normal_effects(model, half_normal=False, title=None):
    """Normal or half-normal plot of fitted effects."""
    effect_df = model["coefficient_df"].query("term != 'Intercept'").copy()
    effect_df["effect"] = effect_df["coefficient"]
    if half_normal:
        effect_df["plot_effect"] = effect_df["effect"].abs()
        effect_df = effect_df.sort_values("plot_effect")
    else:
        effect_df["plot_effect"] = effect_df["effect"]
        effect_df = effect_df.sort_values("plot_effect")

    normal = NormalDist()
    n = len(effect_df)
    probabilities = (np.arange(1, n + 1) - 0.5) / n
    if half_normal:
        theoretical = [normal.inv_cdf((p + 1) / 2) for p in probabilities]
        x_label = "Half-normal theoretical quantile"
    else:
        theoretical = [normal.inv_cdf(p) for p in probabilities]
        x_label = "Normal theoretical quantile"

    colors = [factor_color_for_term(term) for term in effect_df["term"]]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.scatter(
        theoretical,
        effect_df["plot_effect"],
        color=colors,
        edgecolor="black",
        linewidth=0.4,
        zorder=3,
    )
    for x, (_, row) in zip(theoretical, effect_df.iterrows()):
        ax.annotate(short_label_for_term(row["term"]), (x, row["plot_effect"]), xytext=(4, 4),
                    textcoords="offset points", fontsize=8)
    ax.axhline(0, color="#333333", linewidth=1)
    ax.set_xlabel(x_label)
    unit = response_unit(model["response"])
    effect_label = "Absolute effect" if half_normal else "Signed effect"
    if unit:
        effect_label = f"{effect_label} ({unit})"
    ax.set_ylabel(effect_label)
    ax.set_title(title or (
        f"Copper DSD - {'Half-normal' if half_normal else 'Normal'} plot for "
        f"{response_label(model['response'])}"
    ))
    ax.grid(alpha=0.25)
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor="black",
            markeredgewidth=0.4,
            label=factor,
        )
        for factor, color in FACTOR_COLORS.items()
        if factor in model["factor_names"]
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=8, loc="best")
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_main_effects(
    coded_with_response_df,
    levels,
    response,
    factor_names,
    max_columns=None,
):
    """Plot raw response means at low, center, and high factor levels."""
    plotted_means = []
    for factor in factor_names:
        means = (
            coded_with_response_df.dropna(subset=[response])
            .groupby(factor)[response]
            .mean()
            .reindex([-1, 0, 1])
        )
        plotted_means.extend(means.dropna().tolist())

    y_min = min(plotted_means)
    y_max = max(plotted_means)
    y_padding = 0.12 * (y_max - y_min) if y_max > y_min else 0.1

    n_columns = min(max_columns or len(factor_names), len(factor_names))
    n_rows = int(np.ceil(len(factor_names) / n_columns))
    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(4.3 * n_columns + 1.0, 5.0 * n_rows),
        sharey=True,
        squeeze=False,
    )
    axes = axes.ravel()
    response_label_y = -0.23 if n_rows > 1 else -0.32

    for ax, factor in zip(axes, factor_names):
        means = (
            coded_with_response_df.dropna(subset=[response])
            .groupby(factor)[response]
            .mean()
            .reindex([-1, 0, 1])
        )
        factor_color = FACTOR_COLORS.get(factor, "#08519c")
        response_diff = means.max() - means.min()
        x_values = [levels[factor][coded] for coded in [-1, 0, 1]]
        ax.set_xticks(x_values)
        ax.set_xticklabels([f"{value:g}" for value in x_values])
        factor_label = FACTOR_LABELS.get(factor, factor)
        ax.set_title(factor_label.split(" (")[0], fontsize=15)
        ax.set_xlabel(factor_label, fontsize=14)
        ax.set_ylabel(response_label(response), fontsize=14)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
        y_bottom, y_top = ax.get_ylim()
        y_safe = 0.09 * (y_top - y_bottom)

        ax.plot(x_values, means, marker="o", color=factor_color)
        valid_points = [
            (x_value, mean_value)
            for x_value, mean_value in zip(x_values, means)
            if not pd.isna(mean_value)
        ]
        line_slope = valid_points[-1][1] - valid_points[0][1] if len(valid_points) > 1 else 0
        for point_index, (x_value, mean_value) in enumerate(valid_points):
            if pd.isna(mean_value):
                continue
            if len(valid_points) == 1:
                label_x_axes = 0.5
                ha = "center"
            elif point_index == 0:
                label_x_axes = 0.08
                ha = "left"
            elif point_index == len(valid_points) - 1:
                label_x_axes = 0.92
                ha = "right"
            else:
                label_x_axes = 0.5
                ha = "center"

            if line_slope > 0:
                preferred_above = point_index == 0
            elif line_slope < 0:
                preferred_above = point_index == len(valid_points) - 1
            else:
                preferred_above = True

            if mean_value > y_top - 1.8 * y_safe:
                label_y = mean_value - y_safe
                va = "top"
            elif mean_value < y_bottom + 1.8 * y_safe:
                label_y = mean_value + y_safe
                va = "bottom"
            elif preferred_above:
                label_y = mean_value + y_safe
                va = "bottom"
            else:
                label_y = mean_value - y_safe
                va = "top"

            ax.text(
                label_x_axes,
                label_y,
                f"{mean_value:.2f}",
                transform=ax.get_yaxis_transform(),
                ha=ha,
                va=va,
                fontsize=13,
                color=factor_color,
                clip_on=True,
                bbox={
                    "boxstyle": "round,pad=0.15",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.78,
                },
            )

        ax.tick_params(axis="both", labelsize=12)
        ax.tick_params(axis="y", labelleft=True)
        ax.set_box_aspect(1)
        ax.grid(alpha=0.25)
        ax.text(
            0.5,
            response_label_y,
            (
                rf"Response range: $\Delta$ = {response_diff:.2f} "
                f"{response_unit(response)}"
            ).rstrip(),
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=13,
            color=factor_color,
        )

    for ax in axes[len(factor_names):]:
        ax.set_visible(False)

    fig.suptitle(
        f"Copper DSD - Main effects on {response_label(response)}",
        fontsize=18,
    )
    fig.subplots_adjust(
        left=0.06,
        right=0.985,
        bottom=0.16 if n_rows > 1 else 0.29,
        top=0.90 if n_rows > 1 else 0.82,
        wspace=0.15,
        hspace=0.45 if n_rows > 1 else 0.2,
    )
    plt.show(block=False)


def plot_copper_curvature(model):
    """Plot individual quadratic coefficients in the reduced model."""
    curvature_df = model["coefficient_df"].copy()
    curvature_df = curvature_df[curvature_df["term"].str.endswith("^2")]
    if curvature_df.empty:
        print("No quadratic terms are included in the current model.")
        return

    colors = [factor_color_for_term(term) for term in curvature_df["term"]]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(curvature_df["term"], curvature_df["coefficient"], color=colors)
    ax.axhline(0, color="#333333", linewidth=1)
    unit = response_unit(model["response"])
    y_label = "Quadratic coefficient"
    if unit:
        y_label = f"{y_label} ({unit})"
    ax.set_ylabel(y_label)
    ax.set_title(
        "Copper DSD - Individual curvature terms for "
        f"{response_label(model['response'])}"
    )
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_interactions(
    coded_with_response_df,
    response,
    interaction_pairs,
    figsize=None,
    font_size=None,
):
    """Plot raw means for selected physically plausible interactions."""
    plot_df = coded_with_response_df.dropna(subset=[response]).copy()
    if plot_df.empty:
        print(f"No measured {response} values are available for interaction plots.")
        return

    valid_pairs = [
        (factor_x, factor_line)
        for factor_x, factor_line in interaction_pairs
        if factor_x in plot_df.columns and factor_line in plot_df.columns
    ]
    if not valid_pairs:
        print("No valid interaction pairs are available for the current data.")
        return

    n_pairs = len(valid_pairs)
    if n_pairs > 3:
        fig = plt.figure(figsize=figsize or (13.2, 9.5))
        grid = fig.add_gridspec(2, 6)
        spans = [(0, 0, 2), (0, 2, 4), (0, 4, 6)]
        remaining_count = n_pairs - 3
        if remaining_count == 1:
            spans.append((1, 2, 4))
        else:
            spans.extend([(1, 1, 3), (1, 3, 5)])

        axes = []
        for row, column_start, column_stop in spans[:n_pairs]:
            axes.append(
                fig.add_subplot(
                    grid[row, column_start:column_stop],
                    sharey=axes[0] if axes else None,
                )
            )
    else:
        fig, axes = plt.subplots(
            1,
            n_pairs,
            figsize=figsize or (4.2 * n_pairs, 5.4),
            sharey=True,
        )
        if n_pairs == 1:
            axes = [axes]

    for ax, (factor_x, factor_line) in zip(axes, valid_pairs):
        factor_x_label = INTERACTION_FACTOR_LABELS.get(factor_x, factor_x)
        factor_line_label = INTERACTION_FACTOR_LABELS.get(factor_line, factor_line)
        for line_level in [-1, 0, 1]:
            subset = plot_df[plot_df[factor_line] == line_level]
            means = subset.groupby(factor_x)[response].mean().reindex([-1, 0, 1])
            if means.notna().any():
                ax.plot(
                    [-1, 0, 1],
                    means,
                    marker="o",
                    label=f"{factor_line_label}: {line_level}",
                )

        ax.set_xlabel(
            f"{factor_x_label} coded level",
            fontsize=font_size if font_size is not None else None,
        )
        ax.set_title(
            f"{factor_x_label} x {factor_line_label}",
            fontsize=font_size + 1 if font_size is not None else None,
        )
        ax.set_xticks([-1, 0, 1])
        if font_size is not None:
            ax.tick_params(axis="both", labelsize=font_size - 1)
        ax.grid(alpha=0.25)
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.30),
            frameon=False,
            fontsize=font_size - 1 if font_size is not None else 10,
            ncol=1,
        )

    axes[0].set_ylabel(
        response_label(response),
        fontsize=font_size if font_size is not None else None,
    )
    if n_pairs > 3:
        axes[3].set_ylabel(
            response_label(response),
            fontsize=font_size if font_size is not None else None,
        )
    fig.suptitle(
        f"Copper DSD - Selected interaction plots for {response_label(response)}",
        fontsize=font_size + 3 if font_size is not None else None,
    )
    if n_pairs > 3:
        fig.subplots_adjust(
            left=0.07,
            right=0.985,
            bottom=0.16,
            top=0.88,
            wspace=0.55,
            hspace=0.72,
        )
    else:
        fig.subplots_adjust(
            left=0.055,
            right=0.99,
            bottom=0.32,
            top=0.86,
            wspace=0.25,
        )
    plt.show(block=False)


def build_copper_paired_measurement_table(
    measurement_table,
    measurement_repeats_col="weight_measurements_g",
    exclude_runs=None,
):
    """Return one row per specimen having at least two weight measurements."""
    excluded = {str(run) for run in (exclude_runs or [])}
    paired_rows = []
    for _, row in measurement_table.iterrows():
        if str(row["Run"]) in excluded:
            continue
        values = row.get(measurement_repeats_col)
        if not isinstance(values, (list, tuple, np.ndarray, pd.Series)):
            continue
        numeric_values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
        if len(numeric_values) < 2:
            continue
        measurement_1 = float(numeric_values.iloc[0])
        measurement_2 = float(numeric_values.iloc[1])
        paired_rows.append(
            {
                "Run": row["Run"],
                "measurement_1_g": measurement_1,
                "measurement_2_g": measurement_2,
                "offset_2_minus_1_g": measurement_2 - measurement_1,
            }
        )
    return pd.DataFrame(paired_rows)


def plot_copper_measurement_correlation(
    measurement_table,
    measurement_repeats_col="weight_measurements_g",
    exclude_runs=None,
    value_digits=3,
):
    """Plot agreement and offsets between the first two copper weighings."""
    paired = build_copper_paired_measurement_table(
        measurement_table,
        measurement_repeats_col=measurement_repeats_col,
        exclude_runs=exclude_runs,
    )
    if len(paired) < 2:
        raise ValueError("At least two specimens with paired measurements are required.")

    measurement_1 = paired["measurement_1_g"].to_numpy(dtype=float)
    measurement_2 = paired["measurement_2_g"].to_numpy(dtype=float)
    offsets = paired["offset_2_minus_1_g"].to_numpy(dtype=float)
    mean_offset = offsets.mean()
    offset_std = offsets.std(ddof=1)
    correlation = np.corrcoef(measurement_1, measurement_2)[0, 1]
    slope, intercept = np.polyfit(measurement_1, measurement_2, 1)

    fig, (ax_corr, ax_offset) = plt.subplots(
        1,
        2,
        figsize=(13, 4.8),
        gridspec_kw={"width_ratios": [1, 1.35]},
    )

    low = min(measurement_1.min(), measurement_2.min()) - 0.4
    high = max(measurement_1.max(), measurement_2.max()) + 0.4
    fit_x = np.array([low, high])
    ax_corr.scatter(measurement_1, measurement_2, color="#2171b5", s=45, zorder=3)
    for run, x_value, y_value in zip(paired["Run"], measurement_1, measurement_2):
        ax_corr.annotate(
            run,
            (x_value, y_value),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    ax_corr.plot(fit_x, fit_x, color="#777777", linestyle="--", label="no offset: y = x")
    ax_corr.plot(
        fit_x,
        fit_x + mean_offset,
        color="#31a354",
        linestyle=":",
        linewidth=2,
        label=f"mean offset: +{mean_offset:.{value_digits}f} g",
    )
    ax_corr.plot(
        fit_x,
        slope * fit_x + intercept,
        color="#cb181d",
        linewidth=1.8,
        label="linear fit",
    )
    ax_corr.set_xlim(low, high)
    ax_corr.set_ylim(low, high)
    ax_corr.set_aspect("equal", adjustable="box")
    ax_corr.set_xlabel("Measurement 1 (g)")
    ax_corr.set_ylabel("Measurement 2 (g)")
    ax_corr.set_title("Measurement 1 vs measurement 2")
    ax_corr.text(
        0.04,
        0.96,
        f"Pearson r = {correlation:.3f}\n$R^2$ = {correlation ** 2:.3f}",
        transform=ax_corr.transAxes,
        ha="left",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.9},
    )
    ax_corr.grid(alpha=0.25)
    ax_corr.legend(frameon=False, fontsize=8, loc="lower right")

    run_x = np.arange(len(paired), dtype=float)
    ax_offset.plot(run_x, offsets, color="#756bb1", marker="o", linewidth=1.6)
    ax_offset.axhline(
        mean_offset,
        color="#31a354",
        linestyle="--",
        linewidth=1.8,
        label=f"mean = {mean_offset:.{value_digits}f} g",
    )
    ax_offset.fill_between(
        [-0.5, len(paired) - 0.5],
        mean_offset - offset_std,
        mean_offset + offset_std,
        color="#31a354",
        alpha=0.12,
        label=f"mean ± 1 SD ({offset_std:.{value_digits}f} g)",
    )
    ax_offset.set_xticks(run_x)
    ax_offset.set_xticklabels(paired["Run"], rotation=45, ha="right")
    ax_offset.set_ylabel("Measurement 2 − measurement 1 (g)")
    ax_offset.set_title("Offset consistency by specimen")
    ax_offset.grid(axis="y", alpha=0.25)
    ax_offset.legend(frameon=False, fontsize=8)

    fig.suptitle("Copper repeated mass measurements")
    plt.tight_layout()
    plt.show(block=False)
    return paired


def plot_copper_weight_measurements(
    measurement_table,
    response="weight_g",
    measurement_repeats_col="weight_measurements_g",
    measurement_1_offset_g=None,
    offset_excluded_runs=None,
    base_runs=None,
    repeat_runs=None,
    extra_runs=None,
    y_min=None,
    y_max=None,
    value_digits=2,
    reference_weight_g=None,
    reference_label="m_geom",
    show_weight_line=True,
    measurement_indices=None,
    title="Copper DSD - Mass measurements",
):
    """Plot individual weighings, print repeats, and the mean for each copper run."""
    df = measurement_table.copy()
    df[response] = pd.to_numeric(df[response], errors="coerce")
    df = df.dropna(subset=[response])

    if base_runs is None:
        base_runs = [f"T{i}" for i in range(1, 14)]
    if repeat_runs is None:
        repeat_runs = [run for run in df["Run"] if "." in str(run)]
    if extra_runs is None:
        extra_runs = [
            run
            for run in df["Run"]
            if str(run).startswith("T")
            and "." not in str(run)
            and str(run) not in base_runs
        ]

    base_df = (
        df[df["Run"].isin(base_runs)]
        .set_index("Run")
        .reindex(base_runs)
        .reset_index()
        .dropna(subset=[response])
    )
    repeat_df = df[df["Run"].isin(repeat_runs)].copy()
    extra_df = (
        df[df["Run"].isin(extra_runs)]
        .set_index("Run")
        .reindex(extra_runs)
        .reset_index()
        .dropna(subset=[response])
    )

    x = np.arange(len(base_df), dtype=float)
    bar_width = 0.22
    tick_positions = x.tolist()
    tick_labels = base_df["Run"].tolist()
    extra_x = np.arange(
        len(base_df) + 1,
        len(base_df) + 1 + len(extra_df),
        dtype=float,
    )
    original_run_positions = {
        **{run: position for position, run in enumerate(base_df["Run"])},
        **dict(zip(extra_df["Run"], extra_x)),
    }
    original_measurements = pd.concat(
        [base_df[["Run", response]], extra_df[["Run", response]]],
        ignore_index=True,
    )

    original_plot_df = pd.concat(
        [
            base_df.assign(_plot_x=x),
            extra_df.assign(_plot_x=extra_x),
        ],
        ignore_index=True,
    )

    def _measurement_values(row):
        values = row.get(measurement_repeats_col)
        if isinstance(values, (list, tuple, np.ndarray, pd.Series)):
            numeric_values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
        elif pd.isna(values):
            numeric_values = pd.Series(dtype=float)
        else:
            numeric_values = pd.to_numeric(pd.Series([values]), errors="coerce").dropna()
        if numeric_values.empty:
            return [float(row[response])]
        return numeric_values.astype(float).tolist()

    original_plot_df["_measurements"] = original_plot_df.apply(
        _measurement_values,
        axis=1,
    )
    offset_excluded = {str(run) for run in (offset_excluded_runs or [])}
    applied_offset = (
        float(measurement_1_offset_g)
        if measurement_1_offset_g is not None and pd.notna(measurement_1_offset_g)
        else 0.0
    )

    def _corrected_average(row):
        values = list(row["_measurements"])
        if (
            applied_offset
            and len(values) >= 2
            and str(row["Run"]) not in offset_excluded
        ):
            values[0] += applied_offset
        return float(np.mean(values))

    original_plot_df["_plot_average"] = original_plot_df.apply(
        _corrected_average,
        axis=1,
    )
    max_measurements = max(original_plot_df["_measurements"].map(len), default=1)
    if measurement_indices is None:
        plotted_measurement_indices = list(range(max_measurements))
    else:
        plotted_measurement_indices = sorted(
            {
                int(index)
                for index in measurement_indices
                if 0 <= int(index) < max_measurements
            }
        )
        if not plotted_measurement_indices:
            raise ValueError("measurement_indices does not select an available measurement")

    fig, ax = plt.subplots(figsize=(16, 4.3))
    if reference_weight_g is not None and pd.notna(reference_weight_g):
        ax.axhline(
            reference_weight_g,
            color="#888888",
            linestyle="-",
            linewidth=1.4,
            alpha=0.55,
            zorder=0,
            label=reference_label,
        )
        ax.annotate(
            f"{reference_label} = {reference_weight_g:.{value_digits}f} g",
            xy=(0.01, reference_weight_g),
            xycoords=("axes fraction", "data"),
            xytext=(0, 5),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=9,
            color="#555555",
        )

    measurement_colors = ["#9ecae1", "#6baed6", "#3182bd", "#08519c"]
    measurement_offsets = (
        np.arange(len(plotted_measurement_indices), dtype=float)
        - (len(plotted_measurement_indices) - 1) / 2
    ) * bar_width
    for plotted_index, measurement_index in enumerate(plotted_measurement_indices):
        measurement_rows = original_plot_df[
            original_plot_df["_measurements"].map(len) > measurement_index
        ]
        measurement_x = (
            measurement_rows["_plot_x"].to_numpy(dtype=float)
            + measurement_offsets[plotted_index]
        )
        measurement_y = measurement_rows["_measurements"].map(
            lambda values: values[measurement_index]
        )
        color = measurement_colors[min(measurement_index, len(measurement_colors) - 1)]
        ax.bar(
            measurement_x,
            measurement_y,
            width=bar_width * 0.88,
            color=color,
            label=f"measurement {measurement_index + 1}",
            alpha=0.9,
        )
        for xi, mass in zip(measurement_x, measurement_y):
            ax.annotate(
                f"{mass:.{value_digits}f}",
                xy=(xi, mass),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                color="#333333",
            )

        if measurement_index == 0 and applied_offset:
            correction_rows = measurement_rows[
                (measurement_rows["_measurements"].map(len) >= 2)
                & (~measurement_rows["Run"].astype(str).isin(offset_excluded))
            ]
            correction_x = (
                correction_rows["_plot_x"].to_numpy(dtype=float)
                + measurement_offsets[plotted_index]
            )
            correction_bottom = correction_rows["_measurements"].map(
                lambda values: values[measurement_index]
            ).to_numpy(dtype=float)
            ax.bar(
                correction_x,
                np.full(len(correction_rows), applied_offset),
                bottom=correction_bottom,
                width=bar_width * 0.88,
                color="#fdae6b",
                edgecolor="#e6550d",
                hatch="///",
                alpha=0.75,
                label=f"M1 offset added (+{applied_offset:.{value_digits}f} g)",
            )
            for xi, bottom in zip(correction_x, correction_bottom):
                ax.annotate(
                    f"+{applied_offset:.{value_digits}f}",
                    xy=(xi, bottom + applied_offset / 2),
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="#7f2704",
                    rotation=90,
                )

    if show_weight_line:
        ax.plot(
            original_plot_df["_plot_x"],
            original_plot_df["_plot_average"],
            color="#000000",
            marker="o",
            linewidth=2,
            markersize=5,
            label=(
                "offset-corrected average"
                if applied_offset
                else "average measurement"
            ),
            zorder=4,
        )

    repeat_summary_lines = []
    if not repeat_df.empty:
        repeat_df["original_run"] = repeat_df["Run"].astype(str).str.split(".").str[0]
        repeat_df["x"] = repeat_df["original_run"].map(original_run_positions)
        repeat_df = repeat_df.dropna(subset=["x", response])
        if not repeat_df.empty:
            repeat_df["_repeat_index"] = repeat_df.groupby("original_run").cumcount()
            repeat_x = (
                repeat_df["x"].to_numpy(dtype=float)
                + bar_width * (1.6 + repeat_df["_repeat_index"].to_numpy(dtype=float))
            )
            ax.bar(
                repeat_x,
                repeat_df[response],
                width=bar_width,
                color="#31a354",
                label="repeat",
                alpha=0.85,
            )
            for xi, mass in zip(repeat_x, repeat_df[response]):
                ax.annotate(
                    f"{mass:.{value_digits}f}",
                    xy=(xi, mass),
                    xytext=(0, 8),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

            repeat_summary_df = (
                repeat_df.groupby("original_run", as_index=False)
                .agg({response: "mean"})
            )
            comparison_df = repeat_summary_df.merge(
                original_measurements,
                left_on="original_run",
                right_on="Run",
                suffixes=("_repeat", "_original"),
            )
            if not comparison_df.empty:
                abs_diff = (
                    comparison_df[f"{response}_repeat"]
                    - comparison_df[f"{response}_original"]
                ).abs()
                repeat_summary_lines.append(
                    "Repeatability: "
                    f"mean abs diff {abs_diff.mean():.3f} g; "
                    f"max {abs_diff.max():.3f} g"
                )

    if not extra_df.empty:
        tick_positions.extend(extra_x.tolist())
        tick_labels.extend(extra_df["Run"].tolist())

    if repeat_summary_lines:
        ax.text(
            0.99,
            0.04,
            "\n".join(repeat_summary_lines),
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

    ax.set_xlabel("Test plate")
    ax.set_ylabel("Mass (g)")
    ax.set_title(title)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45)
    if y_min is not None or y_max is not None:
        ax.set_ylim(bottom=y_min, top=y_max)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_actual_vs_predicted(model):
    """Actual versus model-predicted response plot."""
    fitted_df = model["fitted_df"]
    response = model["response"]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(fitted_df[response], fitted_df["predicted_response"], color="#08519c")
    for _, row in fitted_df.iterrows():
        ax.annotate(row["Run"], (row[response], row["predicted_response"]),
                    xytext=(4, 4), textcoords="offset points", fontsize=8)
    low = min(fitted_df[response].min(), fitted_df["predicted_response"].min())
    high = max(fitted_df[response].max(), fitted_df["predicted_response"].max())
    ax.plot([low, high], [low, high], color="#333333", linestyle="--")
    ax.set_xlabel(f"Actual {response_label(response)}")
    ax.set_ylabel(f"Predicted {response_label(response)}")
    ax.set_title(
        f"Copper DSD - Actual vs predicted {response_label(response)}"
    )
    ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_residuals(model, factor_names):
    """Residual diagnostics: predicted, normal plot, print order, and factors."""
    fitted_df = model["fitted_df"]
    residual = fitted_df["residual"]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].scatter(fitted_df["predicted_response"], residual, color="#08519c")
    axes[0].axhline(0, color="#333333", linewidth=1)
    axes[0].set_xlabel(f"Predicted {response_label(model['response'])}")
    axes[0].set_ylabel(residual_label(model["response"]))
    axes[0].set_title("Residuals vs predicted")

    sorted_residual = np.sort(residual.to_numpy(dtype=float))
    n = len(sorted_residual)
    normal = NormalDist()
    theoretical = [normal.inv_cdf((i - 0.5) / n) for i in range(1, n + 1)]
    axes[1].scatter(theoretical, sorted_residual, color="#08519c")
    axes[1].axhline(0, color="#333333", linewidth=1)
    axes[1].set_xlabel("Normal theoretical quantile")
    axes[1].set_title("Normal plot of residuals")

    axes[2].plot(fitted_df["print_order"], residual, marker="o", color="#08519c")
    axes[2].axhline(0, color="#333333", linewidth=1)
    axes[2].set_xlabel("Print order")
    axes[2].set_title("Residuals vs print order")

    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle(
        "Copper DSD - Residual diagnostics for "
        f"{response_label(model['response'])}"
    )
    plt.tight_layout()
    plt.show(block=False)

    fig, axes = plt.subplots(1, len(factor_names), figsize=(3.2 * len(factor_names), 3.4))
    if len(factor_names) == 1:
        axes = [axes]
    for ax, factor in zip(axes, factor_names):
        ax.scatter(fitted_df[factor], residual, color=FACTOR_COLORS.get(factor, "#08519c"))
        ax.axhline(0, color="#333333", linewidth=1)
        ax.set_xlabel(f"{factor} coded")
        ax.set_ylabel(residual_label(model["response"]))
        ax.grid(alpha=0.25)
    fig.suptitle(
        "Copper DSD - Residuals vs factors for "
        f"{response_label(model['response'])}"
    )
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_contours(
    model,
    levels,
    factor_pairs,
    fixed_level=0,
    figsize=None,
    font_size=14,
):
    """Contour plots for selected two-factor slices of a reduced final model."""
    terms = model["terms"]
    response = model["response"]
    response_text = response_label(response)
    factor_names = list(levels)
    grid = np.linspace(-1, 1, 50)
    factor_pairs = list(factor_pairs)
    if not factor_pairs:
        return

    max_cols = 3
    n_cols = min(max_cols, len(factor_pairs))
    n_rows = int(np.ceil(len(factor_pairs) / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize or (6.4 * n_cols, 6.4 * n_rows),
        squeeze=False,
    )
    axes_flat = axes.ravel()

    for ax, (factor_x, factor_y) in zip(axes_flat, factor_pairs):
        xx, yy = np.meshgrid(grid, grid)
        coded_grid = pd.DataFrame({factor: fixed_level for factor in factor_names},
                                  index=np.arange(xx.size))
        coded_grid[factor_x] = xx.ravel()
        coded_grid[factor_y] = yy.ravel()
        model_matrix = build_copper_model_matrix(coded_grid, terms)
        zz = (model_matrix @ model["coefficients"]).to_numpy().reshape(xx.shape)

        x_real = np.interp(grid, [-1, 0, 1], [levels[factor_x][i] for i in [-1, 0, 1]])
        y_real = np.interp(grid, [-1, 0, 1], [levels[factor_y][i] for i in [-1, 0, 1]])

        contour = ax.contourf(x_real, y_real, zz, levels=18, cmap="viridis")
        ax.contour(x_real, y_real, zz, levels=8, colors="white", linewidths=0.5)
        colorbar = fig.colorbar(
            contour,
            ax=ax,
            fraction=0.046,
            pad=0.04,
        )
        colorbar.set_label(
            f"Predicted {response_text}",
            fontsize=font_size - 1,
            labelpad=7,
        )
        colorbar.ax.tick_params(labelsize=font_size - 3)

        factor_x_label = INTERACTION_FACTOR_LABELS.get(factor_x, factor_x)
        factor_y_label = INTERACTION_FACTOR_LABELS.get(factor_y, factor_y)

        measured_points = model["fitted_df"][[factor_x, factor_y]].dropna()
        measured_x = np.interp(
            measured_points[factor_x],
            [-1, 0, 1],
            [levels[factor_x][i] for i in [-1, 0, 1]],
        )
        measured_y = np.interp(
            measured_points[factor_y],
            [-1, 0, 1],
            [levels[factor_y][i] for i in [-1, 0, 1]],
        )
        projected_points = np.column_stack([measured_x, measured_y])
        ax.scatter(
            projected_points[:, 0],
            projected_points[:, 1],
            s=58,
            facecolor="white",
            edgecolor="#202020",
            linewidth=1.2,
            alpha=0.95,
            zorder=6,
            clip_on=False,
        )

        ax.set_xlabel(factor_x_label, fontsize=font_size)
        ax.set_ylabel(factor_y_label, fontsize=font_size)
        ax.set_title(
            f"{factor_x_label} vs {factor_y_label}",
            fontsize=font_size + 1,
            pad=16,
        )
        ax.tick_params(axis="both", labelsize=font_size - 2)
        ax.set_box_aspect(1)

    for ax in axes_flat[len(factor_pairs):]:
        ax.set_visible(False)

    fig.suptitle(
        f"Copper DSD - Response contours for {response_text}",
        fontsize=font_size + 4,
        y=0.995,
    )
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="white",
                markeredgecolor="#202020",
                markeredgewidth=1.2,
                markersize=7,
                label="Projected measured design point",
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        ncol=1,
        frameon=True,
        framealpha=0.9,
        fontsize=font_size - 1,
    )
    fig.subplots_adjust(
        left=0.055,
        right=0.985,
        bottom=0.105,
        top=0.91,
        wspace=0.31,
        hspace=0.24,
    )
    plt.show(block=False)


def build_copper_desirability_table(
    measurement_table,
    density_col="density_g_cm3",
    porosity_col="porosity_percent",
    leakage_col="log_leakage_rate",
    clumping_col="visible_clumping_score",
):
    """Build a simple multiresponse desirability table from available responses."""
    df = measurement_table.copy()

    def maximize(series):
        valid = series.dropna()
        if valid.empty or valid.max() == valid.min():
            return pd.Series(np.nan, index=series.index)
        return (series - valid.min()) / (valid.max() - valid.min())

    def minimize(series):
        valid = series.dropna()
        if valid.empty or valid.max() == valid.min():
            return pd.Series(np.nan, index=series.index)
        return (valid.max() - series) / (valid.max() - valid.min())

    desirability_cols = []
    if density_col in df:
        df["d_density"] = maximize(pd.to_numeric(df[density_col], errors="coerce"))
        desirability_cols.append("d_density")
    if porosity_col in df:
        df["d_porosity"] = minimize(pd.to_numeric(df[porosity_col], errors="coerce"))
        desirability_cols.append("d_porosity")
    if leakage_col in df:
        df["d_leakage"] = minimize(pd.to_numeric(df[leakage_col], errors="coerce"))
        desirability_cols.append("d_leakage")
    if clumping_col in df:
        df["d_clumping"] = minimize(pd.to_numeric(df[clumping_col], errors="coerce"))
        desirability_cols.append("d_clumping")

    if desirability_cols:
        df["overall_desirability"] = df[desirability_cols].apply(
            lambda row: np.exp(np.log(row.dropna()).mean()) if row.notna().any() else np.nan,
            axis=1,
        )
    else:
        df["overall_desirability"] = np.nan

    return df


def plot_copper_desirability(desirability_table):
    """Plot overall desirability by run."""
    plot_df = desirability_table.dropna(subset=["overall_desirability"]).copy()
    if plot_df.empty:
        print("No desirability plot yet: add density, porosity, leakage, or clumping responses.")
        return

    plot_df = plot_df.sort_values("overall_desirability")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(plot_df["Run"], plot_df["overall_desirability"], color="#2b8cbe")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Overall desirability (0 = poor, 1 = best)")
    ax.set_title("Copper DSD - Multiresponse desirability")
    ax.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)
