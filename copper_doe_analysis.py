"""Analysis helpers for the copper definitive screening design."""

from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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


def plot_copper_pareto(model, title=None):
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
    ax.set_title(title or f"Copper DSD - Pareto chart ({model['response']})")
    ax.grid(axis="x", alpha=0.25)
    plt.tight_layout()
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
    ax.set_ylabel("Absolute effect" if half_normal else "Signed effect")
    ax.set_title(title or (
        f"Copper DSD - {'Half-normal' if half_normal else 'Normal'} plot ({model['response']})"
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


def plot_copper_main_effects(coded_with_response_df, levels, response, factor_names):
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

    fig, axes = plt.subplots(
        1,
        len(factor_names),
        figsize=(3.2 * len(factor_names) + 1.0, 4.2),
        sharey=True,
    )
    if len(factor_names) == 1:
        axes = [axes]

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
        ax.set_title(factor)
        ax.set_xlabel("Real factor value")
        ax.set_ylabel(response)
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
                fontsize=8,
                color=factor_color,
                clip_on=True,
                bbox={
                    "boxstyle": "round,pad=0.15",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.78,
                },
            )

        ax.tick_params(axis="y", labelleft=True)
        ax.grid(alpha=0.25)
        ax.text(
            0.5,
            -0.26,
            f"Response range: Δ = {response_diff:.2f} g",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=9,
            color=factor_color,
        )

    fig.suptitle(f"Copper DSD - Main effects ({response})")
    fig.subplots_adjust(
        left=0.06,
        right=0.985,
        bottom=0.24,
        top=0.82,
        wspace=0.35,
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
    ax.set_ylabel(f"Quadratic coefficient on {model['response']}")
    ax.set_title(f"Copper DSD - Individual curvature terms ({model['response']})")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_interactions(coded_with_response_df, response, interaction_pairs):
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
    fig, axes = plt.subplots(
        1,
        n_pairs,
        figsize=(4.2 * n_pairs, 5.4),
        sharey=True,
    )
    if n_pairs == 1:
        axes = [axes]

    for ax, (factor_x, factor_line) in zip(axes, valid_pairs):
        for line_level in [-1, 0, 1]:
            subset = plot_df[plot_df[factor_line] == line_level]
            means = subset.groupby(factor_x)[response].mean().reindex([-1, 0, 1])
            if means.notna().any():
                ax.plot([-1, 0, 1], means, marker="o", label=f"{factor_line}: {line_level}")

        ax.set_xlabel(f"{factor_x} coded level")
        ax.set_title(f"{factor_x} x {factor_line}")
        ax.set_xticks([-1, 0, 1])
        ax.grid(alpha=0.25)
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.30),
            frameon=False,
            fontsize=10,
            ncol=1,
        )

    axes[0].set_ylabel(response)
    fig.suptitle(f"Copper DSD - Selected interaction plots ({response})")
    fig.subplots_adjust(
        left=0.055,
        right=0.99,
        bottom=0.32,
        top=0.86,
        wspace=0.25,
    )
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
    ax.set_xlabel(f"Actual {response}")
    ax.set_ylabel(f"Predicted {response}")
    ax.set_title(f"Copper DSD - Actual vs predicted ({response})")
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
    axes[0].set_xlabel("Predicted response")
    axes[0].set_ylabel("Residual")
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
    fig.suptitle(f"Copper DSD - Residual diagnostics ({model['response']})")
    plt.tight_layout()
    plt.show(block=False)

    fig, axes = plt.subplots(1, len(factor_names), figsize=(3.2 * len(factor_names), 3.4))
    if len(factor_names) == 1:
        axes = [axes]
    for ax, factor in zip(axes, factor_names):
        ax.scatter(fitted_df[factor], residual, color=FACTOR_COLORS.get(factor, "#08519c"))
        ax.axhline(0, color="#333333", linewidth=1)
        ax.set_xlabel(f"{factor} coded")
        ax.set_ylabel("Residual")
        ax.grid(alpha=0.25)
    fig.suptitle(f"Copper DSD - Residuals vs factors ({model['response']})")
    plt.tight_layout()
    plt.show(block=False)


def plot_copper_contours(model, levels, factor_pairs, fixed_level=0):
    """Contour plots for selected two-factor slices of a reduced final model."""
    terms = model["terms"]
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
        figsize=(5.4 * n_cols, 4.9 * n_rows),
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
        fig.colorbar(
            contour,
            ax=ax,
            label=f"Predicted {model['response']}",
            fraction=0.046,
            pad=0.04,
        )
        ax.set_xlabel(factor_x)
        ax.set_ylabel(factor_y)
        ax.set_title(
            f"{factor_x} vs {factor_y}",
            pad=14,
        )

    for ax in axes_flat[len(factor_pairs):]:
        ax.set_visible(False)

    fig.subplots_adjust(
        left=0.055,
        right=0.985,
        bottom=0.08,
        top=0.91,
        wspace=0.34,
        hspace=0.38,
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
