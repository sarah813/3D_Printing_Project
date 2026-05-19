import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Numeric literature values
# -----------------------------

numeric_data = [
    # parameter, value, unit, label, effect
    ("Layer height", 0.09, "mm", "0.09 mm", "better vacuum tightness"),
    ("Layer height", 0.10, "mm", "0.10 mm", "better hole accuracy"),
    ("Layer height", 0.15, "mm", "0.15 mm", "good surface quality"),

    ("Nozzle temperature", 215, "°C", "215 °C", "high density"),
    ("Nozzle temperature", 230, "°C", "230 °C", "minimum porosity"),
    ("Nozzle temperature", 240, "°C", "240 °C", "best PETG channel accuracy"),

    ("Bed temperature", 50, "°C", "50 °C", "dimensional stability"),
    ("Bed temperature", 60, "°C", "60 °C", "dimensional stability"),
    ("Bed temperature", 70, "°C", "70 °C", "best PETG channel accuracy"),

    ("Print speed", 20, "mm/s", "20 mm/s", "good surface quality"),
    ("Print speed", 25, "mm/s", "25 mm/s", "minimum porosity"),
    ("Print speed", 30, "mm/s", "30 mm/s", "better vacuum tightness"),

    ("Flow rate", 100, "%", "100 %", "best channel accuracy"),
    ("Flow rate", 110, "%", "110 %", "better vacuum tightness"),

    ("Line width / raster width", 0.5, "mm", "0.5 mm", "minimum porosity"),

    ("Wall thickness", 6, "mm", "6 mm", "better vacuum tightness"),

]

df_num = pd.DataFrame(
    numeric_data,
    columns=["Parameter", "Value", "Unit", "Label", "Effect"]
)

parameter_order = [
    "Layer height",
    "Nozzle temperature",
    "Bed temperature",
    "Print speed",
    "Flow rate",
    "Line width / raster width",
    "Wall thickness",
    #"Channel diameter"
]

colors = {
    "Layer height": "#4E79A7",
    "Nozzle temperature": "#F28E2B",
    "Bed temperature": "#E15759",
    "Print speed": "#76B7B2",
    "Flow rate": "#59A14F",
    "Line width / raster width": "#EDC948",
    "Wall thickness": "#B07AA1",
    #"Channel diameter": "#FF9DA7",
}

# -----------------------------
# Plot numeric parameters
# -----------------------------

fig, axes = plt.subplots(
    nrows=len(parameter_order),
    ncols=1,
    figsize=(9, 10),
    sharex=False
)

for ax, parameter in zip(axes, parameter_order):
    group = df_num[df_num["Parameter"] == parameter].sort_values("Value")

    y = np.arange(len(group))

    ax.barh(
        y,
        group["Value"],
        height=0.35,
        color=colors[parameter],
        edgecolor="black",
        linewidth=0.4
    )

    ax.set_yticks(y)
    ax.set_yticklabels(group["Label"], fontsize=9)

    ax.set_title(parameter, loc="left", fontsize=10, fontweight="bold", pad=4)
    ax.set_xlabel(group["Unit"].iloc[0], fontsize=8)

    # Add effect text at bar end
    xmax = group["Value"].max()
    for i, row in enumerate(group.itertuples()):
        ax.text(
            row.Value + xmax * 0.03,
            i,
            row.Effect,
            va="center",
            fontsize=8
        )

    ax.set_xlim(0, xmax * 1.55)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.suptitle(
    "Literature Derived FDM Parameter Values Linked to Improved Properties",
    fontsize=12,
    fontweight="bold",
    y=0.995
)

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()