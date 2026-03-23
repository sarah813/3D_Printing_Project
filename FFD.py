import pandas as pd
from pyDOE import fracfact
import tkinter as tk
from tkinter import ttk


# 16-run fractional factorial design for 6 factors
design = fracfact('a b c d ab ac')

coded_df = pd.DataFrame(design, columns=[
    "A_flow_ratio",
    "B_nozzle_temp",
    "C_speed",
    "D_line_width",
    "E_layer_height",
    "F_walls"
])

coded_df.index += 1
coded_df.index.name = "Run"

levels = {
    "A_flow_ratio":   {-1: 1.00,  1: 1.08},
    "B_nozzle_temp":  {-1: 220,   1: 230},
    "C_speed":    {-1: 50,    1: 200},
    "D_line_width":   {-1: 0.62,  1: 0.75},
    "E_layer_height": {-1: 0.12,  1: 0.18},
    "F_walls":        {-1: 4,     1: 10},
}

real_df = coded_df.copy()
for col, mapping in levels.items():
    real_df[col] = real_df[col].map(mapping)

display_names = {
    "A_flow_ratio": "Flow Ratio (-)",
    "B_nozzle_temp": "Nozzle Temp (deg C)",
    "C_speed": "Print Speed (mm/s)",
    "D_line_width": "Line Width (mm)",
    "E_layer_height": "Layer Height (mm)",
    "F_walls": "Walls (-)",
}

# still print in terminal if you want
print("CODED MATRIX (-1 / +1):")
print(coded_df)

print("\nREAL VALUE MATRIX:")
print(real_df)


def show_dataframe_popup(df, title):
    window = tk.Toplevel()
    window.title(title)
    window.geometry("1080x520")
    window.configure(bg="#05070d")

    style = ttk.Style(window)
    style.theme_use("clam")

    style.configure(
        "Dark.Treeview",
        background="#08111f",
        fieldbackground="#08111f",
        foreground="#7fd2ff",
        rowheight=34,
        borderwidth=0,
        font=("Segoe UI", 10),
    )
    style.configure(
        "Dark.Treeview.Heading",
        background="#0d1b2a",
        foreground="#b8ecff",
        relief="flat",
        font=("Segoe UI Semibold", 10),
        padding=8,
    )
    style.map(
        "Dark.Treeview",
        background=[("selected", "#123a63")],
        foreground=[("selected", "#e8f7ff")],
    )
    style.configure(
        "Dark.Vertical.TScrollbar",
        troughcolor="#05070d",
        background="#123a63",
        arrowcolor="#b8ecff",
    )
    style.configure(
        "Dark.Horizontal.TScrollbar",
        troughcolor="#05070d",
        background="#123a63",
        arrowcolor="#b8ecff",
    )

    outer = tk.Frame(window, bg="#05070d", padx=18, pady=18)
    outer.pack(fill="both", expand=True)

    title_label = tk.Label(
        outer,
        text=title,
        bg="#05070d",
        fg="#d7f5ff",
        font=("Segoe UI Semibold", 16),
        anchor="w",
    )
    title_label.pack(fill="x", pady=(0, 12))

    subtitle_label = tk.Label(
        outer,
        text="Fractional factorial design table",
        bg="#05070d",
        fg="#63b9ea",
        font=("Segoe UI", 10),
        anchor="w",
    )
    subtitle_label.pack(fill="x", pady=(0, 14))

    table_frame = tk.Frame(outer, bg="#0a1220", highlightbackground="#1a4f7a", highlightthickness=1)
    table_frame.pack(fill="both", expand=True)

    tree = ttk.Treeview(table_frame, style="Dark.Treeview")
    tree.grid(row=0, column=0, sticky="nsew")

    scrollbar_y = ttk.Scrollbar(
        table_frame, orient="vertical", command=tree.yview, style="Dark.Vertical.TScrollbar"
    )
    scrollbar_y.grid(row=0, column=1, sticky="ns")
    tree.configure(yscrollcommand=scrollbar_y.set)

    scrollbar_x = ttk.Scrollbar(
        table_frame, orient="horizontal", command=tree.xview, style="Dark.Horizontal.TScrollbar"
    )
    scrollbar_x.grid(row=1, column=0, sticky="ew")
    tree.configure(xscrollcommand=scrollbar_x.set)

    table_frame.rowconfigure(0, weight=1)
    table_frame.columnconfigure(0, weight=1)

    columns = [df.index.name if df.index.name else "Index"] + list(df.columns)
    tree["columns"] = columns
    tree["show"] = "headings"

    for col in columns:
        heading_text = "Run" if col == "Run" else display_names.get(col, col)
        width = 85 if col == "Run" else 155
        tree.heading(col, text=heading_text)
        tree.column(col, anchor="center", width=width, stretch=True)

    tree.tag_configure("oddrow", background="#08111f", foreground="#7fd2ff")
    tree.tag_configure("evenrow", background="#0d1727", foreground="#7fd2ff")

    for row_number, (idx, row) in enumerate(df.iterrows()):
        tag = "evenrow" if row_number % 2 == 0 else "oddrow"
        tree.insert("", "end", values=[idx] + list(row), tags=(tag,))


root = tk.Tk()
root.withdraw()  # hide main empty window

show_dataframe_popup(coded_df, "Coded Matrix (-1 / +1)")
show_dataframe_popup(real_df, "Real Value Matrix")

root.mainloop()
