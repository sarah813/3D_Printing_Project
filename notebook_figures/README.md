# DoE notebook figures

Figures are saved automatically whenever a plot is displayed in the project DoE notebooks.

## Folder order

- `01_FD_I`: first fractional factorial series
- `02_FD_II`: second fractional factorial / RR series
- `03_Optimization`: extra runs and Box-Behnken analysis
- `04_Ironing`: ironing series
- `05_Channel_geometry`: channel-dimension calibration and target-area model
- `06_Copper_VirtualFoundry2`: second copper series using a three-factor Box-Behnken design

Each series is divided by response or analysis type, such as `01_mass`, `02_volume`,
`03_density`, `04_correlations`, `05_vacuum`, `06_quality`, and
`08_combined_responses`.

Run the notebook from its setup cell to refresh the PNG files. Existing figures with
the same title are overwritten during a notebook session.
