# 3D Printing, Density, and SEM Analysis

This project contains the experimental designs, measurements, analysis tools,
and figures used to study PLA and copper-filled material extrusion. The work
connects printing parameters to mass, volume, density, relative density,
surface quality, sintering behavior, and SEM-observed microstructure.

## What I did in this project

The project covers five related areas:

1. **PLA design of experiments**
   - Built fractional-factorial screening designs for printing parameters.
   - Compared rectilinear and concentric infill patterns.
   - Evaluated mass, hydrostatic volume, density, repeatability, and visual
     defects.
   - Added Box-Behnken optimization and separate ironing experiments.

2. **Copper-filament design of experiments**
   - Analyzed Virtual Foundry VF1 with a definitive screening design.
   - Analyzed Virtual Foundry VF2 with a factorial design.
   - Summarized the available VF3 sintered-density measurements.
   - Built exploratory Markforged mass and density response models.
   - Produced main-effect, interaction, Pareto, ANOVA, residual, contour,
     correlation, and desirability plots.

3. **Green and sintered density**
   - Calculated Archimedes density from repeated measurements.
   - Compared green VF1, VF2, and VF3 filament-reference densities.
   - Compared sintered reference density and relative density across thermal
     cycles, including measurement-specific error bars.

4. **SEM and EDS analysis**
   - Analyzed pore measurements exported from Fiji/ImageJ.
   - Calculated number-weighted and area-weighted pore distributions.
   - Added Feret X/Y spans, equivalent circular diameter, and D10/D50/D90.
   - Summarized copper-particle count, area fraction, circularity, and solidity
     in green filament cross-sections.
   - Extracted spectra and element maps from Thermo Fisher Axia ChemiSEM PAD
     projects.
   - Created semi-quantitative EDS weight-percent estimates using T3 as a
     transfer reference.

5. **Geometry and thermal-profile support**
   - Modeled printed channel geometry and inverse dimensional correction.
   - Created thermal-profile and CO/CO2 analysis helpers used alongside the
     separate oven project.

## Setup

Python 3.10 or newer is recommended. From this folder, create or activate a
virtual environment and install the analysis packages:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy pandas matplotlib scipy pyDOE jupyter
```

`tkinter` is used by the desktop tools and is normally included with the
standard Windows Python installation.

## Project structure

I did a notebook for each series: PLA series (FDI, FDII), Copper series (VF1. VF2, MF)

### PLA DoE notebooks

| File | Purpose |
| --- | --- |
| [`DoE_PLA_FD_I.ipynb`](DoE_PLA_FD_I.ipynb) | First PLA screening design, mass measurements, repeatability, hydrostatic density, and regression. |
| [`DoE_PLA_FD_II.ipynb`](DoE_PLA_FD_II.ipynb) | Second screening design and rectilinear/concentric comparisons. |
| [`DoE_PLA_Box_Behnken.ipynb`](DoE_PLA_Box_Behnken.ipynb) | Additional runs, optimization, and response-surface analysis. |
| [`DoE_PLA_Ironing.ipynb`](DoE_PLA_Ironing.ipynb) | Ironing trials, confirmation runs, and visual-quality analysis. |
| [`DoE_pla_notebook.ipynb`](DoE_pla_notebook.ipynb) | Legacy combined PLA notebook. The split notebooks above are preferred. |

PLA measurements are stored in
[`pla_doe_measurements.json`](pla_doe_measurements.json) and loaded through
`pla_doe_measurement_loader.py` and `pla_doe_measurements.py`.

### Copper DoE and density notebooks

| File | Purpose |
| --- | --- |
| [`DoE_copper_VirtualFoundry1.ipynb`](DoE_copper_VirtualFoundry1.ipynb) | VF1 definitive-screening design, repeated measurements, density, ANOVA, contours, and desirability. |
| [`DoE_copper_VirtualFoundry2.ipynb`](DoE_copper_VirtualFoundry2.ipynb) | VF2 factorial design, mass/density effects, interactions, and ANOVA. |
| [`DoE_copper_VirtualFoundry3.ipynb`](DoE_copper_VirtualFoundry3.ipynb) | Available VF3 sintered-density measurements. |
| [`DoE_copper_Markforged.ipynb`](DoE_copper_Markforged.ipynb) | Markforged mass/density analysis and exploratory quadratic contours. |
| [`Copper_filament_density_comparison.ipynb`](Copper_filament_density_comparison.ipynb) | Green and sintered filament-reference density comparison by thermal cycle. |

The shared copper fitting and plotting functions are in
[`copper_doe_analysis.py`](copper_doe_analysis.py). Factor levels and editable
measurement data are stored in [`copper_doe_levels.json`](copper_doe_levels.json)
and [`copper_doe_measurements.json`](copper_doe_measurements.json).

### Shared DoE tools

| File or folder | Purpose |
| --- | --- |
| [`ffd_utils.py`](ffd_utils.py) | Design generation, mapping coded levels, density calculations, statistics, and plots. |
| [`FFD.py`](FFD.py) | Small desktop viewer for an earlier fractional-factorial design. |
| [`DoE_old_scale.py`](DoE_old_scale.py) | Archived measurements from the old scale; kept separate from current analysis. |
| [`doe_matrices`](doe_matrices/) | Exported design matrices. |
| [`ffd_data`](ffd_data/) | Exported measurements and supporting tables. |

### SEM and EDS tools

| File or folder | Purpose |
| --- | --- |
| [`SEM_analysis.ipynb`](SEM_analysis.ipynb) | Main combined pore, particle, and EDS reporting notebook. |
| [`sem_pore_analysis.py`](sem_pore_analysis.py) | Fiji pore-data validation, geometry statistics, and histograms. |
| [`sem_particle_analysis.py`](sem_particle_analysis.py) | Green-filament copper-particle summary with D10, D50, and D90. |
| [`sem_eds_analysis.py`](sem_eds_analysis.py) | Extracts PAD acquisition metadata, summed spectra, maps, and montages. |
| [`sem_eds_wt_estimates.py`](sem_eds_wt_estimates.py) | Creates T3-reference-calibrated EDS estimates. |
| [`SEM_analysis_outputs`](SEM_analysis_outputs/) | Input CSV/PAD files and generated specimen-specific outputs. |

### Other tools

| File | Purpose |
| --- | --- |
| [`channel_geometry.ipynb`](channel_geometry.ipynb) | Channel calibration and inverse geometry recommendation. |
| [`notebook_plot_saver.py`](notebook_plot_saver.py) | Automatically saves Matplotlib figures shown by the notebooks. |
| [`oven_temperature_profile_tool.py`](oven_temperature_profile_tool.py) | Desktop two-stage debinding/sintering profile generator. |
| [`oven_parameters.py`](oven_parameters.py) | Temperature-profile and PLA mass-loss helpers. |
| [`oven_CO.py`](oven_CO.py) | CO/CO2 throughput conversions. |
| [`literature_visualization.py`](literature_visualization.py) | Literature-based process visualization. |

## How to use the DoE notebooks

Start Jupyter from this project folder:

```powershell
python -m jupyter lab
```

Open the required notebook and run it from the first setup cell to the end.
Running only the final plotting cell can leave old variables or measurements in
memory.

When adding measurements:

1. Update the appropriate JSON file or the clearly marked measurement cell.
2. Keep unavailable results empty or `NaN`; do not enter zero for missing data.
3. Check the displayed measurement table and units.
4. Run the notebook from the top so every model and plot is recalculated.
5. Review residuals, model rank, repeated measurements, and data coverage before
   interpreting effect sizes.

Figures are automatically written to [`notebook_figures`](notebook_figures/).
Its folder structure is documented in
[`notebook_figures/README.md`](notebook_figures/README.md).

## How to analyze SEM pores

Export Fiji/ImageJ `Analyze Particles` results as CSV. The required pore-size
column is `Feret`. Extended geometry uses `Area`, `MinFeret`, `FeretX`,
`FeretY`, and `FeretAngle` when available.

Example for M08:

```powershell
python sem_pore_analysis.py `
  ".\SEM_analysis_outputs\M08\results_all.csv" `
  --output-dir ".\SEM_analysis_outputs\M08" `
  --sample-name M08 `
  --input-unit um
```

Show all available options with:

```powershell
python sem_pore_analysis.py --help
```

Use the complete analyzed image area when reporting pore area fraction. The
number-weighted median describes the middle detected pore by count. The
area-weighted distribution shows which pore sizes contribute most of the total
projected pore area.

The green-filament particle table is created with
`create_green_filament_particle_summary()` in `sem_particle_analysis.py` and is
also used by `SEM_analysis.ipynb`. D10, D50, and D90 are number-weighted
percentiles of maximum Feret diameter.

## How to analyze EDS PAD files

Extract one PAD project:

```powershell
python sem_eds_analysis.py `
  ".\SEM_analysis_outputs\T3\T3.pad" `
  --output-dir ".\SEM_analysis_outputs\T3\EDS_analysis"
```

The exact PAD filename may differ between specimen folders. Pass the real path
shown in the folder.

Create combined estimates for all unique PAD datasets found below
`SEM_analysis_outputs`:

```powershell
python sem_eds_wt_estimates.py
```

This scans the specimen folders, detects duplicated PAD datasets, uses the
supplied T3 composition as a transfer reference, and writes per-specimen and
combined CSV/PNG tables.

The estimated EDS wt% values are semi-quantitative. They are not a replacement
for vendor ZAF or Phi-Rho-Z corrected quantification. Carbon can include
binder, coating, mounting material, or contamination. A displayed `n.d.` means
the element was not reliably determined; it does not mean an exact zero wt%.

## Density conventions

- Archimedes density must be calculated from the individual measurement values
  before averaging densities.
- Error bars represent the spread of the measurements belonging to that bar.
- Relative density is calculated against the theoretical density configured in
  the notebook.
- A result above 100% indicates measurement/input uncertainty or an unsuitable
  reference density and must not be interpreted as physically pore-free.
- Keep air weight, water/displaced-volume value, density, and units together so
  corrected measurements can be traced.

## Generated outputs

- `notebook_figures/`: notebook plots grouped by experiment and response;
- `SEM_analysis_outputs/<specimen>/`: pore and particle summaries;
- `SEM_analysis_outputs/<specimen>/EDS_analysis/`: spectra, maps, diagnostics,
  and estimated-composition tables;
- `doe_matrices/`: exported experimental designs;
- `ffd_data/`: exported measurement tables; and
- `oven_outputs/`: generated temperature-profile figures.

Raw instrument exports, PAD files, and measurement JSON files should be kept as
source data. Regenerated figures and tables belong in the existing output
folders rather than beside the analysis scripts.
