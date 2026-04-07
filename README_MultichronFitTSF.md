# MultichronFitTSF

**Predicting bedrock age–elevation transects from multi-chronometer detrital thermochronology data**

A two-script MATLAB toolkit that translates detrital thermochronologic datasets into explicit, spatially referenced age–elevation constraints for thermal-kinematic models such as Pecube and A2E. Developed by J. Giblin, Arizona State University.

If you use this code, please cite:
> Giblin, J. et al. (in prep). 
> Gallagher, K., & Parra, M. (2020). A new approach to thermal history modelling with detrital thermochronological data. *Earth and Planetary Science Letters*, 529, 115872.

---

## Overview

Detrital thermochronologic grains collected at a catchment outlet have no explicit source elevation label. To reconstruct a bedrock age–elevation profile A(z) from such data, this toolkit uses the catchment hypsometry as a **topographic sampling function (TSF)** — the assumption that the probability of a grain being sourced from a given elevation is proportional to the fractional catchment area at that elevation (Gallagher & Parra, 2020).

For a candidate age–elevation function A(z) and dispersion parameter τ, the predicted probability of observing a grain age *a* is:

**p(a_i) = Σ_k [ p(z_k) · N(a_i | A(z_k), σ_i² + τ²) ]**

where p(z_k) is the hypsometric weight of elevation bin k, σ_i is the analytical uncertainty of grain i, and τ absorbs unresolved scatter from kinetic variability, sediment mixing, and source heterogeneity.

This differs from QTQt's detrital implementation (Gallagher & Parra, 2020) in that it solves directly for a statistically optimal age–elevation transect rather than inverting for a full thermal history. It is designed as a controlled intermediate step for incorporating detrital datasets into Pecube-style forward models.

---

## Two-script workflow

| Script | What it does |
|--------|-------------|
| `MultichronFitTSF.m` | **Step 1.** Fits a monotonic bedrock age–elevation curve A(z) for each chronometer using a hypsometry-weighted Gaussian mixture likelihood. All chronometers are optimized jointly with a closure-temperature-aware ordering penalty. |
| `MultichronFitTSF_Georef.m` | **Step 2.** Reads a DEM and flow accumulation raster to assign lat/lon and projected coordinates to each elevation bin in the A(z) output. Produces spatially referenced transect CSVs ready for Pecube input. |

Both scripts use the same `catchment_name` / `base_dir` two-line convention and operate on the same catchment subfolder.

---

## Key features (Step 1 — MultichronFitTSF)

- All chronometers optimized **jointly** using a quasi-Newton solver (`fminunc`), not independently
- A **closure-temperature-scaled ordering penalty** enforces the physical requirement that lower-Tc systems yield younger ages than higher-Tc systems, without requiring explicit kinetic models
- **Bootstrap uncertainty** propagated jointly across all chronometers (each resample reruns the full joint solve)
- Per-grain **posterior source elevation distributions** computed for each chronometer
- Config-file driven: **only two lines change** between catchment runs

---

## Requirements

**Step 1 (`MultichronFitTSF.m`):**
- MATLAB R2019b or later
- Optimization Toolbox (for `fminunc`)

**Step 2 (`MultichronFitTSF_Georef.m`):**
- MATLAB Mapping Toolbox (R2020b+ recommended for `readgeoraster`; falls back to `geotiffread` for older versions)
- DEM must have embedded CRS metadata for true WGS84 lat/lon output

---

## Repository structure

```
MultichronFitTSF/
├── MultichronFitTSF.m           ← Step 1: fit age-elevation transects
├── MultichronFitTSF_Georef.m    ← Step 2: georeference transects using DEM
├── README.md
├── WP/                          ← example catchment subfolder
│   ├── WP_config.csv            ← all catchment-specific settings
│   ├── WP_Hypso.csv             ← catchment hypsometry
│   ├── WP_ApHe.csv              ← detrital grain ages and errors
│   ├── WP_ZHe.csv
│   ├── WP_ApPb.csv
│   ├── WP_Hbl.csv
│   ├── WP_DEM.tif               ← clipped DEM GeoTIFF (Step 2)
│   ├── WP_flowacc.tif           ← flow accumulation GeoTIFF (Step 2)
│   └── figures_svg/             ← created automatically on first run
└── example/
    └── EX/                      ← minimal working example
```

Each catchment has its own subfolder. **Only two lines in each script change between catchment runs** (`catchment_name` and `base_dir`).

---

## Quick start

### Step 1: Fit age–elevation transects

1. Clone or download this repository
2. Create a subfolder for your catchment (e.g. `WP/`)
3. Place your hypsometry CSV, grain data CSVs, and config CSV in that subfolder
4. Open `MultichronFitTSF.m` and set:
   ```matlab
   catchment_name = "WP";
   base_dir       = "/path/to/your/project/folder";
   ```
5. Run the script. All outputs are written into the catchment subfolder.

A minimal working example with synthetic data is provided in `example/EX/`.

### Step 2: Georeference transects

1. Place your clipped DEM and flow accumulation GeoTIFFs in the catchment subfolder, named `<catchment>_DEM.tif` and `<catchment>_flowacc.tif` (or override paths in the script)
2. Open `MultichronFitTSF_Georef.m` and set the same `catchment_name` and `base_dir`
3. Set `chron_names` to match the chronometers you ran in Step 1
4. Run the script. Georeferenced CSVs and a diagnostic map are written to `<catchment>/georef_outputs/`

**Tip:** Always inspect the diagnostic map before using georeferenced outputs as Pecube input. Check that the channel network looks physically reasonable and that representative points span the full elevation range of the catchment.

---

## Input file formats

### Hypsometry file (`<catchment>_Hypso.csv`)

| Column | Description |
|--------|-------------|
| `Elevation_m` or `Elevation` | Elevation in meters (any order; auto-sorted) |
| `RelArea` | Cumulative relative area (0–1); used directly as CDF |
| `Area` | Cumulative area (any units); normalized to CDF internally |

The script auto-detects which CDF column is present.

### Grain data files (`<catchment>_<Chron>.csv`)

One file per chronometer. Must contain:

| Column | Description |
|--------|-------------|
| `Date_Ma` | Grain age in Ma |
| `Error_Ma` | 1σ analytical uncertainty in Ma |

Rows with NaN ages, NaN errors, or non-positive errors are automatically excluded.

### Config file (`<catchment>_config.csv`)

The main file you edit between catchments. One row per chronometer plus one row for the hypsometry.

> **Using fewer than 4 chronometers?** No code changes needed. Simply include only the chronometers you have and leave the rest out. The joint solver automatically builds ordering constraints from whatever is present — 1 pair for 2 chronometers, 3 pairs for 3, 6 pairs for 4. With only one chronometer the script runs normally with no ordering penalty applied.

```
Chronometer,File,TauMin,AgeMargin,Lambda,AgeMinFilter,AgeMaxFilter
Hypsometry,WP_Hypso.csv,,,,,
ApHe,WP_ApHe.csv,0.5,5,0.0,0,Inf
ZHe,WP_ZHe.csv,0.5,5,0.0,0,Inf
ApPb,WP_ApPb.csv,0.1,15,0.5,0,98
Hbl,WP_Hbl.csv,0.1,5,1.0,83,Inf
```

**Optional columns on the Hypsometry row** (override global defaults):
- `w_order`: ordering penalty weight (default 5.0)
- `delta_min`: minimum age separation scaling in Ma (default 1.0)

#### Config column descriptions

| Column | Description |
|--------|-------------|
| `TauMin` | Floor on dispersion parameter τ (Ma). Use 0.5 for AHe/ZHe; 0.1 for ApPb/HblAr |
| `AgeMargin` | Padding beyond min/max observed age for A(z) search bounds (Ma). Use 5 for most; 15 for ApPb |
| `Lambda` | Curvature smoothness regularization weight. Use 0.0 for AHe/ZHe; 0.5 for ApPb; 1.0–2.0 for HblAr |
| `AgeMinFilter` | Exclude grains younger than this (Ma). Use 0 for none |
| `AgeMaxFilter` | Exclude grains older than this (Ma). Use Inf for none |

Filters should be applied with geological justification only (e.g., to exclude grains from older magmatic sources). All filtering decisions should be documented in your methods.

---

## Closure temperatures

Hardcoded defaults from Hodges (2014) Table 2:

| Chronometer | Tc (°C) |
|-------------|---------|
| AHe | 70 |
| ZHe | 170 |
| ApPb | 460 |
| HblAr | 570 |

Edit the `TC_DEFAULTS` struct at the top of `MultichronFitTSF.m` to override for non-standard systems.

---

## Outputs

### Step 1 outputs (written to `<base_dir>/<catchment_name>/`)

| File | Description |
|------|-------------|
| `predicted_bedrock_transect_<Chron>.csv` | Best-fit A(z): elevation, hypsometric weight, predicted age per bin |
| `predicted_bedrock_transect_<Chron>_CI.csv` | Same plus bootstrap median, 16th/84th percentile CI, ±1σ columns |
| `grain_expected_source_<Chron>.csv` | Per-grain posterior source elevation (mean, median, P05–P95, SD) |
| `grain_posteriors_<Chron>.csv` | Full posterior matrix P(z\|age) — one column per grain, one row per elevation bin |
| `summary_fit_params.csv` | One row per chronometer: grain count, NLL, τ, A_min, A_max, all settings |
| `ordering_violations_preFit.csv` | Pre-fit ordering violation log |
| `ordering_penalty_contributions.csv` | Post-fit per-pair penalty contributions |
| `figures_svg/` | SVG diagnostic figures: observed vs predicted PDF, A(z) curve, hypsometry CDF comparison, grain source elevation scatter, all-chronometer joint plot |

### Step 2 outputs (written to `<base_dir>/<catchment_name>/georef_outputs/`)

| File | Description |
|------|-------------|
| `predicted_bedrock_transect_<Chron>_georef.csv` | Transect with Lat, Lon, Easting, Northing, Channel_Elev_m appended |
| `grain_expected_source_<Chron>_georef.csv` | Per-grain source with median and P16/P84 coordinates appended |
| `<catchment>_dem_coord_assignment.svg/pdf/png` | Diagnostic map showing channel network and representative bin points |

---

## Interpreting results

**τ (dispersion parameter)** is the most informative single output number from Step 1.
- τ near τ_min: data are internally consistent; elevation structure explains most age spread.
- τ much larger than analytical errors: real source complexity (multiple lithologies, magmatic pulses, kinetic scatter). This is a geological signal, not a model failure — interpret it.
- τ exactly at the floor: the optimizer wanted to go lower. Check whether the predicted PDF looks too narrow.

**A(z) shape:** A gently increasing monotonic profile is expected. A completely flat profile means the chronometer has no resolvable elevation signal in this catchment (common for high-Tc systems in rapidly exhumed terranes). A staircase pattern is a normal consequence of the monotonic constraint where the likelihood surface is flat.

**Ordering penalty diagnostics:** Check `ordering_penalty_contributions.csv` after each run. Large residual penalties after the joint fit indicate the data are in genuine conflict with the expected Tc ordering — this is a geologically interesting result that warrants investigation.

---

## Common issues

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Profile completely flat; A_min ≈ A_max | No elevation signal, or search window too wide | Lower τ_min; check AgeMargin |
| Upper/lower bins all same age | Profile hitting AgeMargin bound | Increase AgeMargin |
| τ very large (>> analytical errors) | Source heterogeneity, outlier grains, bimodal distribution | Inspect age distribution; filter with geological justification; document elevated τ |
| Jagged/staircase A(z) | Flat likelihood surface | Increase Lambda (try 0.5, 1, 2) |
| < 10 grains after filtering | Too few data | Check filter settings; chronometer skipped automatically |
| Bootstrap CI very wide | Few grains or flat likelihood | Increase n_boot; inspect grain distribution |
| No channel pixels found (Step 2) | flow_acc_threshold too high, or DEM/flowacc extent mismatch | Lower threshold; check raster extents match |
| Lat/Lon output as projected X/Y (Step 2) | DEM GeoTIFF missing embedded CRS metadata | Re-export DEM from GIS with CRS set |

---

## Key assumptions and limitations

- **Hypsometric TSF**: Sediment production and transport efficiency are assumed uniform across elevation. Lithologic heterogeneity, glacial modification, and channel routing may violate this. Compare the implied source CDF (from grain posteriors) against the hypsometric CDF as a diagnostic.
- **Monotonicity**: A(z) is constrained non-decreasing with elevation. Appropriate for steady-state exhumation through horizontal isotherms; may be inappropriate in structurally complex settings.
- **Non-uniqueness**: Multiple A(z) functions can reproduce similar detrital distributions. Bootstrap CIs capture grain sampling uncertainty but not this fundamental non-uniqueness.
- **Single τ per chronometer**: τ is a scalar that absorbs all unresolved variance. It cannot distinguish kinetic dispersion from lithologic mixing from model mismatch.
- **Channel representative point (Step 2)**: The highest-flow-accumulation pixel at each elevation is used as the representative coordinate. This approximates the trunk stream routing path but may not be appropriate in catchments with complex drainage geometry.

---

## Contact

Jackie Giblin — jlgiblin@asu.edu | GitHub: jlgiblin
