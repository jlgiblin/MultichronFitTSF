# MultichronFitTSF

**Predicting bedrock age–elevation transects from multi-chronometer detrital thermochronology data**

A two-script MATLAB toolkit that translates detrital thermochronologic datasets into explicit, spatially referenced age–elevation constraints for thermal-kinematic models such as Pecube and A2E. Developed by J. Giblin, Arizona State University.

If you use this code, please cite:
> Gallagher, K., & Parra, M. (2020). A new approach to thermal history modelling with detrital thermochronological data. *Earth and Planetary Science Letters*, 529, 115872.

---

## Overview

Detrital thermochronologic grains collected at a catchment outlet have no explicit source elevation label. To reconstruct a bedrock age–elevation profile A(z) from such data, this toolkit uses the catchment hypsometry as a **topographic sampling function (TSF)** — the assumption that the probability of a grain being sourced from a given elevation is proportional to the fractional catchment area at that elevation (Gallagher & Parra, 2020).

For a candidate age–elevation function A(z) and dispersion parameter τ, the predicted probability of observing a grain age *a* is:

**p(a_i) = Σ_k [ p(z_k) · N(a_i | A(z_k), σ_i² + τ²) ]**

where p(z_k) is the source weight of elevation bin k, σ_i is the analytical uncertainty of grain i, and τ absorbs unresolved scatter from kinetic variability, sediment mixing, and source heterogeneity. The selected source-weighting mode determines p(z_k).

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
- A clear **fixed or iterative source-weighting choice** analogous to the
  fixed/flexible distinction in detrital QTQt workflows
- Explicit, user-defined **source-weight groups**; the code does not infer
  mineral behavior from chronometer names
- Iterative-mode **bootstrap uncertainty** re-estimates the configured source
  weights inside every grain resample
- Optional **parallel bootstrap execution** with reproducible per-resample
  random streams and automatic serial fallback
- Config-file driven: **only two lines change** between catchment runs

---

## Requirements

**Step 1 (`MultichronFitTSF.m`):**
- MATLAB R2019b or later
- Optimization Toolbox (for `fminunc`)
- Parallel Computing Toolbox (optional, for parallel bootstrap execution)

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
├── SampleA/                     ← example catchment subfolder
│   ├── SampleA_config.csv       ← all catchment-specific settings
│   ├── SampleA_Hypsometry.csv   ← catchment hypsometry
│   ├── SampleA_ApHe.csv         ← detrital grain ages and errors
│   ├── SampleA_ZHe.csv
│   ├── SampleA_ApPb.csv
│   ├── SampleA_Hbl.csv
│   ├── SampleA_DEM.tif          ← clipped DEM GeoTIFF (Step 2)
│   ├── SampleA_flowacc.tif      ← flow accumulation GeoTIFF (Step 2)
│   └── figures_svg/             ← created automatically on first run
└── Example/                     ← minimal working example files (EX_*)
```

Each catchment has its own subfolder. **Only two lines in each script change between catchment runs** (`catchment_name` and `base_dir`).

---

## Quick start

### Step 1: Fit age–elevation transects

1. Clone or download this repository
2. Create a subfolder for your catchment (e.g. `SampleA/`)
3. Place your hypsometry CSV, grain data CSVs, and config CSV in that subfolder
4. Open `MultichronFitTSF.m` and set:
   ```matlab
   catchment_name = "SampleA";
   base_dir       = "/path/to/your/project/folder";
   ```
5. Run the script. All outputs are written into the catchment subfolder.

The default first run leaves bootstrap off so the fit and source-weight
convergence can be checked quickly. Then set `do_bootstrap = true` and use a
small pilot (for example, `n_boot = 20`) before a final uncertainty run.

A minimal working example with synthetic data is provided in `Example/`.

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

> **Using fewer than 4 chronometers?** No code changes are needed. Include only the chronometers you have. The solver constructs ordering constraints between rows with distinct, finite closure temperatures. With one chronometer—or with closure temperatures left blank—the script runs without the corresponding ordering penalties.

```csv
Chronometer,File,TauMin,AgeMargin,Lambda,AgeMinFilter,AgeMaxFilter,ClosureTemperature_C,TSFMode,TSFGroup,EstimateTSF
Hypsometry,SampleA_Hypsometry.csv,,,,,,,iterative,,
ApHe,SampleA_ApHe.csv,0.5,5,0.0,0,Inf,70,,apatite,true
ZHe,SampleA_ZHe.csv,0.5,5,0.0,0,Inf,170,,,true
ApPb,SampleA_ApPb.csv,0.1,15,0.5,0,Inf,460,,apatite,true
Hbl,SampleA_Hbl.csv,0.1,5,1.0,0,Inf,570,,hornblende,false
```

`TSFMode` is required on the Hypsometry row and has no default. Choose
`fixed` to retain the measured hypsometric weights throughout the fit, or
`iterative` to start from measured hypsometry and repeatedly estimate
regularized effective source weights while refitting the transects. Iterative
mode can produce different results and requires more computing time, so the
choice should reflect the user's question and be reported with the results.

**Optional columns on the Hypsometry row** (override numerical defaults):
- `w_order`: ordering penalty weight (default 5.0)
- `delta_min`: minimum age separation scaling in Ma (default 1.0)
- `TSFUpdateFraction`: per-iteration relaxation step toward the newly
  estimated source weights (default 0.4; 0 = no update, 1 = full update)
- `TSFSmoothSpan`: moving-mean span in equal-area bins (default 3; 1 = none)

**Optional chronometer-row columns:**

- `ClosureTemperature_C`: nominal closure temperature used to construct the
  ordering constraints. It does not otherwise change the age model.
- `TSFGroup`: matching labels share one effective source-weight curve in
  iterative mode. A blank group gives that chronometer an independent curve.
  Labels are entirely user-defined and case-insensitive.
- `EstimateTSF`: `true` estimates the group's weights; `false` holds that
  group to the measured hypsometry. Blank defaults to `true`. All rows in one
  named group must agree.

These columns are explicit because sharing weights is a scientific choice,
not something the code should infer from a system abbreviation. For example,
ApHe and ApPb may share an `apatite` group, while a blank ZHe group remains
independent. A group with little resolvable elevation structure can be held
to measured hypsometry with `EstimateTSF=false`. Other groupings and
chronometers work without editing the source code.

For an initial iterative run, keep bootstrap disabled until the convergence
history has been inspected:

```matlab
tsf_update_fraction         = 0.4;
tsf_smooth_span             = 3;
iterative_max_outer          = 10;
iterative_hypsometry_pull    = 0.25;
iterative_smoothness         = 1.0;
iterative_no_improve_patience = 3;
do_bootstrap                = false;
```

Because the source-weight density objective is not identical to the grain likelihood,
the solver retains the lowest-NLL iteration and stops after the configured
number of consecutive likelihood deteriorations. The selected iteration is
identified in `source_weighting_convergence.csv`.

After the unbootstrapped iterative run has been checked for convergence, use
a small uncertainty pilot before a final run:

```matlab
iterative_max_outer         = 40;  % main-data solution
iterative_boot_max_outer    = 20;  % cap within each resample
do_bootstrap               = true;
n_boot                     = 20;
bootstrap_random_seed      = 1;   % repeatable resampling
use_parallel_bootstrap     = true;
parallel_worker_count      = 4;   % reduce this on memory-limited systems
```

For a final analysis, `n_boot = 200` is a practical starting point. More
resamples may be useful when interval bounds remain unstable, but they are not
automatically more appropriate; report the value used and check that the
resulting confidence intervals are sufficiently stable for the application.

Every resample first obtains its own fixed-hypsometry fit and then performs a
guarded iterative source-weight inversion. This is intentionally more
expensive than a fixed bootstrap. Inspect
`source_weighting_bootstrap_summary.csv` before increasing `n_boot`. Keep the
same `bootstrap_random_seed` to reproduce a run exactly, or record a different
nonnegative integer when generating an independent set of resamples. Parallel
mode uses the requested process-worker count when it starts a pool. If a pool
is already open, the code uses that pool; if the toolbox or pool is unavailable,
it reports the issue and continues serially.

#### Config column descriptions

| Column | Description |
|--------|-------------|
| `TauMin` | Floor on dispersion parameter τ (Ma). Use 0.5 for AHe/ZHe; 0.1 for ApPb/HblAr |
| `AgeMargin` | Padding beyond min/max observed age for A(z) search bounds (Ma). Use 5 for most; 15 for ApPb |
| `Lambda` | Curvature smoothness regularization weight. Use 0.0 for AHe/ZHe; 0.5 for ApPb; 1.0–2.0 for HblAr |
| `AgeMinFilter` | Exclude grains younger than this (Ma). Use 0 for none |
| `AgeMaxFilter` | Exclude grains older than this (Ma). Use Inf for none |
| `ClosureTemperature_C` | User-defined nominal closure temperature used only for ordering |
| `TSFGroup` | Optional source-weight group; blank means independent, matching labels share weights |
| `EstimateTSF` | `true` estimates the group; `false` keeps it fixed; blank defaults to `true` |

Filters should be applied with geological justification only (e.g., to exclude grains from older magmatic sources). All filtering decisions should be documented in your methods.

---

## Closure temperatures

Enter the nominal closure temperature appropriate to each dataset in
`ClosureTemperature_C`. These values determine only the order and relative
spacing of the age–elevation constraints; the code does not use them as a
diffusion model. Because closure temperature can depend on mineral kinetics,
grain properties, cooling rate, and the chosen calibration, users should
select and document values appropriate to their application.

For backward compatibility, recognized blank entries fall back to these
built-in values from Hodges (2014) Table 2:

| Chronometer | Tc (°C) |
|-------------|---------|
| AHe | 70 |
| ZHe | 170 |
| ApPb | 460 |
| HblAr | 570 |

Unrecognized chronometers with a blank closure temperature are still fitted,
but are omitted from the ordering penalty.

---

## Outputs

### Step 1 outputs (written to `<base_dir>/<catchment_name>/`)

| File | Description |
|------|-------------|
| `predicted_bedrock_transect_<Chron>.csv` | Best-fit A(z): elevation, hypsometric weight, predicted age per bin |
| `predicted_bedrock_transect_<Chron>_CI.csv` | Same plus bootstrap median, 16th/84th percentile CI, ±1σ columns |
| `source_weighting.csv` | Hypsometric and selected source weights, group membership, raw updates, relative yield, and cumulative distributions |
| `source_weighting_convergence.csv` | Iterative-mode history: NLL, weight change, distance from hypsometry, and source-weight fit SSE |
| `source_weighting_bootstrap_CI.csv` | Group-specific best-fit source weights and bootstrap median/16th/84th percentiles |
| `source_weighting_bootstrap_summary.csv` | Per-resample iterative inversion success, NLL change, selected iteration, and stopping reason |
| `grain_expected_source_<Chron>.csv` | Per-grain posterior source elevation (mean, median, P05–P95, SD) |
| `grain_posteriors_<Chron>.csv` | Full posterior matrix P(z\|age) — one column per grain, one row per elevation bin |
| `summary_fit_params.csv` | One row per chronometer: grain count, NLL, τ, A_min, A_max, all settings |
| `ordering_violations_preFit.csv` | Pre-fit ordering violation log |
| `ordering_penalty_contributions.csv` | Post-fit per-pair penalty contributions |
| `figures_svg/` | SVG diagnostic figures for the selected mode: observed vs predicted PDF, A(z) curve, source-weighting CDF, grain source elevation scatter, and all-chronometer joint plot |

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

**Iterative source-weight diagnostics:** `RelativeYieldVsHypsometry = 1` means
the selected effective source weight matches area-proportional sourcing in
that bin. Values above or below 1 indicate relative over- or
under-representation. Iterative bootstrap runs re-estimate every group marked
`EstimateTSF=true` within each grain resample.

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
| Source weights shift but NLL does not improve | Source weights are weakly resolved, often because the age profile is narrow or flat | Keep that group fixed with `EstimateTSF=false`; compare fixed and iterative diagnostics |
| No channel pixels found (Step 2) | flow_acc_threshold too high, or DEM/flowacc extent mismatch | Lower threshold; check raster extents match |
| Lat/Lon output as projected X/Y (Step 2) | DEM GeoTIFF missing embedded CRS metadata | Re-export DEM from GIS with CRS set |

---

## Key assumptions and limitations

- **Hypsometric source weighting**: In fixed mode, sediment production,
  mineral fertility, preservation, and transport efficiency are assumed
  uniform per unit catchment area. Iterative mode estimates regularized
  effective source weights by elevation. It does not change the measured
  hypsometry and is not a reproduction of QTQt's transdimensional
  thermal-history inversion.
- **Monotonicity**: A(z) is constrained non-decreasing with elevation. It may
  be nonlinear and does not require a constant exhumation rate through time,
  but it may be inappropriate in structurally complex catchments.
- **Non-uniqueness**: Multiple A(z) functions can reproduce similar detrital distributions. Allowing both A(z) and the TSF to vary increases this tradeoff. Bootstrap CIs capture grain sampling uncertainty but not this fundamental non-uniqueness.
- **Single τ per chronometer**: τ is a scalar that absorbs all unresolved variance. It cannot distinguish kinetic dispersion from lithologic mixing from model mismatch.
- **Channel representative point (Step 2)**: The highest-flow-accumulation pixel at each elevation is used as the representative coordinate. This approximates the trunk stream routing path but may not be appropriate in catchments with complex drainage geometry.

---

## Contact

Jackie Giblin — jlgiblin@asu.edu | GitHub: jlgiblin
