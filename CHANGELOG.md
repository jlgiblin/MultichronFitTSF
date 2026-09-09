# Changelog

All notable changes to MultichronFitTSF are documented here.

## 0.1.0 — 2026-09-09

First stable public release.

- Jointly fits monotonic age–elevation transects for a user-configured set of
  chronometers.
- Requires an explicit choice between fixed and iterative source weighting.
- Supports user-defined closure temperatures and optional source-weight
  groups.
- Provides reproducible serial or parallel bootstrap uncertainty estimates.
- Includes guarded iterative source-weight estimation and convergence
  diagnostics.
- Includes optional DEM-based coordinate assignment with explicit QA fields.
- Aligns the included synthetic example's folder and filename prefixes so it
  follows the documented two-line run convention directly.
- Corrects fixed-source-weight bootstrap initialization so uncertainty runs
  work in both supported source-weighting modes.
