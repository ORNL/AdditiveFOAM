---
layout: versioned
title: Utilities
parent: User Guide
nav_order: 11
permalink: /docs/utilities/
redirect_from:
  - /reference/utilities/
---

# Utilities

AdditiveFOAM includes command-line utilities for calibrating projected heat sources, generating raster scan paths, converting measured beam profiles, running and reconstructing multi-layer cases, combining parallel Function Object output, and plotting power, melt-pool dimensions, and solidification conditions.

| Workflow | Utilities |
|---|---|
| Projected-source calibration | `calibrateHeatSource` |
| Scan-path generation | `createScanPath` |
| Measured beam conversion and inspection | `primesToAdditiveFoam`, `tabulatedProfileInfo` |
| Multi-layer simulation | `runLayers`, `reconstructLayers` |
| Function Object reconstruction | `reconstructExaCAData`, `reconstructSolidificationData` |
| Plotting | `plotPower`, `plotDimensions`, `plotCET` |

Sourcing `etc/bashrc` adds `$ADDITIVEFOAM_BIN` to `PATH`, making these commands available in the current shell.

## `calibrateHeatSource`

Run a campaign of rendered AdditiveFOAM cases and infer the `nSlope` and `nIntercept` coefficients used by the `projected` source's exponential axial projection:

```bash
calibrateHeatSource --config config.yml
```

The command requires Python 3.10–3.12 and the packages pinned in the repository `requirements.txt`. It reads named beam profiles, experimental depths, and a case template; evaluates local Bayesian posteriors with deterministic numerical quadrature; stores the AdditiveFOAM response trials; performs an uncertainty-weighted global fit with a `soft_l1` loss; and writes YAML, CSV, and PDF products. See the [calibration model and configuration reference]({{ '/docs/heat-source-calibration/' | relative_url }}) and [worked tutorial]({{ '/tutorials/heat-source-calibration/' | relative_url }}).

## `createScanPath`

Use `createScanPath` when a case requires one or more rotated raster patterns. Define the rectangular scan region, hatch spacing, rotation sequence, power, speed, and inter-track dwell in `constant/createScanPathDict`, then run `createScanPath` from the case directory. It writes `constant/scanPath_0`, `scanPath_1`, and so on for direct use by a source-level `path` entry or by `runLayers`:

```foam
minPoint    (0 -1e-4);
maxPoint    (0.002 1e-4);
angle       180;
hatch       1e-4;
nRotations  2;
power       195;
speed       0.8;
dwellTime   5e-4;
biDirection true;
```

The rectangular region and hatch spacing are metres; `angle` is the rotation increment in degrees. `nRotations`, power, speed, and dwell time are required. `biDirection` defaults to `true`. The utility clips rotated raster lines to the rectangle and inserts zero-power dwell/move rows between hatches.

## `primesToAdditiveFoam`

Use `primesToAdditiveFoam` before running a `tabulated` projected profile when the planar intensity was exported by PRIMES LaserDiagnosticsSoftware. The command converts the vendor CSV into AdditiveFOAM's regular-grid table and reports its calculated beam statistics.

```bash
primesToAdditiveFoam input.csv constant/beamProfile.txt
```

The converter subtracts `Nullvalue`, clamps remaining negative intensities to zero, and treats `SNR` as diagnostic metadata rather than an intensity threshold. For an LDS table export with calculation ROI enabled at the standard 0.5 fill factor, it reconstructs the rectangular calculation ROI iteratively from the pixel centroid and x/y second moments. It then crops to the nonzero support while retaining one zero-valued perimeter node, preserves the original coordinate system, converts coordinates to metres, and normalizes the exact bilinear planar integral.

The report includes the converted centroid, principal `D4Sigma` values ordered `(major minor)`, and azimuth. When compatible rotated-moments metadata are present, it also compares the area-equivalent `D4Sigma` and the converted a, b, x, and y second-moment radii with the PRIMES values. Window-centre metadata supplies the centroid in PRIMES coordinates. The converted table retains asymmetric centroid offsets, and the ROI bounds are reconstructed from the measurement statistics. The solver uses these [`profileMetrics`]({{ '/docs/heat-sources/#beam-plane-profile-metrics' | relative_url }}).

## `tabulatedProfileInfo`

Inspect any AdditiveFOAM tabulated profile without running a solver case:

```bash
tabulatedProfileInfo constant/beamProfile.txt
```

The utility reports grid counts, x/y bounds, spacing, the exact bilinear planar integral, centroid, principal `D4Sigma` values ordered `(major minor)`, and azimuth in radians. An area-equivalent width is the geometric mean of those two reported diameters. `calibrateHeatSource` calls the utility automatically for every named profile before selecting the requested `D4Sigma` component and constructing normalized experimental and simulation depths.

## Multi-layer workflow

Use `runLayers` for sequential powder-layer deposition when each layer can be simulated as a separate case derived from the preceding layer. Run it from a prepared base case containing the mesh, material, solver dictionaries, and scan-path sequence. It creates and runs `layerN` directories; `reconstructLayers` then combines their reconstructed fields into the base case for visualization.

```bash
runLayers -nLayers 2 -layerThickness 40e-6 -nCellsPerLayer 4
reconstructLayers
```

`runLayers` cycles through `constant/scanPath_*` files (or reuses `scanPath`), creates `layerN` cases, extrudes the mesh, maps the previous layer, initializes new powder, advances absolute case time, and runs each layer in parallel. Only a single repeating scan-pattern sequence is supported. Existing `layerN` directories are replaced, so run it only in a disposable case copy.

`reconstructLayers` combines reconstructed layer results into the base case.

## Function Object Reconstruction

Run these commands after a parallel simulation that enabled `ExaCA` or `solidificationData`. They combine the per-process event files into one CSV that can be plotted or passed to another application.

```bash
reconstructExaCAData [-case path]
reconstructSolidificationData [-case path]
```

Both helpers detect `layer*/` cases, merge per-rank CSV files, remove the per-rank inputs after successful reconstruction, and exit quietly when no data exist.

## Plotting

Use the plotting utilities after the corresponding solver or Function Object output exists. Each command accepts one or more cases, allowing results from different layers or parameter studies to be overlaid in the same figure.

```bash
plotPower [case-or-log ...] [-l log.additiveFoam] [-o power.png]
plotDimensions [case ...] [-o melt_pool_dimensions.png]
plotCET [case ...] [-o CET_curve.png]
```

`plotPower` requires Matplotlib. `plotDimensions` requires Pandas and Matplotlib. `plotCET` requires Pandas, NumPy, and Matplotlib. Multiple case arguments overlay results for layer or parameter comparisons.
