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
| Measured beam conversion | `primesToAdditiveFoam` |
| Multi-layer simulation | `runLayers`, `reconstructLayers` |
| Function Object reconstruction | `reconstructExaCAData`, `reconstructSolidificationData` |
| Plotting | `plotPower`, `plotDimensions`, `plotCET` |

Sourcing `etc/bashrc` adds `$ADDITIVEFOAM_BIN` to `PATH`, making these commands available in the current shell.

## `calibrateHeatSource`

Run a campaign of rendered AdditiveFOAM cases and infer the `A` and `B` coefficients used by the common projected-source depth distribution:

```bash
calibrateHeatSource --config config.yml
```

The command requires the Python packages pinned in the repository `requirements.txt`. It reads experimental depths and a complete case template, caches every simulated response, performs local Bayesian inversions and a robust global fit, and writes YAML, CSV, and PDF products. See the [calibration model and configuration reference]({{ '/docs/heat-source-calibration/' | relative_url }}) and [worked tutorial]({{ '/tutorials/heat-source-calibration/' | relative_url }}).

## `createScanPath`

Use `createScanPath` when a case requires one or more rotated raster patterns. Define the rectangular scan region, hatch spacing, rotation sequence, power, speed, and inter-track dwell in `constant/createScanPathDict`, then run `createScanPath` from the case directory. It writes `constant/scanPath_0`, `scanPath_1`, and so on for direct use by `pathName` or `runLayers`:

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

Use `primesToAdditiveFoam` before running a `tabulated` heat source when the lateral intensity was exported by PRIMES LaserDiagnosticsSoftware. The command converts the vendor CSV into AdditiveFOAM's regular-grid table and reports dimensions that can be transferred to the heat-source dictionary.

```bash
primesToAdditiveFoam input.csv constant/beamProfile.txt
```

Converts a PRIMES LaserDiagnosticsSoftware export to the `tabulated` model format. It reads grid dimensions and pixel pitches, subtracts an optional null value, clamps negative intensity, converts micrometres to metres, centres the table, and normalizes its planar integral. It reports measured radii and suggested lateral dimensions when metadata is present.

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
