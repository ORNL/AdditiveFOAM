---
layout: versioned
title: Tutorials
nav_order: 5
has_children: true
permalink: /tutorials/
---

# Tutorials

Copy tutorials from `$ADDITIVEFOAM_TUTORIALS` into `$FOAM_RUN`; their `Allrun` scripts create processor and result directories in place.

| Tutorial | Material | Demonstrates |
|---|---|---|
| [AMB2018-02-B]({{ '/tutorials/amb2018/' | relative_url }}) | IN625 | Calibrated single track and sub-cell source integration |
| [Multi-beam]({{ '/tutorials/multi-beam/' | relative_url }}) | IN625 | Two simultaneous moving sources |
| [Multi-layer PBF]({{ '/tutorials/multi-layer/' | relative_url }}) | IN625 | Powder deposition, layer mapping, dynamic source depth, and ExaCA |
| [nLight AFX]({{ '/tutorials/nlight-afx/' | relative_url }}) | SS316L | Characterized ring-beam modes |
| [Tabulated profile]({{ '/tutorials/tabulated/' | relative_url }}) | AlSi10Mg | Measured PRIMES profile conversion and interpolation |
| [Heat-source calibration]({{ '/tutorials/heat-source-calibration/' | relative_url }}) | SS316L | Bayesian calibration of the projected depth distribution |

## Common workflow

Single-layer tutorials run `blockMesh`, `decomposePar`, `additiveFoam`, `reconstructPar`, and the optional Function Object reconstruction helpers. The multi-layer tutorial uses `createScanPath`, `runLayers`, and `reconstructLayers` instead.

The heat-source calibration tutorial is a campaign rather than one case. `calibrateHeatSource` renders and runs a trial-case grid, extracts liquidus melt-pool depths, fits the projected-source closure, and writes its own report and cache files.

## Post-processing workflows

The AMB2018-02-B tutorial provides complete workflows for generating and processing the standard quantitative and microstructure outputs:

- [Plot absorbed power]({{ '/tutorials/amb2018/#plot-absorbed-power' | relative_url }})
- [Write and plot melt-pool dimensions]({{ '/tutorials/amb2018/#write-and-plot-melt-pool-dimensions' | relative_url }})
- [Write solidification data and plot CET curves]({{ '/tutorials/amb2018/#write-solidification-data-and-plot-cet-curves' | relative_url }})
- [Create ExaCA temperature histories]({{ '/tutorials/amb2018/#create-exaca-temperature-histories' | relative_url }})
