---
layout: versioned
title: Multi-Beam
parent: Tutorials
nav_order: 2
permalink: /tutorials/multi-beam/
---

# Multi-Beam

The `multiBeam` case applies two simultaneous IN625 heat sources separated by 100 µm in the hatch direction.

<figure class="documentation-figure">
  <video controls autoplay loop muted playsinline poster="{{ '/assets/images/visualizations/multi-beam-temperature.png' | relative_url }}" aria-label="Animated plan view of the two-beam IN625 temperature field and phase boundaries">
    <source src="{{ '/assets/images/visualizations/multi-beam-temperature.mp4' | relative_url }}" type="video/mp4">
    <img src="{{ '/assets/images/visualizations/multi-beam-temperature.png' | relative_url }}" alt="Plan view of the two-beam IN625 temperature field with solidus and liquidus contours.">
  </video>
  <figcaption>IN625 top-surface temperature during simultaneous scanning by two heat sources. The black and white contours denote the 1410 K solidus and 1620 K liquidus, respectively.</figcaption>
</figure>

## Physical setup

- IN625 properties from `$ADDITIVEFOAM_ETC/materials/IN625.cfg`.
- Two simultaneous `superGaussian` heat sources.
- Each uses `radius (85 85)` µm, `depth 30` µm, `k 2`, and `definition e2`.
- Independent paths separated by 100 µm in the hatch direction.
- Fluid flow disabled by `nOuterCorrectors 0` in the tutorial case.
- Each source uses `tolerance 1e-3` and `nPoints (10 10 10)`.
- Eight MPI ranks.

## Run

```bash
cp -r "$ADDITIVEFOAM_TUTORIALS/multiBeam" "$FOAM_RUN/multiBeam"
cd "$FOAM_RUN/multiBeam"
./Allrun
```

## Important inputs

| File | Purpose |
|---|---|
| `constant/heatSourceDict` | Defines `beam1`, `beam2`, their models, and source-specific AMR buffers |
| `constant/scanPath_1` | Position, time, and power for `beam1` |
| `constant/scanPath_2` | Position, time, and power for `beam2` |
| `constant/dynamicMeshDict` | OpenFOAM topology changer |
| `system/controlDict` | Run controls and optional Function Objects |
| `system/fvSolution` | Flow, temperature, and phase-coupling controls |
| `system/decomposeParDict` | Eight-rank domain decomposition |

## Workflow

`Allrun` generates the mesh, decomposes the case, runs `additiveFoam` in parallel, and reconstructs the fields. Each named source reads its own scan path; the time and power columns determine when it is active.

When both beams are active, `qDot` and the absorbed power reported in the solver log include both sources.

## Outputs

Open the reconstructed case in ParaView and visualize `T`, `alpha.solid`, and `qDot` to inspect the interaction between the two moving sources and their melt pools.

The AMB2018-02-B tutorial documents the quantitative and microstructure output workflows:

- [Plot absorbed power]({{ '/tutorials/amb2018/#plot-absorbed-power' | relative_url }})
- [Write and plot melt-pool dimensions]({{ '/tutorials/amb2018/#write-and-plot-melt-pool-dimensions' | relative_url }})
- [Write solidification data and plot CET curves]({{ '/tutorials/amb2018/#write-solidification-data-and-plot-cet-curves' | relative_url }})
- [Create ExaCA temperature histories]({{ '/tutorials/amb2018/#create-exaca-temperature-histories' | relative_url }})
