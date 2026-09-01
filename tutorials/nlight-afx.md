---
layout: versioned
title: nLight AFX
parent: Tutorials
nav_order: 4
permalink: /tutorials/nlight-afx/
---

# nLight AFX Heat Source

The SS316L case uses the ORNL-characterized nLight AFX Index 6 profile with a `projected` heat source and `exponential` axial projection.

<figure class="documentation-figure">
  <video controls autoplay loop muted playsinline poster="{{ '/assets/images/visualizations/nlight-afx-temperature.png' | relative_url }}" aria-label="Animated plan view of the nLight AFX SS316L temperature field and phase boundaries">
    <source src="{{ '/assets/images/visualizations/nlight-afx-temperature.mp4' | relative_url }}" type="video/mp4">
    <img src="{{ '/assets/images/visualizations/nlight-afx-temperature.png' | relative_url }}" alt="Plan view of the nLight AFX SS316L temperature field with solidus and liquidus contours.">
  </video>
  <figcaption>SS316L top-surface temperature for the nLight AFX profile. The black and white contours denote the 1471 K solidus and 1709 K liquidus, respectively.</figcaption>
</figure>

## Physical setup

- SS316L properties from `$ADDITIVEFOAM_ETC/materials/SS316L.cfg`.
- One `projected` heat source with an `nLightAFX` planar profile.
- Seven characterized beam modes, `Index0` through `Index6`; the tutorial selects `Index6`.
- An `exponential` axial projection with configured `depth 20` µm.
- A liquidus reference depth; the configured depth is the minimum applied depth.
- Eight MPI ranks.

## Run

```bash
cp -r "$ADDITIVEFOAM_TUTORIALS/nLightAFX" "$FOAM_RUN/nLightAFX"
cd "$FOAM_RUN/nLightAFX"
./Allrun
```

## Important inputs

| File | Purpose |
|---|---|
| `$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg` | Shared definitions for characterized modes `Index0` through `Index6` |
| `constant/heatSourceDict` | Profile, projection, absorption, integration, and AMR settings |
| `constant/scanPath` | Beam position, time, and power |
| `constant/transportProperties` | Includes the SS316L material properties |
| `system/decomposeParDict` | Eight-rank domain decomposition |
| `system/controlDict` | Run controls and optional Function Objects |

## Workflow

The heat-source dictionary includes the shared mode definitions and selects one inside the nested `profile` dictionary:

```foam
#include "$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg"

sources
{
    beam
    {
        path            scanPath;
        widthReference  D4Sigma;
        depthReference  isotherm;

        absorption
        {
            model   constant;
            eta     0.33;
        }

        heatSource
        {
            model       projected;
            depth       20.0e-6;
            nPoints     (10 10 10);

            profile
            {
                model   nLightAFX;
                $Index6;
            }

            projection
            {
                model       exponential;
                nSlope      0.0;
                nIntercept  1.0;
            }
        }
    }
}
```

Change only `$Index6` to `$Index0`, `$Index1`, …, or `$Index5` to select another characterized mode.

Each index supplies:

```foam
alpha
r0
sigma0
r1
sigma1
```

`r0`, `sigma0`, `r1`, and `sigma1` are metres. `alpha` is the integrated power fraction assigned to the outer ring. The tutorial uses `nSlope 0`, `nIntercept 1`, and the area-equivalent `D4Sigma` width. See [Heat Source Models]({{ '/docs/heat-sources/#nlightafx-profile' | relative_url }}) for the profile and projection equations.

## Outputs

Open the reconstructed case in ParaView and visualize `T`, `alpha.solid`, and `qDot` to compare the selected ring profile with its temperature and phase response.

Absorbed power, melt-pool dimensions, solidification data, and ExaCA histories are covered by the [baseline output workflows]({{ '/tutorials/amb2018/#outputs' | relative_url }}).
