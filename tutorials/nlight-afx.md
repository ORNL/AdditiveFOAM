---
layout: versioned
title: nLight AFX
parent: Tutorials
nav_order: 4
permalink: /tutorials/nlight-afx/
---

# nLight AFX Heat Source

This SS316L case represents an nLight AFX beam as a weighted combination of inner and outer Gaussian rings with projected axial distributions.

## Physical setup

- SS316L properties from `$ADDITIVEFOAM_ETC/materials/SS316L.cfg`.
- One `nLightAFX` heat source composed of inner and outer Gaussian rings.
- Seven characterized beam modes, `Index0` through `Index6`.
- A projected axial distribution with transient source-depth updating.

## Run

```bash
cp -r "$ADDITIVEFOAM_TUTORIALS/nLightAFX" "$FOAM_RUN/nLightAFX"
cd "$FOAM_RUN/nLightAFX"
./Allrun
```

## Important inputs

| File | Purpose |
|---|---|
| `constant/nLightAFX.cfg` | Defines the seven characterized beam modes |
| `constant/heatSourceDict` | Selects the mode, absorption model, source sampling, and AMR model |
| `constant/scanPath` | Beam position, time, and power |
| `constant/transportProperties` | Includes the SS316L material properties |
| `system/controlDict` | Run controls and optional Function Objects |

## Workflow

`constant/nLightAFX.cfg` defines `Index0` through `Index6`. The heat-source dictionary includes that file and expands one mode:

```foam
#include "nLightAFX.cfg"

nLightAFXCoeffs
{
    $Index6;
    transient true;
    nPoints   (10 10 10);
}
```

Change only `$Index6` to another defined index to select a different characterized profile. Each mode provides lateral dimensions, outer-ring fraction `alpha`, and `radius`, `sigma`, `A`, and `B` for both components.

Lower modes are more centre-weighted; higher modes place a larger share of power in the outer ring. The scan path sets total incident power. `alpha` partitions that power between the profile components and is distinct from the material absorptivity. The transient closure updates source depth using the local liquidus isotherm.

## Outputs

Open the reconstructed case in ParaView and visualize `T`, `alpha.solid`, and `qDot` to compare the selected inner- and outer-ring intensity distribution with the resulting temperature and phase fields.

Absorbed power, melt-pool dimensions, solidification data, and ExaCA histories are covered by the [baseline output workflows]({{ '/tutorials/amb2018/#outputs' | relative_url }}).
