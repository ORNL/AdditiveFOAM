# Multi-Layer PBF Tutorial

This tutorial demonstrates a multi-layer powder bed fusion workflow using AdditiveFOAM.

The purpose of this tutorial is to show how to run repeated layer simulations with a powder layer using a transient volumetric heat source and depth-dependent absorption.

This tutorial uses the IN625 material configuration from
`$ADDITIVEFOAM_ETC/materials/IN625.cfg`.

## Reference

The heat source setup follows the dynamic volumetric heat source formulation cited in the tutorial dictionary:

```text
J. Coleman, G.L. Knapp, B. Stump, M. Rolchigo, K. Kincaid, A. Plotkowski,
A dynamic volumetric heat source model for laser additive manufacturing,
Additive Manufacturing (2024), https://doi.org/10.1016/j.addma.2024.104531.
```

## File structure

The important files for this tutorial are:

```text
constant/heatSourceDict
constant/scanPath
system/extrudeMeshDict
```

```text
constant/heatSourceDict
```

Defines the moving heat source and absorption model.

```text
constant/scanPath
```

Defines the laser path, laser power, and scan speed or dwell time for the current layer.

```text
system/extrudeMeshDict
```

Defines the mesh extrusion behavior used by the multi-layer workflow.

## Heat source model

The tutorial uses one heat source:

```foam
sources (beam);
```

The source uses the `Kelly` absorption model and a transient `modifiedSuperGaussian` heat source:

```foam
beam
{
    pathName            scanPath;

    absorptionModel     Kelly;

    KellyCoeffs
    {
        geometry        cone;
        eta0            0.28;
        etaMin          0.35;
    }

    heatSourceModel     modifiedSuperGaussian;

    modifiedSuperGaussianCoeffs
    {
        dimensions      (40.0e-6 40.0e-6 30e-6);
        m               2.72;
        k               7.95;
        transient       true;
        nPoints         (10 10 10);
    }
}
```

### Coefficients

`geometry`

Absorption geometry used by the `Kelly` absorption model.

`eta0` and `etaMin`

Parameters for the depth-dependent absorption model.

`dimensions`

Sets the heat source dimensions. The third component is the initial projected depth.

`m` and `k`

Shape parameters for the `modifiedSuperGaussian` heat source.

`transient`

When `true`, AdditiveFOAM updates the heat source depth using the material
liquidus from `constant/transportProperties`.

`nPoints`

Controls sub-cell sampling resolution used when integrating the heat source over mesh cells.

## Notes

This tutorial is the main example for the dynamic volumetric heat source formulation in a multi-layer process setting. It is useful for studying layer-by-layer thermal history, melt pool depth evolution, and coupling to downstream solidification or microstructure workflows.
