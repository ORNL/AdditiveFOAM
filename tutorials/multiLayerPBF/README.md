# Multi-Layer PBF Tutorial

## Case description

This tutorial demonstrates a multi-layer powder bed fusion workflow using AdditiveFOAM.

The purpose of this tutorial is to show how to run repeated layer simulations with a powder layer using a transient volumetric heat source and depth-dependent absorption.

This tutorial uses the IN625 material configuration from
`$ADDITIVEFOAM_ETC/materials/IN625.cfg`.

## Running the tutorial

Build AdditiveFOAM against OpenFOAM-14, source both environments, and run:

```sh
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
cd "$ADDITIVEFOAM_TUTORIALS/multiLayerPBF"
./Allrun
```

Use `./Allclean` to remove generated layer and result files.

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

`eta0`

Fresnel absorption of the liquid metal.

`etaMin`

Effective absorption in the conduction regime, used when the melt-pool aspect
ratio is less than `aspectRatioSwitch`.

`aspectRatioSwitch`

Optional aspect-ratio cutoff for switching between `etaMin` and the Kelly
multiple-reflection model. If omitted, the cutoff is `1.0`.

`dimensions`

Sets the heat source dimensions. The third component is the initial projected depth.

`m` and `k`

Shape parameters for the `modifiedSuperGaussian` heat source.

`transient`

When `true`, AdditiveFOAM updates the heat source depth using the material
liquidus from `constant/transportProperties`.

`nPoints`

Controls sub-cell sampling resolution used when integrating the heat source over mesh cells.

## Post-processing

Optional AdditiveFOAM function objects are listed in `system/controlDict` and
are controlled by their `enabled` entries. To write additional data, set the
selected function object entry to:

```foam
enabled true;
```

To disable a function object, set:

```foam
enabled false;
```

`meltPoolDimensions` writes melt-pool length, width, and depth data.
`solidificationData` writes solidification events for CET analysis.
`ExaCA` writes temperature history data for ExaCA input files.

The `Allrun` script calls the reconstruction helpers after all layers finish:

```sh
reconstructExaCAData

reconstructSolidificationData
```

These commands detect `layer*/` directories and exit quietly when no matching
function object data were written.

After reconstructing the temperature data, run ExaCA from the base case
directory so that the layerwise relative paths in `ExaCA/input.json` resolve
correctly:

```sh
mpirun -np <nProcs> <path-to-ExaCA> ExaCA/input.json
```

Plot absorbed power from the layer solver logs with:

```sh
plotPower layer0 layer1
```

Plot melt-pool dimensions:

```sh
plotDimensions layer0 layer1
```

Plot CET data:

```sh
plotCET layer0 layer1
```
