# nLight AFX Heat Source Tutorial

## Case description

This tutorial demonstrates the `nLightAFX` heat source model in AdditiveFOAM. The model represents an nLight AFX beam as a linear combination of inner and outer Gaussian-ring components with the same projected axial distribution used by the `projectedGaussian` heat source.

The purpose of this tutorial is to show how ORNL-characterized nLight AFX beam profiles can be selected from dictionary inputs.

This tutorial uses the SS316L material configuration from
`$ADDITIVEFOAM_ETC/materials/SS316L.cfg`.

## Running the tutorial

Build AdditiveFOAM against OpenFOAM-14, source both environments, and run:

```sh
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
cd "$ADDITIVEFOAM_TUTORIALS/nLightAFX"
./Allrun
```

Use `./Allclean` to remove generated mesh, decomposition, and result files.

## File structure

The important files for this tutorial are:

```text
$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg
```

Defines the ORNL-characterized AFX mode parameters.

```text
constant/heatSourceDict
```

Defines the moving heat source model, absorption model, projected axial closure, and selected nLight AFX mode.

```text
constant/scanPath
```

Defines the laser path, laser power, and scan speed or dwell time.

## Heat source model

The tutorial uses:

```foam
heatSourceModel nLightAFX;
```

The corresponding coefficient dictionary is:

```foam
depth   5.0e-5;
innerA  0.0;
innerB  1.0;
outerA  0.0;
outerB  1.0;

#include "$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg"

nLightAFXCoeffs
{
    $Index6;

    transient   true;
    nPoints     (10 10 10);
}
```

The selected mode can be changed by replacing:

```foam
$Index6;
```

with one of:

```foam
$Index0;
$Index1;
$Index2;
$Index3;
$Index4;
$Index5;
$Index6;
```

## Characterized AFX modes

The shared `nLightAFX-1000.cfg` file contains the ORNL-characterized AFX beam
parameters for modes 0 through 6. Each mode defines:

```foam
dimensions
alpha

inner
{
    radius
    sigma
    A
    B
}

outer
{
    radius
    sigma
    A
    B
}
```

`alpha`

Fraction of the laser power applied to the outer Gaussian-ring component.

`inner` and `outer`

Define the radial beam-shape parameters for the inner and outer components. The `radius` and `sigma` values are specified in meters.

`dimensions`

Sets the lateral heat source dimensions from the characterized beam size. The third component is the initial projected depth and is updated by the transient heat source logic when `transient` is enabled.

`A` and `B`

Define the projected axial shape closure for each component:

```text
n = A*log2(x) + B
```

where `x` is the ratio between the current heat source depth and lateral heat source size. The implementation clamps this internal numerical exponent consistently with the `projectedGaussian` model.

## Example mode

A typical mode block from `nLightAFX-1000.cfg` looks like:

```foam
Index3
{
    dimensions  (1.09290e-4 1.09290e-4 5.0e-5);

    alpha       0.483;

    inner
    {
        radius  14.39e-6;
        sigma   20.78e-6;
        A       0.0;
        B       1.0;
    }

    outer
    {
        radius  100.98e-6;
        sigma   16.92e-6;
        A       0.0;
        B       1.0;
    }
}
```

This mode has approximately half of the power in the outer ring. Lower modes are more center-weighted, while higher modes place more power in the outer ring.

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

The `Allrun` script calls the reconstruction helpers after the solver finishes:

```sh
reconstructExaCAData

reconstructSolidificationData
```

These commands exit quietly when no matching function object data were written.

After reconstructing the temperature data, run ExaCA from the case directory so
that the relative paths in `ExaCA/input.json` resolve correctly:

```sh
mpirun -np <nProcs> <path-to-ExaCA> ExaCA/input.json
```

Plot absorbed power from the solver log with:

```sh
plotPower
```

Plot melt-pool dimensions:

```sh
plotDimensions
```

Plot CET data:

```sh
plotCET
```
