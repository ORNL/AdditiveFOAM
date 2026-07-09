# nLight AFX Heat Source Tutorial

This tutorial demonstrates the `nLightAFX` heat source model in AdditiveFOAM. The model represents an nLight AFX beam as a linear combination of inner and outer Gaussian-ring components with the same projected axial distribution used by the `projectedGaussian` heat source.

The purpose of this tutorial is to show how ORNL-characterized nLight AFX beam profiles can be selected from dictionary inputs without hard-coding mode values in the heat source model.

This tutorial uses the IN625 material configuration from
`$ADDITIVEFOAM_ETC/materials/IN625.cfg`.

## File structure

The important files for this tutorial are:

```text
constant/nLightAFX.cfg
constant/heatSourceDict
constant/scanPath
```

```text
constant/nLightAFX.cfg
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
#include "nLightAFX.cfg"

nLightAFXCoeffs
{
    $Index6

    transient   true;
    nPoints     (10 10 10);
}
```

The selected mode can be changed by replacing:

```foam
$Index3
```

with one of:

```foam
$Index0
$Index1
$Index2
$Index3
$Index4
$Index5
$Index6
```

## Characterized AFX modes

The included `nLightAFX.cfg` file contains the ORNL-characterized AFX beam parameters for modes 0 through 6. Each mode defines:

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

A typical mode block from `nLightAFX.cfg` looks like:

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

## Power normalization

The nLight AFX mode parameters define only the relative beam shape. The laser power should still be set in `constant/scanPath`, exactly as with the other heat source models.

The heat source model normalizes the inner and outer ring components so that the total volumetric source integrates to the absorbed scan-path power.

## Notes

The model uses one shared transient heat source depth from the existing AdditiveFOAM heat source logic. This keeps the implementation consistent with the other projected heat source models while allowing the radial AFX beam shape to vary by mode.

The characterized mode values are provided in `nLightAFX.cfg` for convenience. Users can copy a mode block and modify `alpha`, `radius`, `sigma`, `A`, or `B` to represent their own measured beam profiles.
