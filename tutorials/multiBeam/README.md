# AMB2018-02-B multiBeam Tutorial

This tutorial demonstrates a simple two-beam extension of the AdditiveFOAM
AMB2018-02-B single-track tutorial.

The baseline AMB2018-02-B tutorial is calibrated to the AMBenchmark 2018
single-track data using a Gaussian heat source model. This tutorial keeps the
same calibrated beam parameters, scan speed, power, absorption, and heat source dimensions, but adds a second beam running parallel to the first beam.

The two beams are separated by a 100 micron center-to-center offset in the
hatch direction.

This tutorial is intended to show how neighboring beams can alter the thermal
field and solidification conditions relative to a single-beam setup.

This tutorial uses the IN625 material configuration from
`$ADDITIVEFOAM_ETC/materials/IN625.cfg`.

## Reference

The baseline case is based on the AMBenchmark 2018 AMB2018-02 description:

```text
https://www.nist.gov/ambench/amb2018-02-description
```

The model coefficients used by the AMB2018-02-B tutorial, including absorption
and heat source dimensions, were calibrated in:

```text
G.L. Knapp, J. Coleman, M. Rolchigo, M. Stoyanov, A. Plotkowski,
Calibrating uncertain parameters in melt pool simulations of additive
manufacturing (2023), https://doi.org/10.1016/j.commatsci.2022.11190.
```

## Case description

The original AMB2018-02-B tutorial contains one moving heat source. This tutorial
uses two moving heat sources:

```foam
sources (beam1 beam2);
```

The two beams follow the same scan direction, power, speed, and heat source
parameters of the original AMB2018-02-B scan path, but are translated `+/- 50`
microns in the hatch direction.

Because each beam uses the original AMB2018-02-B laser power, the total applied
laser power is twice that of the single-beam baseline.

## File structure

The important files for this tutorial are:

```text
constant/heatSourceDict
constant/scanPath_0
constant/scanPath_1
constant/dynamicMeshDict
system/blockMeshDict
```

```text
constant/heatSourceDict
```

Defines the two moving heat sources, their absorption models, heat source
models, and mesh refinement model.

```text
constant/scanPath_0
```

Defines the first beam path.

```text
constant/scanPath_1
```

Defines the second beam path.

```text
constant/dynamicMeshDict
```

Defines the dynamic mesh/refinement settings.

```text
system/blockMeshDict
```

Defines the base computational mesh.

## Heat source model

This tutorial uses two heat sources:

```foam
sources (beam1 beam2);
```

Both beams use the same heat source parameters as the calibrated AMB2018-02-B
single-beam case.

```foam
beam1
{
    pathName            scanPath_0;

    absorptionModel     constant;

    constantCoeffs
    {
        eta             0.33;
    }

    heatSourceModel     superGaussian;

    superGaussianCoeffs
    {
        k               2.0;
        dimensions      (85.0e-6 85.0e-6 30e-6);
        nPoints         (10 10 10);
    }
}

beam2
{
    pathName            scanPath_1;

    absorptionModel     constant;

    constantCoeffs
    {
        eta             0.33;
    }

    heatSourceModel     superGaussian;

    superGaussianCoeffs
    {
        k               2.0;
        dimensions      (85.0e-6 85.0e-6 30e-6);
        nPoints         (10 10 10);
    }
}
```

### Coefficients

`eta`

Constant absorptivity applied to the laser power from each scan path.

`k`

Super-Gaussian shape exponent. In this tutorial, `k = 2.0`, giving a
Gaussian-like source.

`dimensions`

Sets the heat source dimensions used by the moving heat source integration,
taken as `2sigma`.

`nPoints`

Controls the sub-cell sampling resolution used when integrating each heat source
over mesh cells.

## Refinement model

This tutorial can use the same `targetCellLoad` refinement model as the
single-beam AMB2018-02-B tutorial. The refinement buffer is applied relative to
each path interval for each beam.

A representative refinement setup is:

```foam
refinementModel
{
    refinementModel         targetCellLoad;

    refine                  true;
    nLevels                 1;
    refinementTemperature   1000;

    buffer                  (85.0e-6 200.0e-6 100e-6);

    minRefineVolumeFactor   4;

    volumeSearchMaxIter        10;
    volumeSearchTimeTolerance  1e-4;

    targetCellLoadCoeffs
    {
        cellsPerProc                   5000;
        targetVolumeSafetyFactor       1.0;
        maxTargetVolumeGrowth          1.25;
        maxTargetVolumeShrink          0.8;
        postScanUpdateIntervalFactor   10;
    }
}
```
