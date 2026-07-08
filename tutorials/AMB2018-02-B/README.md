# AMB2018-02-B Tutorial

This tutorial demonstrates an AdditiveFOAM single-track case calibrated to the AMBenchmark 2018 AMB2018-02-B single-track data.

The purpose of this tutorial is to provide a calibrated AMBench baseline case using a Gaussian heat source model.

## Reference

This case is based on the AMBenchmark 2018 AMB2018-02 description:

```text
https://www.nist.gov/ambench/amb2018-02-description
```

The model coefficients (absorption and heat source dimensions) were calibrated in:
```text
    G.L. Knapp, J. Coleman, M. Rolchigo, M. Stoyanov, A. Plotkowski,
    Calibrating uncertain parameters in melt pool simulations of additive 
    manufacturing (2023), https://doi.org/10.1016/j.commatsci.2022.11190.
```

## File structure

The important files for this tutorial are:

```text
constant/heatSourceDict
constant/scanPath
constant/dynamicMeshDict
system/blockMeshDict
```

```text
constant/heatSourceDict
```

Defines the moving heat source, absorption model, and mesh refinement model.

```text
constant/scanPath
```

Defines the laser path, laser power, and scan speed or dwell time.

```text
constant/dynamicMeshDict
```

Defines the dynamic mesh/refinement settings.

```text
system/blockMeshDict
```

Defines the base computational mesh.

## Heat source model

The tutorial uses one heat source:

```foam
sources (beam);
```

The source uses constant absorptivity and a `superGaussian` heat source:

```foam
beam
{
    pathName            scanPath;

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

Constant absorptivity applied to the laser power from the scan path.

`k`

Super-Gaussian shape exponent. In this tutorial, `k = 2.0`, giving a Gaussian-like source.

`dimensions`

Sets the heat source dimensions used by the moving heat source integration, taken as `2sigma`

`nPoints`

Controls sub-cell sampling resolution used when integrating the heat source over mesh cells.

## Refinement model

This tutorial has the option to use the `targetCellLoad` refinement model:

```foam
refinementModel
{
    refinementModel         targetCellLoad;
    
    refine                  true;
    nLevels                 1;
    refinementTemperature   1000;

    buffer                  (85.0e-6 85.0e-6 100e-6);

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

The refinement region projects along the heat source path and targets a desired cell load per processor.
