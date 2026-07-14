# AMB2018-02-B multiBeam Tutorial

## Case description

This tutorial demonstrates a simple two-beam extension of the AdditiveFOAM
AMB2018-02-B single-track tutorial.

This tutorial is intended to show how neighboring beams can alter the thermal
field and solidification conditions relative to a single-beam setup.

This tutorial uses the IN625 material configuration from
`$ADDITIVEFOAM_ETC/materials/IN625.cfg`.

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

## Running the tutorial

Build AdditiveFOAM against OpenFOAM-14, source both environments, and run:

```sh
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
cd "$ADDITIVEFOAM_TUTORIALS/multiBeam"
./Allrun
```

Use `./Allclean` to remove generated mesh, decomposition, and result files.

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
manufacturing (2023), https://doi.org/10.1016/j.commatsci.2022.111904.
```

## File structure

The important files for this tutorial are:

```text
constant/heatSourceDict
```

Defines the two moving heat sources, their absorption models, heat source
models, and mesh refinement model.

```text
constant/scanPath_1
```

Defines the first beam path.

```text
constant/scanPath_2
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

beam2
{
    pathName            scanPath_2;

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
single-beam AMB2018-02-B tutorial. The `buffers` entries are keyed by source
name and are applied relative to each path interval for each beam.

A representative refinement setup is:

```foam
refinementModel
{
    refinementModel         none;

    //refinementModel         targetCellLoad;

    refinementTemperature   1000;

    buffers
    {
        beam1               (85.0e-6 85.0e-6 100e-6);
        beam2               (85.0e-6 85.0e-6 100e-6);
    }

    targetCellLoadCoeffs
    {
        targetCellsPerProc  5000;
        nBufferVolumes      4;
        maxSearchIter       10;
        timeTolerance       1e-4;
    }
}
```

### Coefficients

`refinementModel`

Selects the refinement model. `targetCellLoad` projects refinement ahead along
the scan paths and adjusts the projected volume to target a cell count per
processor.

`refinementTemperature`

Temperature threshold used to mark hot cells for refinement during AMR updates.

`buffers`

Source-specific scan-path projection buffers. Entries are keyed by heat source
name and are applied to the corresponding beam path. The vector components
define the buffer size in the scan direction, transverse direction, and build
direction.

`targetCellsPerProc`

Target cell count per MPI processor. The `targetCellLoad` model adjusts the
projected refinement volume to keep the total mesh size near this load.

`nBufferVolumes`

Minimum projected scan-path volume expressed as a multiple of the combined
source-buffer volume. This keeps the projected refinement region from shrinking
below the local scan-path coverage.

`maxSearchIter`

Maximum number of bisection iterations used when searching for the scan-path
time interval that gives the target refinement volume.

`timeTolerance`

Stopping tolerance, in seconds, for the scan-path time-interval search.

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
