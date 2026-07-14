# AMB2018-02-B Tutorial

## Case description

This tutorial demonstrates an AdditiveFOAM single-track case calibrated to the AMBenchmark 2018 AMB2018-02-B single-track data.

The purpose of this tutorial is to provide a calibrated AMBench baseline case using a Gaussian heat source model.

This tutorial uses the IN625 material configuration from
`$ADDITIVEFOAM_ETC/materials/IN625.cfg`.

## Running the tutorial

Build AdditiveFOAM against OpenFOAM-14, source both environments, and run:

```sh
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
cd "$ADDITIVEFOAM_TUTORIALS/AMB2018-02-B"
./Allrun
```

Use `./Allclean` to remove generated mesh, decomposition, and result files.

## Reference

This case is based on the AMBenchmark 2018 AMB2018-02 description:

```text
https://www.nist.gov/ambench/amb2018-02-description
```

The model coefficients (absorption and heat source dimensions) were calibrated in:
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
    refinementModel         none;

    //refinementModel         targetCellLoad;

    refinementTemperature   1000;

    buffers
    {
        beam                (85.0e-6 85.0e-6 100e-6);
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

The refinement region projects along the heat source path and targets a desired cell load per processor.

### Coefficients

`refinementModel`

Selects the refinement model. `targetCellLoad` projects refinement ahead along
the scan path and adjusts the projected volume to target a cell count per
processor.

`refinementTemperature`

Temperature threshold used to mark hot cells for refinement during AMR updates.

`buffers`

Source-specific scan-path projection buffers. Entries are keyed by heat source
name. The vector components define the buffer size in the scan direction,
transverse direction, and build direction.

`targetCellsPerProc`

Target cell count per MPI processor. The `targetCellLoad` model adjusts the
projected refinement volume to keep the total mesh size near this load.

`nBufferVolumes`

Minimum projected scan-path volume expressed as a multiple of the source-buffer
volume. This keeps the projected refinement region from shrinking below the
local scan-path coverage.

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
