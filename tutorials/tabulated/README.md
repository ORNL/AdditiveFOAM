# Tabulated Heat Source Tutorial

## Case description

This tutorial demonstrates the `tabulated` heat source model in AdditiveFOAM
using a measured nLight AFX Index 3 laser profile. The profile was exported from
PRIMES LaserDiagnosticsSoftware and converted to the AdditiveFOAM tabulated
heat-source format with `primesToAdditiveFoam`.

The purpose of this tutorial is to show how measured beam profiles can be used
directly in AdditiveFOAM without adding a new analytic heat-source model for
each laser shape.

This tutorial uses the AlSi10Mg material configuration from
`$ADDITIVEFOAM_ETC/materials/AlSi10Mg.cfg`.

The tabulated beam profile defines the measured two-dimensional laser intensity
distribution, while the scan path defines the applied laser power, scan speed,
and beam motion.

## File structure

The important files for this tutorial are:

`Allrun`

Converts the PRIMES beam-profile export to the AdditiveFOAM tabulated format,
then runs the mesh generation, decomposition, solver, and reconstruction steps.

`constant/heatSourceDict`

Defines the moving heat source and selects the `tabulated` heat-source model.

`constant/scanPath`

Defines the laser path, laser power, and scan speed. The tabulated profile only
defines the relative beam shape. The applied laser power still comes from the
scan path.

`constant/primes-export.csv`

PRIMES LaserDiagnosticsSoftware CSV export for the measured nLight AFX Index 3
laser profile.

`constant/beamProfile.txt`

AdditiveFOAM tabulated beam-profile file generated from
`constant/primes-export.csv`.

`constant/transportProperties`

Defines the thermophysical properties used for the AlSi10Mg simulation.

## Running the tutorial

Run the tutorial with:

```sh
./Allrun
```

The `Allrun` script first converts the PRIMES beam-profile export:

```sh
primesToAdditiveFoam constant/primes-export.csv constant/beamProfile.txt
```

It then runs the standard OpenFOAM workflow:

```sh
runApplication blockMesh
runApplication decomposePar
runParallel $application
runApplication reconstructPar
```

The converted file `constant/beamProfile.txt` is regenerated each time `Allrun`
is executed. This keeps the tutorial reproducible from the original PRIMES CSV
export.

## PRIMES to AdditiveFOAM conversion

The converter reads the measured PRIMES beam-profile export:

```text
constant/primes-export.csv
```

and writes the AdditiveFOAM tabulated heat-source profile:

```text
constant/beamProfile.txt
```

## Heat source model

The tutorial uses:

```foam
heatSourceModel tabulated;
```

The corresponding coefficient dictionary is:

```foam
tabulatedCoeffs
{
    file        "beamProfile.txt";
    dimensions  (2.50e-4 2.50e-4 5.0e-5);
    A           0;
    B           1;
    nPoints     (10 10 10);
}
```

### Coefficients

`file`

Name of the tabulated two-dimensional beam file. Relative paths are interpreted
relative to the `constant/` directory.

`dimensions`

Sets the heat-source dimensions used by the base moving heat-source integration.
For this tutorial, the lateral dimensions should cover the support of the
measured nLight AFX Index 3 beam profile. The third component is the initial
projected heat-source depth.

`transient`

Enables transient heat-source depth adjustment based on the local melt-pool
response.

`isoValue`

Optional temperature isovalue used by the transient projected heat-source
closure. If omitted, the material liquidus from `constant/transportProperties`
is used.

`A` and `B`

Define the projected axial shape closure:

```text
n = A*log2(x) + B
```

where `x` is the ratio between the current heat-source depth and lateral
heat-source size.

`nPoints`

Controls the sub-cell sampling resolution used when integrating the heat source
over mesh cells.

## Tabulated beam file format

The tabulated beam file is a headerless ASCII file. It must not contain comment
lines.

The required format is:

```text
nx ny
x0 y0
dx dy
f00 f10 f20 ... f(nx-1,0)
f01 f11 f21 ... f(nx-1,1)
...
f0(ny-1) ... f(nx-1,ny-1)
```

where:

| Entry | Meaning |
|---|---|
| `nx`, `ny` | Number of grid points in the x and y directions |
| `x0`, `y0` | Minimum x and y coordinates of the table, in meters |
| `dx`, `dy` | Uniform grid spacing in x and y, in meters |
| `f` | Relative planar intensity values |

The intensity values are stored in row-major order with `i` varying fastest:

```text
f_[i + nx*j]
```

The table is interpreted as nodal data. Bilinear interpolation is used between
nodes.

The valid interpolation domain is:

```text
x0 <= x <= x0 + (nx - 1)*dx
y0 <= y <= y0 + (ny - 1)*dy
```

Outside this domain, the heat source returns zero.

## Example table

The included tabulated profile:

```text
constant/beamProfile.txt
```

is generated from the PRIMES LaserDiagnosticsSoftware export:

```text
constant/primes-export.csv
```

The PRIMES file contains the measured nLight AFX Index 3 laser profile used in
this tutorial. The converted table is normalized so that the planar integral is
approximately 1.0. This normalization makes the tabulated profile a relative
beam-shape definition rather than an absolute-power input.

The scan-path laser power still controls the total applied power:

```text
constant/scanPath
```

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
