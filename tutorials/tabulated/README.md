# Tabulated Heat Source Tutorial

## Case description

This tutorial demonstrates a `projectedHeatSource` with a measured `tabulated`
planar profile and an `exponential` axial projection. The nLight AFX Index 3
profile was exported from PRIMES LaserDiagnosticsSoftware and converted to the
AdditiveFOAM tabulated format with `primesToAdditiveFoam`.

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

Defines the moving heat source and selects its profile and projection.

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

Build AdditiveFOAM against OpenFOAM-14, source both environments, and run:

```sh
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
cd "$ADDITIVEFOAM_TUTORIALS/tabulated"
./Allrun
```

Use `./Allclean` to remove generated mesh, decomposition, and result files.

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

The supported input is an LDS table export with the calculation Region of
Interest (ROI) enabled at a fill factor of 0.5. A rotated-moments export provides
the `Radius a`, `Radius b`, and `Azimuth angle φ` values used for the
beam-statistics comparison.

`primesToAdditiveFoam` subtracts `Nullvalue`, sets samples below the resulting
zero level to zero, and updates the rectangular calculation ROI from the pixel
centroid and x/y second moments until the pixel mask converges. `SNR` is
diagnostic metadata and is not an intensity threshold.

The normalized output retains one zero-valued point outside each cropped edge.
This boundary preserves the bilinear profile, its coordinate system, centroid,
asymmetry, and beam statistics. The conversion report compares the profile
radii and azimuth with the corresponding PRIMES header values.

## Heat source model

The tutorial uses:

```foam
heatSourceModel projectedHeatSource;
```

The corresponding coefficient dictionary is:

```foam
projectedHeatSourceCoeffs
{
    profile         tabulated;
    projection      exponential;

    minimumDepth    5.0e-5;

    tabulatedCoeffs
    {
        file        "beamProfile.txt";
    }

    exponentialCoeffs
    {
        nSlope      0;
        nIntercept  1;
    }

    tolerance       1.0e-3;
    nPoints         (10 10 10);
}
```

### Coefficients

`file`

Name of the tabulated two-dimensional beam file. Relative paths are interpreted
relative to the `constant/` directory.

`minimumDepth`

Sets the minimum projected heat-source depth. The profile file determines the
lateral characteristic size, interpolation bounds, normalization, and D4Sigma
metrics.

The profile integral and D4Sigma metrics use the exact raw moments of the
bilinear interpolant, `Mpq = integral(x^p*y^q*I(x,y) dx dy)`.

`transient`

Enables transient heat-source depth adjustment based on the local melt-pool
response.

`isoValue`

Optional temperature isovalue used by the transient projected heat-source
closure. If omitted, the material liquidus from `constant/transportProperties`
is used.

`nSlope` and `nIntercept`

Define the projected axial shape closure:

```text
n = clip(nSlope*log2(a) + nIntercept, 0, 9)
k = 2^n
```

where `a = 2*depth/D4Sigma`. Here `D4Sigma` is the area-equivalent diameter,
`sqrt(D4SigmaMajor*D4SigmaMinor)`. The source is applied below the beam
plane with axial weight `exp(-3*(-z/depth)^k)`; its normalization uses the
matching one-sided axial integral.

`tolerance`

Optional maximum fraction of analytical source power outside the integration
bounds. The default `1e-3` retains at least 99.9% of the source power.

`nPoints`

Sets the target sub-cell spacing by dividing the retained source bounds by
`nPoints`.

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
