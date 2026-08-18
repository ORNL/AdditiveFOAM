# nLight AFX Heat Source Tutorial

## Case description

This tutorial demonstrates the `nLightAFX` heat source model in AdditiveFOAM.
The model represents an nLight AFX beam as a linear combination of two
Gaussian-ring components with the same projected axial distribution used by
the `projectedGaussian` heat source.

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
#include "$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg"

nLightAFXCoeffs
{
    $Index6;

    nSlope      0.0;
    nIntercept  1.0;

    transient   true;
    nPoints     (10 10 10);
}
```

The shared configuration sets the initial heat-source depth to 20 microns.

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
r0
sigma0
r1
sigma1
```

`alpha`

Integrated fraction of laser power assigned to component 1, which represents
the outer ring in the characterized AFX modes.

`r0`, `sigma0`, `r1`, and `sigma1`

Define the radial beam-shape parameters for Gaussian components 0 and 1. The
values are specified in meters.

`dimensions`

Sets the lateral heat source dimensions to half of D4sigma. The third component
is the initial projected depth and is updated by the transient heat source
logic when `transient` is enabled.

`nSlope` and `nIntercept`

Define the projected axial shape closure shared by both components:

```text
n = clip(nSlope*log2(x) + nIntercept, 0, 9)
k = 2^n
```

where `x` is the ratio between the current heat source depth and lateral heat
source size, and `k` is the exponent in the axial decay. This is the same
notation and closure used by the `projectedGaussian` model.

`transient`

Enables heat-source depth adjustment using the depth of the selected
temperature isotherm. The solver reports this calculated value as `isoDepth`.

`isoValue`

Optional temperature isovalue used to calculate `isoDepth`. If omitted, the
material liquidus from `constant/transportProperties` is used.

## Example mode

A typical mode block from `nLightAFX-1000.cfg` looks like:

```foam
Index3
{
    dimensions  (1.09290e-4 1.09290e-4 20.0e-6);

    alpha       0.483;

    r0          14.39e-6;
    sigma0      20.78e-6;
    r1          100.98e-6;
    sigma1      16.92e-6;
}
```

The selected mode is combined with the shared projection coefficients:

```foam
nLightAFXCoeffs
{
    $Index3;

    nSlope      0.0;
    nIntercept  1.0;
}
```

This mode has approximately half of the power in component 1. Lower modes are
more center-weighted, while higher modes place more power in the outer ring.

The component integrals are accounted for independently so `alpha` remains the
integrated power fraction in component 1. If `a0` and `a1` are the one-sided
volume integrals of the two Gaussian components with their shared axial
distribution, the implemented weight and normalization are

```text
gi(r) = exp(-0.5*((r - ri)/sigmai)^2)
      + exp(-0.5*((r + ri)/sigmai)^2)
s(z) = exp(-3*(-z/d)^k), z <= 0
w = s(z)*((1 - alpha)*g0(r) + alpha*g1(r)*a0/a1)
V0 = a0
```

The analytic component integral is

```text
ai = 2*pi*sigmai*d*Gamma(1/k)/(k*3^(1/k))
   * (2*sigmai*exp(-0.5*(ri/sigmai)^2)
      + sqrt(2*pi)*ri*erf(ri/(sqrt(2)*sigmai)))
```

Therefore, the normalized component contributions integrate to `1 - alpha`
and `alpha`, and the total normalized heat source integrates to one.

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
