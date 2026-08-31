# nLight AFX Heat Source Tutorial

## Case description

This tutorial demonstrates a `projected` heat source with an `nLightAFX` planar
profile and an `exponential` axial projection. The profile represents an
nLight AFX beam as a linear combination of two Gaussian-ring components.

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

The heat-source dictionary is:

```foam
#include "$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg"

widthReference D4Sigma;
depthReference isotherm;

heatSource
{
    model       projected;
    depth       20.0e-6;
    nPoints     (10 10 10);

    profile
    {
        model   nLightAFX;
        $Index6;
    }

    projection
    {
        model       exponential;
        nSlope      0.0;
        nIntercept  1.0;
    }
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

The profile calculates its radial D4Sigma from these coefficients. The
reported D4Sigma values in the shared configuration provide reference values
for the characterized modes.

`nSlope` and `nIntercept`

Define the projected axial shape closure shared by both components:

```text
n = clip(nSlope*log2(a) + nIntercept, 0, 9)
k = 2^n
```

where `a` is twice the current heat-source depth divided by the selected
D4Sigma reference width, and `k` is the exponent in the axial decay.

`tolerance`

Optional maximum fraction of analytical source power outside the integration
bounds. The default `1e-3` retains at least 99.9% of the source power.

`nPoints`

Sets the target sub-cell spacing by dividing the retained source bounds by
`nPoints`.

`widthReference`

`D4Sigma` selects the beam-plane reference width. `D4Sigma` optionally selects
`areaEquivalent` (default), `major`, or `minor` from the profile metrics.

`depthReference`

Selects `constant` (default) or `isotherm`. The solver reports the calculated
reference depth when `isotherm` is selected.

`isotherm`

Optional temperature used by the isotherm depth reference. If omitted, the
material liquidus from `constant/transportProperties` is used.

## Example mode

A typical mode block from `nLightAFX-1000.cfg` looks like:

```foam
Index3
{
    alpha       0.483;

    r0          14.39e-6;
    sigma0      20.78e-6;
    r1          100.98e-6;
    sigma1      16.92e-6;
}
```

The selected profile is combined with the projection coefficients:

```foam
projection
{
    model       exponential;
    nSlope      0.0;
    nIntercept  1.0;
}
```

This mode has approximately half of the power in component 1. Lower modes are
more center-weighted, while higher modes place more power in the outer ring.

The component integrals are accounted for independently so `alpha` remains the
integrated power fraction in component 1. If `J0` and `J1` are the planar
component integrals without the common factor of `2*pi`, the implemented
weight and normalization are

```text
Ii(r) = exp(-0.5*((r - ri)/sigmai)^2)
      + exp(-0.5*((r + ri)/sigmai)^2)
p(z) = exp(-3*(-z/d)^k), z <= 0
w = p(z)*((1 - alpha)*I0(r) + alpha*I1(r)*J0/J1)
V0 = 2*pi*J0*d*Gamma(1/k)/(k*3^(1/k))
```

where

```text
Ji = 2*sigmai^2*exp(-0.5*(ri/sigmai)^2)
   + sqrt(2*pi)*ri*sigmai*erf(ri/(sqrt(2)*sigmai))
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
