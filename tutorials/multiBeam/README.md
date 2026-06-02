# Multi-Beam Tutorial

This tutorial demonstrates how to use multiple moving heat sources in a single AdditiveFOAM case.

The purpose of this tutorial is to show that separate beams can use different scan paths, absorption models, and heat source models in the same simulation.

## References

The first beam uses the calibrated fixed heat source formulation:

```text
G.L. Knapp, J. Coleman, M. Rolchigo, M. Stoyanov, A. Plotkowski,
Calibrating uncertain parameters in melt pool simulations of additive
manufacturing (2023), https://doi.org/10.1016/j.commatsci.2022.11190.
```

The second beam uses the dynamic volumetric heat source formulation:

```text
J. Coleman, G.L. Knapp, B. Stump, M. Rolchigo, K. Kincaid, A. Plotkowski,
A dynamic volumetric heat source model for laser additive manufacturing,
Additive Manufacturing (2024), https://doi.org/10.1016/j.addma.2024.104531.
```

## File structure

The important files for this tutorial are:

```text
constant/heatSourceDict
constant/scanPath_1
constant/scanPath_2
```

```text
constant/heatSourceDict
```

Defines both moving heat sources.

```text
constant/scanPath_1
```

Defines the path, power, and scan speed or dwell time for the first beam.

```text
constant/scanPath_2
```

Defines the path, power, and scan speed or dwell time for the second beam.

## Heat source models

The tutorial defines two sources:

```foam
sources (beam1 beam2);
```

### Beam 1

`beam1` uses constant absorptivity and a fixed `superGaussian` heat source:

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
        dimensions      (85.0e-6 85.0e-6 30e-6);
        k               2.0;
        transient       false;
        nPoints         (10 10 10);
    }
}
```

This beam is a conventional Gaussian source with constant absorption.

### Beam 2

`beam2` uses the `Kelly` absorption model and a transient `modifiedSuperGaussian` heat source:

```foam
beam2
{
    pathName            scanPath_2;

    absorptionModel     Kelly;

    KellyCoeffs
    {
        geometry        cone;
        eta0            0.28;
        etaMin          0.35;
    }

    heatSourceModel     modifiedSuperGaussian;

    modifiedSuperGaussianCoeffs
    {
        dimensions      (40.0e-6 40.0e-6 30e-6);
        m               2.72;
        k               7.95;
        transient       true;
        isoValue        1620;
        nPoints         (10 10 10);
    }
}
```

This beam demonstrates a dynamic volumetric heat source with depth-dependent absorption.

## Notes

Each beam has its own `pathName`, so each source can follow a different scan path. This tutorial is useful for testing multi-laser workflows.

