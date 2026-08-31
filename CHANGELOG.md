# Release 2.0.0

## Features

- Ported AdditiveFOAM to OpenFOAM-14.
- Added dynamic refinement models and reusable material and heat-source
  configuration files.
- Added a general projected heat source with super-Gaussian, nLight AFX, and
  tabulated planar profiles and an exponential axial projection.
- Added conversion of PRIMES beam measurements with `primesToAdditiveFoam`.
- Added general planar-profile moment metrics, tabulated-profile inspection,
  and automatic `D4Sigma` calculation.
- Added rotated elliptical analytic profiles, `e2` and `secondMoment` radius
  definitions, and selectable area-equivalent, major, or minor `D4Sigma`
  reference widths.
- Added model-specific heat-source bounds that retain at least
  `1-tolerance` of the analytical source power.
- Added support-based heat-source quadrature.
- Added an adaptive Bayesian workflow with deterministic posterior quadrature
  for calibrating the projected heat-source axial closure.
- Expanded the AMB2018, multi-beam, multi-layer, nLight AFX, and tabulated
  tutorials.

## Bug fixes

- Made all heat sources one-sided, consistent with their analytic
  normalization.
- Synchronized heat-source depth, shape parameters, normalization, and bounds
  once before beam subcycling.
- Preserved the active-power fraction when a solver time step crosses the end
  of a scan path.
- Restricted isotherm source-depth measurements to the planar source support
  below the beam plane and removed a modified super-Gaussian endpoint
  singularity.
- Added calculation-ROI reconstruction, coordinate-preserving support cropping,
  and beam-statistics reporting to `primesToAdditiveFoam`.
- Preserved an explicitly configured Kelly `aspectRatioSwitch` value and
  applied `etaMin` as a lower bound on the Kelly absorption curve.

## Upgrade notes

- AdditiveFOAM 2.0 requires OpenFOAM-14. Source the OpenFOAM-14 environment
  before sourcing `etc/bashrc` and rebuild AdditiveFOAM with `./Allwmake`.
- Sources now use nested `absorption`, `heatSource`, `profile`, `projection`,
  and `refinement` dictionaries with a `model` entry instead of `*Coeffs`
  dictionaries.
- The former projected Gaussian, nLight AFX, and tabulated heat-source models
  are profiles selected by the `projected` model. Their shared axial model is
  selected in the `projection` dictionary.
- `superGaussian` remains a volumetric model and is also a planar profile.
  `radius` is two-dimensional and `depth` is a separate scalar.
- The tabulated profile derives its D4Sigma metrics from the beam-profile
  table. No user-supplied lateral dimensions are required.
- Profile metrics expose `D4Sigma` as `(major minor)`. Reference dimensions
  select the scalar width and constant or isotherm depth used by aspect-ratio
  closures.
- Aspect-ratio closures use a `D4Sigma` width reference, with selectable
  `areaEquivalent`, `major`, or `minor` component, and a `constant` or
  `isotherm` depth reference.
- The `A` and `B` projected-depth coefficients are named `nSlope` and
  `nIntercept` in the `projection` dictionary.
- Python 3.10 through 3.12 and the packages in `requirements.txt` are required
  only for the calibration and plotting utilities.

## Full Changelog

For a complete list of changes, see the [full changelog](https://github.com/ORNL/AdditiveFOAM/compare/1.2.0...2.0.0).


# Release 1.2.0

## Features

- Added the projected Gaussian heat-source model.
- Added melt-pool-dimension and solidification-data function objects.
- Added optional implicit temperature limiting through the PIMPLE dictionary.
- Added multi-layer scan-path support and updated the ExaCA tutorials.
- Added the top-level build script for solvers and utilities.

## Bug Fixes

- Corrected thermophysical file handling for cases launched with `-case`.

## Full Changelog

For a complete list of changes, see the [full changelog](https://github.com/ORNL/AdditiveFOAM/compare/1.1.0...1.2.0).


# Release 1.1

## Features
- **CI setup** [\#4](https://github.com/ORNL/AdditiveFOAM/pull/4)
- **Time-parallel event-based interpolation functionality for ExaCA** [\#14](https://github.com/ORNL/AdditiveFOAM/pull/14)
- **Adaptive heat source integration using Riemann sum** [\#15](https://github.com/ORNL/AdditiveFOAM/pull/15)
- **Utility for creating scan paths** [\#42](https://github.com/ORNL/AdditiveFOAM/pull/42)
- **Function object for ExaCA with improved tutorials** [\#39](https://github.com/ORNL/AdditiveFOAM/pull/39)

## Improvements
- **Added SPACK installation details in README** [\#13](https://github.com/ORNL/AdditiveFOAM/pull/13)
- **Version and build info for tracking** [\#50](https://github.com/ORNL/AdditiveFOAM/pull/50) [\#54](https://github.com/ORNL/AdditiveFOAM/pull/54)

## Bug Fixes
- **Corrected Marangoni patch evaluation to a slip condition and added tutorial** [\#11](https://github.com/ORNL/AdditiveFOAM/pull/11)
- **Resolved explicit energy solve to only use forward Euler method** [\#37](https://github.com/ORNL/AdditiveFOAM/pull/37)
- **Corrected the normalization factor V0 in the modifiedSuperGaussian heat source model** [\#20](https://github.com/ORNL/AdditiveFOAM/pull/20)

## Full Changelog
For a complete list of changes, see the [full changelog](https://github.com/ORNL/AdditiveFOAM/compare/1.0.0...1.1.0).


# Release 1.0

Initial Release
