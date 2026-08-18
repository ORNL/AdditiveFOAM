# Release 2.0.0

## Features

- Ported AdditiveFOAM to OpenFOAM-14.
- Added dynamic refinement models and reusable material and heat-source
  configuration files.
- Added `nLightAFX` and tabulated heat-source models, including conversion of
  PRIMES beam measurements with `primesToAdditiveFoam`.
- Added tabulated-profile inspection and automatic D4sigma calculation.
- Added an adaptive Bayesian workflow for calibrating the projected
  heat-source depth distribution.
- Expanded the AMB2018, multi-beam, multi-layer, nLight AFX, and tabulated
  tutorials.

## Bug fixes

- Made the projected Gaussian, nLight AFX, and tabulated heat sources
  one-sided, consistent with their analytic normalization.
- Preserved an explicitly configured Kelly `aspectRatioSwitch` value and
  applied `etaMin` as a lower bound on the Kelly absorption curve.

## Upgrade notes

- AdditiveFOAM 2.0 requires OpenFOAM-14. Source the OpenFOAM-14 environment
  before sourcing `etc/bashrc` and rebuild AdditiveFOAM with `./Allwmake`.
- The `A` and `B` projected-depth coefficients are now named `nSlope` and
  `nIntercept` in `projectedGaussian`, `tabulated`, and `nLightAFX`.
- The tabulated model now derives its lateral dimensions and D4sigma metrics
  from the beam-profile table. Replace the old `dimensions` entry with
  `minimumDepth`; no user-supplied lateral dimensions are required.
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
