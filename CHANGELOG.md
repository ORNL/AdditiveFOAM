# Unreleased

## Compatibility

- Ported AdditiveFOAM to OpenFOAM-14, including custom boundary conditions,
  moving-mesh/MRF integration, utilities, case headers, tutorials, and CI.

## Bug Fixes

- Corrected the legacy Marangoni `thermoPath` test data to use the pair syntax
  required by the `thermoPath` reader.
- Updated the multi-layer powder initialisation to the OpenFOAM-14 `setFields`
  zone syntax.

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
