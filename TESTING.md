# Testing

AdditiveFOAM has two test layers:

- Native C++ unit-style tests under [`tests/`](tests), built with `wmake` and linked against the existing AdditiveFOAM/OpenFOAM libraries.
- The tutorial smoke run in GitHub Actions, which checks end-to-end integration.

## Prerequisites

- An OpenFOAM-13 environment must be sourced before building or running tests.
- AdditiveFOAM must be built first so the native tests can link against `libmovingBeamModels`.

## Build And Run

From the repository root:

```bash
. /path/to/openfoam/etc/bashrc
./Allwmake
./tests/Allwmake
./tests/run
```

`./tests/Allwmake` builds the native test executables without changing the default production build path. `./tests/run` executes the complete native suite.

## Current Coverage

The native suite currently builds four executables:

- `additiveFoamSegmentTests` validates `Foam::segment` default construction and parsing.
- `additiveFoamMovingBeamTests` covers scan-path timing, index selection, interpolation, and timestep adjustment in `Foam::movingBeam`.
- `additiveFoamMovingHeatSourceModelTests` exercises absorption-model and heat-source-model math for the current beam model implementations.
- `additiveFoamUtilityTests` protects `interpolateXY` and the graph utilities used by solver setup and post-processing.

The `movingBeam` and heat-source-model tests use a small file-backed fixture case under [`tests/fixtures/movingHeatSourceCase`](tests/fixtures/movingHeatSourceCase) so the constructors read real OpenFOAM dictionaries and scan-path files.

## Adding A New Native Test

1. Create a new subdirectory under `tests/`.
2. Add a `Make/files` that includes `../shared/testMain.C`, your test source, and an `EXE` target name.
3. Add a `Make/options` file with the required include paths and linked AdditiveFOAM/OpenFOAM libraries.
4. Add the new directory to [`tests/Allwmake`](tests/Allwmake).
5. Add the produced executable to [`tests/run`](tests/run).

The vendored header at [`tests/vendor/doctest/doctest.h`](tests/vendor/doctest/doctest.h) keeps the harness self-contained and avoids extra package dependencies.
