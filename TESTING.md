# Testing

AdditiveFOAM has two test layers:

- Native C++ unit-style tests under [`tests/`](tests), built with `wmake`, linked against the existing AdditiveFOAM/OpenFOAM libraries, and driven by GoogleTest.
- The tutorial smoke run in GitHub Actions, which checks end-to-end integration.

## Prerequisites

- An OpenFOAM-13 environment must be sourced before building or running tests.
- AdditiveFOAM must be built first so the native tests can link against `libmovingBeamModels`.
- GoogleTest must be installed outside the repository. `./tests/Allwmake` auto-detects standard install locations and also honors `GTEST_ROOT`, `GTEST_INCLUDE_DIR`, and `GTEST_LIB_DIR`.

## Install GoogleTest

On Debian/Ubuntu systems:

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends cmake libgtest-dev

if [ ! -f /usr/lib/x86_64-linux-gnu/libgtest.a ] && [ ! -f /usr/local/lib/libgtest.a ]; then
  cmake -S /usr/src/googletest -B /tmp/googletest-build
  cmake --build /tmp/googletest-build -j2
  sudo cp /tmp/googletest-build/lib/libgtest*.a /usr/local/lib/
fi
```

If GoogleTest is installed somewhere else, set `GTEST_ROOT` or both
`GTEST_INCLUDE_DIR` and `GTEST_LIB_DIR` before running `./tests/Allwmake`.

## Build And Run

From the repository root:

```bash
. /path/to/openfoam/etc/bashrc
./Allwmake
./tests/Allwmake
./tests/run
```

If GoogleTest is installed in a non-standard location, export one of these before `./tests/Allwmake`:

```bash
export GTEST_ROOT=/path/to/gtest
# or
export GTEST_INCLUDE_DIR=/path/to/include
export GTEST_LIB_DIR=/path/to/lib
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
