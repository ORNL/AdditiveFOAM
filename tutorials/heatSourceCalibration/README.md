# Heat-source calibration tutorial

This tutorial calibrates the projected depth-distribution closure used by
AdditiveFOAM heat sources. The supplied worked example uses a `tabulated`
source with SS316L and a circular Gaussian planar profile. The measured
D4sigma diameter is 109.69 microns. The finite tabulated profile is constructed
so its bilinear interpolant has that D4sigma exactly.

## Installation and setup

Build AdditiveFOAM against OpenFOAM 14, then source both environments:

```bash
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
```

Copy the tutorial before running it:

```bash
mkdir -p "$FOAM_RUN/AdditiveFOAM"
cp -r "$ADDITIVEFOAM_TUTORIALS/heatSourceCalibration" \
    "$FOAM_RUN/AdditiveFOAM/heatSourceCalibration"
cd "$FOAM_RUN/AdditiveFOAM/heatSourceCalibration"
```

Install the Python dependencies into the environment you intend to use before
running the calibration:

```bash
python -m pip install \
    -r "$ADDITIVEFOAM_PROJECT_DIR/requirements.txt"
python -m pip check
calibrateHeatSource --help
```

Do not run a calibration inside `$ADDITIVEFOAM_TUTORIALS`; the campaign writes
cases and reports beneath the tutorial directory.

## Run the calibration

Review `system/decomposeParDict` in the template before starting. The supplied
configuration uses 8 MPI ranks for each AdditiveFOAM case, evaluates ten trial
values for each of five experiments, and uses 2,000 posterior draws.

```bash
calibrateHeatSource --config config.yml
```

The configuration resolves relative paths from the location of `config.yml`.
Environment variables and `~` are supported in the three `paths` entries.

Generated output is contained in `campaign/`:

```text
campaign/
├── cases/
│   └── P187p5_V500_measured_beam/
│       ├── nIntercept0/
│       ├── nIntercept4/
│       └── nIntercept9/
├── simulations.yml
├── calibration_state.yml
├── calibration_fit.yml
└── reports/
    ├── calibration_report.pdf
    └── calibration_summary.csv
```

Successful cases retain their rendered inputs, solver log, post-processing
output, and `metrics.yml`. With `keep_successful: false`, processor and numeric
time directories are removed after their melt-pool dimensions are recorded.

## Material-derived liquidus

The case includes:

```foam
#include "$ADDITIVEFOAM_ETC/materials/SS316L.cfg"
```

The configuration selects `melt_pool_isovalue: liquidus`. The calibration
command expands `constant/transportProperties` with `foamDictionary`, obtains
the temperature paired with `alpha.solid = 0`, and reads the matching file from
`postProcessing/meltPoolDimensions`. The SS316L liquidus is therefore not
duplicated in `config.yml`, `heatSourceDict`, or `controlDict`.

## Beam profiles and rendered parameters

Beam profiles are named once in `config.yml`:

```yaml
profiles:
  measured_beam:
    file: template/constant/beam_profile.txt
```

Each experimental condition selects a profile by name:

```yaml
parameters:
  Power_W: 187.5
  Speed_mm_s: 500
  Profile: measured_beam
```

The tabulated coefficients require the profile, a minimum projected depth, and
the two closure coefficients. No lateral `dimensions` input is required:

```foam
tabulatedCoeffs
{
    file            "beam_profile.txt";
    minimumDepth    20.0e-6;
    nSlope          0.0;
    nIntercept      <<nIntercept>>;
}
```

The closure uses

```text
n = clip(nSlope*log2(x) + nIntercept, 0, 9)
k = 2^n
```

where `x = 2*depth/D4sigma` and `k` is the exponent in the axial
decay. Calibration trials set `nSlope` to zero, so each trial `nIntercept` is
also its applied `n` value.

The headerless `constant/beam_profile.txt` table contains a normalized circular
Gaussian:

```text
sampled Gaussian sigma = 27.8150918491 microns
I(x,y) proportional to exp(-(x^2 + y^2)/(2 sigma^2))
```

The 67-by-67 grid uses 2.5-micron spacing and spans -82.5 through +82.5
microns in both directions. Its bilinear planar integral is normalized to one.
The sampled Gaussian width differs slightly from `109.69/4` because D4sigma is
calculated from the finite table and its bilinear interpolation, rather than an
infinite analytic Gaussian.

At startup, the calibration uses `tabulatedProfileInfo` to inspect each
configured profile and report its integral, centroid, principal D4sigma
values, area-equivalent D4sigma, and azimuth. The reported `D4sigma` is
`sqrt(D4sigmaMajor*D4sigmaMinor)`, the diameter of the circular beam with the
same second-moment area. The calibration uses half this value to normalize
measured and simulated depths. The profile file and measured melt-pool depths
are therefore the only beam-size and calibration measurements supplied by the
user; there is no separate expected-D4sigma or tolerance input.

## Resuming a campaign

`simulations.yml` caches completed AdditiveFOAM runs. A cached result is reused
only when the beam profile, template case, AdditiveFOAM executable and
libraries, shared AdditiveFOAM configuration, OpenFOAM version, run command,
and melt-pool isovalue are unchanged. The local posterior state additionally
tracks the experimental depths and calibration settings. Changing any of these
inputs marks the affected results as stale. Stale simulation entries remain in
`simulations.yml` for provenance but are ignored; stale posterior entries are
replaced after recalibration.

The calibration consumes melt-pool depth only. Power, speed, end time, and
write interval are rendered for each case. Width is left available for a later
workflow extension and length is not included in the calibration output.
