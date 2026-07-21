# Heat-source calibration tutorial

This tutorial calibrates the projected depth-distribution closure used by
AdditiveFOAM heat sources. The supplied worked example uses a
`projectedGaussian` source with SS316L. The laser D4sigma diameter was measured
to be 109.69 microns, giving a 2sigma heat-source radius of 54.845 microns.

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
│   └── P187p5_V500_D109p69/
│       ├── B0/
│       ├── B4p5/
│       └── B9/
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

## Rendered projected-source parameters

The template exposes the lateral source dimension once in
`constant/heatSourceDict`:

```foam
heatSourceRadius <<heatSourceRadius>>;
```

The projected Gaussian coefficients use that value for both lateral dimensions
and render the trial value directly:

```foam
projectedGaussianCoeffs
{
    dimensions ($heatSourceRadius $heatSourceRadius 30.0e-6);
    A 0.0;
    B <<B>>;
}
```

For each experiment, `Spot_Size_microns` is the measured D4sigma diameter. The
calibration command divides it by two to obtain the heat-source radius, then
uses that radius for both lateral source dimensions and to normalize measured
and simulated depth. The calibration coordinate is therefore
`d / (D4sigma / 2)`, matching `d / min(dimensions.x(), dimensions.y())` in the
heat-source models. Power, speed, end time, and write interval are also
rendered for each case.
