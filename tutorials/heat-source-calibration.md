---
layout: versioned
title: Heat-source Calibration
parent: Tutorials
nav_order: 6
permalink: /tutorials/heat-source-calibration/
usemathjax: true
---

# Heat-source Calibration

This tutorial calibrates the projected axial distribution of an SS316L `projectedGaussian` source against melt-pool depths measured from five SS316L single tracks on a bare plate. It demonstrates the complete `calibrateHeatSource` workflow: render and run trial cases, extract maximum liquidus melt-pool depth, infer a local shape parameter for each condition, fit the global $$A,B$$ closure, and generate uncertainty and report products.

Read [Projected Heat-source Calibration]({{ '/docs/heat-source-calibration/' | relative_url }}) for the derivation, likelihood, robust global fit, bootstrap interval, cache behavior, and transfer rules used by the command.

## Physical setup

The calibration data are SS316L single tracks produced on a bare plate. `experiments.yml` contains two measured melt-pool depths at each of five powers; scan speed and measured D4σ spot diameter are constant:

| Power (W) | Speed (mm/s) | D4σ spot diameter (µm) | Measured depth (µm) |
|---:|---:|---:|---:|
| 187.5 | 500 | 109.69 | 87.05, 98.10 |
| 300.0 | 500 | 109.69 | 168.57, 190.67 |
| 412.5 | 500 | 109.69 | 272.19, 264.07 |
| 525.0 | 500 | 109.69 | 359.24, 353.71 |
| 637.5 | 500 | 109.69 | 464.25, 467.01 |

The AdditiveFOAM setup contains:

- SS316L properties are included from `$ADDITIVEFOAM_ETC/materials/SS316L.cfg`.
- The source is a transient `projectedGaussian`. For each experiment, the driver converts the measured D4σ diameter to a 2σ radius and assigns that value to both lateral dimensions. The supplied 109.69 µm D4σ diameter therefore produces `(54.845 54.845)` µm lateral dimensions. Initial depth is 30 µm and `nPoints` is `(10 10 10)`.
- The cylindrical Kelly model uses `eta0 0.27`, `etaMin 0.27`, and an effectively zero aspect-ratio switch.
- A 2 mm linear track is rendered at the experimental power and 0.5 m/s speed.
- `targetCellLoad` refinement follows the beam, and `dynamicMeshDict` permits one refinement level.
- Fluid flow is disabled with `nOuterCorrectors 0`.
- `meltPoolDimensions` records the material-derived SS316L liquidus isotherm.
- `decomposeParDict` uses eight subdomains with `scotch` decomposition.

## Run

Copy the tutorial to a writable run directory, then start or resume the campaign:

```bash
mkdir -p "$FOAM_RUN/AdditiveFOAM"
cp -r "$ADDITIVEFOAM_TUTORIALS/heatSourceCalibration" \
    "$FOAM_RUN/AdditiveFOAM/heatSourceCalibration"
cd "$FOAM_RUN/AdditiveFOAM/heatSourceCalibration"
calibrateHeatSource --config config.yml
```

Do not run the campaign inside `$ADDITIVEFOAM_TUTORIALS`; it writes cases, caches, and reports beneath the copied tutorial. Each AdditiveFOAM trial uses eight MPI ranks.

## Important inputs

| File | Purpose |
|---|---|
| `experiments.yml` | Process conditions and repeated measured melt-pool depths |
| `config.yml` | Case paths, trial grid, posterior sampling, refinement tolerances, and reporting |
| `template/constant/heatSourceDict` | Projected source, absorption, dimensions, sampling, and calibration placeholders |
| `template/constant/scanPath` | Rendered beam power, track, and velocity |
| `template/constant/transportProperties` | SS316L material configuration |
| `template/system/controlDict` | Run time, write interval, and `meltPoolDimensions` output |
| `template/system/decomposeParDict` | Parallel decomposition for each trial |

## Workflow

The calibration driver treats every distinct `parameters` mapping in `experiments.yml` as an independent process condition. `Power_W`, `Speed_mm_s`, and `Spot_Size_microns` are required by the renderer and normalization logic.

### Rendered trial cases

The heat-source radius and trial value appear once in `template/constant/heatSourceDict`:

```foam
heatSourceRadius <<heatSourceRadius>>;

projectedGaussianCoeffs
{
    dimensions ($heatSourceRadius $heatSourceRadius 30.0e-6);
    A          0.0;
    B          <<B>>;
    transient  true;
    nPoints    (10 10 10);
}
```

Setting $$A=0$$ makes every trial exponent independent of transient melt-pool depth:

$$n=B,\qquad k=2^B.$$

The scan path and run controls expose the remaining renderer tokens:

```text
Mode  X      Y      Z  Power      Param
1     0.000  0.000  0  0          0
0     0.002  0.000  0  <<power>>  <<velocity>>
```

```foam
endTime       <<endTime>>;
writeInterval <<writeInterval>>;
```

The driver converts the full D4σ diameter from micrometres to a 2σ radius in metres and assigns it to both lateral dimensions. Measured and simulated depths use the same normalization:

$$x=\frac{d}{D_{4\sigma}/2}
=\frac{d}{\min(\mathrm{dimensions}_x,\mathrm{dimensions}_y)}.$$

It also converts speed to m/s, renders the track, computes its travel time, and uses that time for both run-control placeholders.

### Calibration configuration

The supplied `config.yml` specifies:

- Trial values $$B=0,1,\ldots,9$$ for each experimental condition.
- Bounds $$0\leq B\leq9$$ and at most ten simulations per condition.
- A 1 µm depth-match tolerance and 0.15 local-posterior standard-deviation tolerance.
- A 1,000-point PCHIP surrogate grid.
- 2,000 retained posterior draws after 1,000 tuning steps, using eight sampling cores.
- 1,000 global uncertainty-propagation refits.
- Automatic PDF and CSV report generation.

Because the ten initial values equal the ten-simulation limit, this tutorial runs a fixed grid of 50 AdditiveFOAM cases.

### Campaign stages

For each process condition, `calibrateHeatSource`:

1. Creates `campaign/cases/<condition>/B<trial>/` from `template/`.
2. Renders the 2σ lateral dimension, trial value, power, velocity, end time, and write interval.
3. Runs `blockMesh`, `decomposePar`, and parallel `additiveFoam` through the case `Allrun` script.
4. Reads the liquidus CSV written by `meltPoolDimensions` and records the largest depth in the scan.
5. Updates `campaign/simulations.yml` immediately so an interrupted campaign can resume.
6. Builds the PCHIP response surrogate and samples the local posterior.

After all local calibrations are current, it performs the robust uncertainty-weighted global fit and generates the report.

<figure class="documentation-figure documentation-figure--plot">
  <img src="{{ '/assets/images/visualizations/heat-source-calibration-responses.png' | relative_url }}" alt="Five heat-source calibration response curves showing simulated melt-pool depth, measured replicate ranges, and inferred local shape parameters.">
  <figcaption>The five local response curves produced by the tutorial campaign. Each panel combines the cached AdditiveFOAM trials, the PCHIP surrogate used by the likelihood, the repeated measured depths, and the inferred local posterior mode.</figcaption>
</figure>

### Adapt to another data set

1. Copy the tutorial to a new writable directory and select a new `paths.campaign` location.
2. Replace `experiments.yml` with repeated depth measurements and complete process-condition metadata.
3. Enter the measured full D4σ diameter as `Spot_Size_microns`; the driver renders its 2σ radius into `dimensions`.
4. Update the case template's material, lateral source profile, absorption model, mesh, boundary conditions, and parallel decomposition.
5. Preserve each required placeholder exactly once and keep $$A=0$$ during local trials.
6. Run the campaign and inspect every local response before adopting the fitted source slope and intercept.

## Outputs

The completed directory has this structure:

```text
campaign/
├── cases/
│   └── P187p5_V500_D109p69/
│       ├── B0/
│       ├── B1/
│       └── ...
├── simulations.yml
├── calibration_state.yml
├── calibration_fit.yml
└── reports/
    ├── calibration_report.pdf
    └── calibration_summary.csv
```

| Output | Contents |
|---|---|
| `cases/.../log.run` | Complete solver output for one trial |
| `cases/.../metrics.yml` | Rendered inputs, selected isovalue, time and maximum length, width, and depth |
| `simulations.yml` | Trial values and simulated depths for every process condition |
| `calibration_state.yml` | Local posterior summaries and retained samples |
| `calibration_fit.yml` | Global coefficients, residual diagnostics, calibrated range, covariance, and bootstrap intervals |
| `calibration_summary.csv` | One row per experiment for downstream analysis |
| `calibration_report.pdf` | Response curves, measured depths, local posteriors, final fit, and numerical summary |

With `keep_successful: false`, the rendered dictionaries, log, post-processing data, and metrics remain, while processor and numeric time directories are removed after a successful trial.

### Apply the fitted coefficients

Open `campaign/calibration_fit.yml` and inspect `A`, `B`, the bootstrap intervals, `x_min`, `x_max`, and weighted residuals. Also inspect the local response and posterior plots in the PDF; a narrow-looking global fit does not repair a non-monotonic or poorly resolved AdditiveFOAM response curve.

<figure class="documentation-figure documentation-figure--plot">
  <img src="{{ '/assets/images/visualizations/heat-source-calibration-fit.png' | relative_url }}" alt="Final heat-source calibration fit with local posterior intervals and empirical uncertainty band.">
  <figcaption>The global relation fitted to the five local calibrations using x = d/(D4σ/2), the same normalized depth used by the projected heat-source models. The orange interval is obtained by drawing from the stored local posteriors and repeating the robust global fit.</figcaption>
</figure>

The coefficients are used in a production projected source as:

```foam
projectedGaussianCoeffs
{
    dimensions (54.845e-6 54.845e-6 30e-6);
    A          <source-slope>;
    B          <source-intercept>;
    transient  true;
    nPoints    (10 10 10);
}
```

The command uses `Spot_Size_microns: 109.69` as the full measured D4σ diameter, computes the 54.845 µm 2σ heat-source radius, and renders it into both runtime lateral dimensions. Measured and simulated depths are normalized with this same radius:

$$x=\frac{d}{D_{4\sigma}/2}
=\frac{d}{\min(\mathrm{dimensions}_x,\mathrm{dimensions}_y)}.$$

{: .important }
Insert `A` and `B` from `calibration_fit.yml` directly as `<source-slope>` and `<source-intercept>`.
