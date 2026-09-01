# Documentation visualizations

The documentation figures are rendered from reconstructed AdditiveFOAM tutorial results with ParaView's Python interface. The renderer fixes the camera, color ranges, legends, annotations, image dimensions, and AMR frame rate so the published assets can be regenerated without manual ParaView interaction.

## Required cases

Prepare six case directories and one completed calibration campaign from the AdditiveFOAM tutorials:

1. Run `AMB2018-02-B` without changing its configured refinement model.
2. Copy `AMB2018-02-B`, set `refinement/model` in `constant/heatSourceDict` to `timeStep`, confirm that `refineInterval` is `1` in `constant/dynamicMeshDict`, and run the case.
3. Run `multiLayerPBF` with its two-layer `Allrun` workflow.
4. Run `multiBeam`.
5. Run `nLightAFX`.
6. Run `tabulated`.
7. Run `heatSourceCalibration` with `calibrateHeatSource` so its `campaign/` directory contains `simulations.yml`, `calibration_state.yml`, `calibration_fit.yml`, and the report summary.

Run each tutorial's `Allrun` script so the parallel fields are reconstructed. Create an empty `case.foam` marker in each reconstructed case directory. For the multi-layer input, use the parent `multiLayerPBF` directory produced by `reconstructLayers`; it contains the temperature history across all layers.

## Render

### Regenerate everything

Place the six reconstructed cases and the completed calibration under one directory using their standard tutorial names:

```text
documentation-figures/
├── AMB2018-02-B/
├── AMB2018-02-B-AMR/
├── multiBeam/
├── multiLayerPBF/
├── nLightAFX/
├── tabulated/
└── heatSourceCalibration/
    ├── experiments.yml
    └── campaign/
```

Set the location once and regenerate and verify all 22 generated documentation
assets with no per-case paths:

```bash
export ADDITIVEFOAM_DOC_CASES="$FOAM_RUN/documentation-figures"
python3 scripts/visualization/regenerate_all.py
```

The AdditiveFOAM and OpenFOAM environments must be sourced so that
`plotPower`, `plotDimensions`, and `plotCET` are on `PATH` and the shared
`nLightAFX-1000.cfg` location is available. The driver also requires
`pvbatch`, FFmpeg, Pillow, Matplotlib, NumPy, Pandas, PyYAML, and SciPy. It invokes
each ParaView asset group in a fresh off-screen process, regenerates the
quantitative and analytic plots, writes the results to
`assets/images/visualizations`, and verifies the expected output set. Use
`--output-dir` to write elsewhere or `--dry-run` to inspect every command.
Alternatively, pass `--cases-root` once on the command line. Individual case
arguments are available only as overrides for nonstandard directory layouts. Use
`--amb2018-data-case` if the completed AMB power log, melt-pool CSVs, and
solidification CSV are stored separately from the reconstructed baseline case.
The driver validates these inputs before replacing the three quantitative plots.
Use `--calibration-campaign` and `--calibration-experiments` when those outputs
are stored elsewhere.

The case directories and their computed fields remain external inputs; the
driver does not copy them into the documentation repository.

### Render only the ParaView assets

ParaView 5.10 on Ubuntu installs its Python modules under `/usr/lib/python3/dist-packages`. From the documentation repository, run:

```bash
PYTHONPATH=/usr/lib/python3/dist-packages \
QT_QPA_PLATFORM=offscreen \
pvbatch --force-offscreen-rendering \
  scripts/visualization/render_documentation.py \
  --baseline-case /path/to/AMB2018-02-B \
  --multi-beam-case /path/to/multiBeam \
  --nlight-case /path/to/nLightAFX \
  --tabulated-case /path/to/tabulated \
  --amr-case /path/to/AMB2018-02-B-AMR \
  --multi-layer-case /path/to/multiLayerPBF \
  --output-dir assets/images/visualizations
```

The script writes a temperature animation and 300 dpi poster for every tutorial, together with `multi-layer-fields.png`, `amr-refinement.png`, and `amr-refinement.mp4`. The tutorial videos are `amb2018-temperature.mp4`, `multi-beam-temperature.mp4`, `multi-layer-temperature.mp4`, `nlight-afx-temperature.mp4`, and `tabulated-temperature.mp4`. FFmpeg is required to encode the 1600×600 H.264 videos. Encoding uses the compatible High profile with a near-lossless constant-rate factor so contour lines, labels, and color gradients remain sharp; video resolution is defined in pixels rather than dpi.

On headless systems whose OpenGL driver produces incomplete framebuffer tiles after several render views, add `--only NAME` and run each asset group in a fresh `pvbatch` process. Valid names are `amb2018`, `multi-beam`, `nlight-afx`, `tabulated`, `multi-layer`, `amr`, and `multi-layer-fields`.

## AMB2018 output plots

For the quantitative tutorial figures, enable `meltPoolDimensions` and `solidificationData` in the AMB2018-02-B `system/controlDict`, run `Allrun`, and invoke the plotting utilities:

```bash
plotPower AMB2018-02-B \
  -o assets/images/visualizations/amb2018-power.png
plotDimensions AMB2018-02-B \
  -o assets/images/visualizations/amb2018-melt-pool-dimensions.png
plotCET AMB2018-02-B \
  -o assets/images/visualizations/amb2018-cet.png
```

These figures are produced from `log.additiveFoam`, `postProcessing/meltPoolDimensions/*.csv`, and `solidificationData/solidification-data.csv`, respectively.

## Heat-source model comparison

Generate the normalized model comparison from the documented analytic distributions and the measured table in the `tabulated` tutorial:

```bash
python3 scripts/visualization/plot_heat_source_models.py \
  --tabulated-profile "$ADDITIVEFOAM_TUTORIALS/tabulated/constant/beamProfile.txt" \
  --nlight-config "$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg" \
  --output assets/images/visualizations/heat-source-models.png \
  --projection-output assets/images/visualizations/heat-source-projection.png \
  --metrics-output assets/images/visualizations/profile-metrics.png
```

The script reads the nLight AFX Index 6 coefficients from the shared configuration and the selected tabulated-profile file, which is the measured Index 3 table in the tabulated tutorial. `heat-source-models.png` distinguishes volumetric models from `projected` profile-plus-projection compositions. `profile-metrics.png` shows a rotated elliptical centroid and principal `D4Sigma` axes, the area-equivalent reference width, and the separate reference-aspect-ratio and applied-depth data paths. The outputs are saved at 300 dpi.

Generate the Kelly absorption-model figure from the parameters used in the example dictionary:

```bash
python3 scripts/visualization/plot_kelly_absorption.py \
  --output assets/images/visualizations/kelly-absorption.png
```

## Heat-source calibration plots

Generate the local AdditiveFOAM response curves and final uncertainty-propagated fit from a completed calibration campaign:

```bash
python3 scripts/visualization/plot_heat_source_calibration.py \
  --experiments /path/to/heatSourceCalibration/experiments.yml \
  --simulations /path/to/campaign/simulations.yml \
  --summary /path/to/campaign/reports/calibration_summary.csv \
  --state /path/to/campaign/calibration_state.yml \
  --fit /path/to/campaign/calibration_fit.yml \
  --responses-output assets/images/visualizations/heat-source-calibration-responses.png \
  --fit-output assets/images/visualizations/heat-source-calibration-fit.png
```

The response figure reconstructs each PCHIP surrogate from the stored AdditiveFOAM trials and overlays the measured replicate range and local posterior mode. The fit figure uses `x = 2*Depth / selected diameter`, resamples the stored local posteriors with the campaign seed, repeats the uncertainty-weighted `soft_l1` fit, and plots the empirical 95% response interval. The generated CSV uses `x_2depth_over_diameter`; the fit must declare `x_definition: 2*Depth / selected diameter` and provide `nSlope` and `nIntercept`. The `areaEquivalent` tutorial uses the 109.69 µm area-equivalent `D4Sigma` as its selected diameter.
