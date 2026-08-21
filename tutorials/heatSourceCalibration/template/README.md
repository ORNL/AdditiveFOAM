# Calibration case template

This is the OpenFOAM-14 AdditiveFOAM case copied for every trial value in the
heat-source calibration campaign. It uses a `projectedHeatSource` with a
`tabulated` circular Gaussian profile, an `exponential` projection, and SS316L.
The measured D4Sigma diameter is 109.69 microns.

The calibration command expects this full case structure:

```text
0/
constant/
system/
Allrun
Allclean
```

The required renderer placeholders are:

- `constant/heatSourceDict`: `<<nIntercept>>`
- `constant/scanPath`: `<<power>>`, `<<velocity>>`
- `system/controlDict`: `<<endTime>>`, `<<writeInterval>>`

The normalized 67-by-67 profile is provided directly in
`constant/beam_profile.txt`. It uses 2.5-micron spacing and spans -82.5 through
+82.5 microns in both lateral directions. The model calculates the profile
integral and D4Sigma during initialization and uses `minimumDepth` as the
initial projected source depth.

The projected closure is `n = clip(nSlope*log2(a) + nIntercept, 0, 9)` and the
axial decay exponent is `k = 2^n`.

The SS316L material configuration supplies `thermoPath`, emissivity, and the
Marangoni coefficient. Transient source depth and `meltPoolDimensions` obtain
their default isovalues from that material configuration.
