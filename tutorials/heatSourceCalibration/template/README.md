# Calibration case template

This is the OpenFOAM-14 AdditiveFOAM case copied for every trial value in the
heat-source calibration campaign. It uses a `projectedGaussian` heat source
with SS316L. The laser D4sigma diameter was measured to be 109.69 microns,
giving a 2sigma heat-source radius of 54.845 microns.

The calibration command expects this full case structure:

```text
0/
constant/
system/
Allrun
Allclean
```

The required renderer placeholders are:

- `constant/heatSourceDict`: `<<B>>` and `<<heatSourceRadius>>`
- `constant/scanPath`: `<<power>>`, `<<velocity>>`
- `system/controlDict`: `<<endTime>>`, `<<writeInterval>>`

The projected Gaussian closure reads `<<B>>` directly. The
`$heatSourceRadius` alias supplies both lateral dimensions from one rendered
value.

The calibration command divides `Spot_Size_microns` by two, converts the
result to metres, assigns it to both lateral source dimensions, and uses the
same radius to normalize measured and simulated depth.

The SS316L material configuration supplies `thermoPath`, emissivity, and the
Marangoni coefficient. Transient source depth and `meltPoolDimensions` obtain
their default isovalues from that material configuration.
