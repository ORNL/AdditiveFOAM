# additiveVoFFoam Tutorials

Both cases repeat the [tabulated](../additiveFoam/tabulated) heat-source tutorial, a 2 mm
single track in AlSi10Mg at 400 W and 0.8 m/s with a measured nLight AFX
profile, with `additiveVoFFoam`. The metal and the argon above it are separated
by a volume of fluid interface, and the beam is absorbed at the first metal
surface it reaches. The heat-source and scan-path dictionaries are those of the
tabulated tutorial, and `Allrun` converts its PRIMES export.

`tabulated`

The bare plate, 0.3 mm thick with 0.1 mm of gas above, on a 10 um mesh.

`tabulatedPowderBed`

A 50 um powder layer on the plate, on a 5 um mesh. `createPowderBed` writes
`system/setFieldsDict` with non-overlapping particles of 15 um mean radius.

`smoothSetFields` sets `alpha.metal` from `system/setFieldsDict` with
volume fractions in the cells cut by the plate surface and the particles.

Material properties, including the surface tension and the vaporization data
for the recoil pressure, come from `$ADDITIVEFOAM_ETC/materials/AlSi10Mg.cfg`.
The `gas` dictionary in `constant/transportProperties` sets the properties of
the shielding gas.
