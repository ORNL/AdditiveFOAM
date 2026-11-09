# additiveVoFFoam Tutorials

Both cases run a 2 mm single track in AlSi10Mg at 1000 W and 0.7 m/s with
`additiveVoFFoam`, using the measured nLight AFX profile of the
[tabulated](../additiveFoam/tabulated) heat-source tutorial. The metal and the
argon above it are separated by a volume of fluid interface, and the beam is
absorbed at the first metal surface it reaches. The heat-source dictionary is
that of the tabulated tutorial, and `Allrun` converts its PRIMES export.

`tabulated`

A 1.5 mm wide, 1 mm thick plate with 0.3 mm of gas above. The mesh is 10 um
within 0.35 mm of the track and coarsens towards the sides, the bottom and the
top.

`tabulatedPowderBed`

The same plate and mesh with a 50 um powder layer within 0.3 mm of the track.
`createPowderBed` writes `system/setFieldsDict` with non-overlapping particles
of 15-25 um radius, at least three cells across.

`smoothSetFields` sets `alpha.metal` from `system/setFieldsDict` with
volume fractions in the cells cut by the plate surface and the particles.

Material properties, including the surface tension and the vaporization data
for the recoil pressure, come from `$ADDITIVEFOAM_ETC/materials/AlSi10Mg.cfg`.
The `gas` dictionary in `constant/transportProperties` sets the properties of
the shielding gas.
