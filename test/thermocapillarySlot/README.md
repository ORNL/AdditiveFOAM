# thermocapillarySlot

Verification of the thermocapillary (Marangoni) force against analytic
references, for both solvers, on the Sen and Davis (1982) slot: a 1 mm closed
slot holding 0.2 mm of liquid under 0.3 mm of gas, with a linear wall-to-wall
temperature difference giving a surface traction tau = 0.4 Pa, sigma = 0.01 N/m
and zero gravity (properties of Wirth et al. 2022, table 2).

| case | solver | checks |
|---|---|---|
| `af-100x20` | additiveFoam | flat-surface return flow under the `marangoni` slip condition against the exact finite-aspect Stokes solution |
| `vof-100x240` | additiveVoFFoam | free-surface deflection against the finite-aspect reference and the wall heights of Wirth et al. 2022, table 6 (y = Y/w = 0.188 / 0.213) |

## Running

In an OpenFOAM-14 environment with AdditiveFOAM built:

    python3 run.py <output directory> [--np N] [--endTime T]

This writes both cases, runs them and builds
`<output directory>/thermocapillarySlot.html` comparing AdditiveFOAM with the
analytic predictions. `--endTime 0.002` is a quick smoke test.

## Files

- `reference.py`: thin-slot solution, two-layer factor, and an exact
  finite-aspect Stokes solver (liquid with optional gas layer) giving the
  surface velocity, the normal-stress jump and the deflection.
- `slot.py`: case generation for both solvers, field readers and profiles.
- `plot.py`: the comparison page.
- `run.py`: driver.

## Status

- additiveFoam: second-order convergence to the exact flat-surface solution
  (velocity L2 error 0.9%, 0.2%, 0.05% on 50, 100, 200 cells).
- additiveVoFFoam: **not yet verified**. The free-surface case uses the
  100 x 240 mesh and interfaceCompression (MULES) transport of the original
  marangoni-test tutorial, the closest set-up to the benchmark found so far:
  wall heights 0.1880 / 0.2111 after 20 ms (reference 0.1877 / 0.2122).
  That run had not reached a steady state, and MULES runs on other meshes
  drifted with a widening interface band. The geometric schemes (MPLIC,
  MPLICU, PLIC, PLICU) are steady but give a smaller deflection (rise 0.82
  of the reference at 50 x 120, 0.77 at 100 x 240) that does not improve
  with refinement in x. Square meshes, on which the deflection spans under
  two cells, give 0.63 to 0.67 with either scheme.
