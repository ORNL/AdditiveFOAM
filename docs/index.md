---
layout: versioned
title: User Guide
nav_order: 2
has_children: true
has_toc: false
permalink: /docs/
redirect_from:
  - /reference/
---

# User Guide

This User Guide explains how to install, run, configure, and extend AdditiveFOAM. It documents the governing equations, numerical methods, input files, heat-source models, boundary conditions, mesh refinement, solver controls, outputs, and utilities provided by AdditiveFOAM.

## Getting started

1. [Install AdditiveFOAM]({{ '/docs/installation/' | relative_url }}).
2. [Run the quick start]({{ '/docs/quick-start/' | relative_url }}).
3. Select the closest [tutorial]({{ '/tutorials/' | relative_url }}) and copy it as the starting point for your case.
4. Use the sections below to configure the physical models, numerical methods, and outputs required by the case.

## Table of contents

1. **[Installation]({{ '/docs/installation/' | relative_url }})**

   Install AdditiveFOAM with OpenFOAM 14 and configure its environment variables.

2. **[Quick Start]({{ '/docs/quick-start/' | relative_url }})**

   Copy, run, monitor, and visualize the AMB2018-02-B baseline case.

3. **[Governing Equations]({{ '/docs/theory/' | relative_url }})**

   Review the momentum, temperature, phase, and thermodynamic coupling equations solved by AdditiveFOAM.

4. **[Heat Source Models]({{ '/docs/heat-sources/' | relative_url }})**

   Define scan paths, absorption models, volumetric heat distributions, and spatial and temporal source integration.

5. **[Boundary Conditions]({{ '/docs/boundary-conditions/' | relative_url }})**

   Configure the Marangoni velocity and mixed convection-radiation temperature conditions.

6. **[Adaptive Mesh Refinement]({{ '/docs/amr/' | relative_url }})**

   Refine the mesh around thermal regions and moving heat sources and redistribute it between MPI processes.

7. **[Case Files and Materials]({{ '/docs/case-files/' | relative_url }})**

   Define initial fields, material properties, phase paths, scan paths, and the required case dictionaries.

8. **[Solver and Time Controls]({{ '/docs/solver-controls/' | relative_url }})**

   Select temperature-phase coupling, time discretization, fluid-flow coupling, linear solvers, and adaptive time stepping.

9. **[Function Objects]({{ '/docs/function-objects/' | relative_url }})**

   Write melt-pool dimensions, solidification conditions, and ExaCA temperature histories during a simulation.

10. **[Utilities]({{ '/docs/utilities/' | relative_url }})**

    Generate scan paths, run multi-layer cases, reconstruct parallel output, convert beam profiles, and plot results.
