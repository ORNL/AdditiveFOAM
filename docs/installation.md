---
layout: versioned
title: Installation
parent: User Guide
nav_order: 1
permalink: /docs/installation/
---

# Installation

AdditiveFOAM 2.0 requires the OpenFOAM Foundation release **OpenFOAM 13**.

## Prerequisites

- OpenFOAM 13 built or installed according to the [Foundation instructions](https://openfoam.org/download/13-source/).
- Python 3 with Matplotlib, Pandas, and NumPy.
- ParaView for visualizing OpenFOAM fields.
- Optional: Zoltan support in OpenFOAM for [dynamic load-balanced AMR]({{ '/docs/amr/#mesh-refinement-and-load-balancing' | relative_url }}).

## Install from main

Source OpenFOAM first:

```bash
source /path/to/OpenFOAM-13/etc/bashrc
echo "$WM_PROJECT $WM_PROJECT_VERSION"
```

The output must identify OpenFOAM version 13. Then clone AdditiveFOAM and load its environment:

```bash
git clone https://github.com/ORNL/AdditiveFOAM.git
cd AdditiveFOAM
source etc/bashrc
```

`etc/bashrc` verifies the OpenFOAM version and sets:

- `ADDITIVEFOAM_PROJECT_DIR`
- `ADDITIVEFOAM_APPLICATIONS`
- `ADDITIVEFOAM_TUTORIALS`
- `ADDITIVEFOAM_ETC`
- `$ADDITIVEFOAM_BIN` on `PATH`

Build all libraries, the solver, and utilities from the repository root:

```bash
./Allwmake
```

## Verify the installation

```bash
command -v additiveFoam
command -v createScanPath
command -v primesToAdditiveFoam
echo "$ADDITIVEFOAM_VERSION"
```

Run `additiveFoam -help`. Solver startup also prints the AdditiveFOAM version, Git description, and commit SHA recorded at build time.

## Shell startup

Source both environments in this order in every new shell:

```bash
source /path/to/OpenFOAM-13/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
```

{: .custom }
Add those commands to the appropriate shell startup file only after confirming they work interactively.
