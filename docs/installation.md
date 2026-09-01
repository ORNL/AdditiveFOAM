---
layout: versioned
title: Installation
parent: User Guide
nav_order: 1
permalink: /docs/installation/
---

# Installation

AdditiveFOAM 2.0.0 requires the OpenFOAM Foundation release **OpenFOAM 14**.

## Prerequisites

- OpenFOAM 14 built or installed according to the [Foundation instructions](https://openfoam.org/download/14-source/).
- Python 3.10 through 3.12 when using the optional calibration and plotting utilities.
- ParaView for visualizing OpenFOAM fields.
- Optional: Zoltan support in OpenFOAM for [dynamic load-balanced AMR]({{ '/docs/amr/#mesh-refinement-and-load-balancing' | relative_url }}).

## Build AdditiveFOAM

Source OpenFOAM first:

```bash
source /path/to/OpenFOAM-14/etc/bashrc
echo "$WM_PROJECT $WM_PROJECT_VERSION"
```

The output must identify OpenFOAM version 14. Then clone AdditiveFOAM and load its environment:

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

To use the Python calibration and reporting utilities, create or activate a Python environment and install the pinned dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip check
```

`requirements.txt` pins Matplotlib, NumPy, Pandas, PyYAML, ReportLab, SciPy, and Seaborn. These packages are optional for the compiled solver and C++ utilities; they are required by `calibrateHeatSource` and the plotting commands.

The website-asset workflow has separate tooling requirements: ParaView's `pvbatch` renders reconstructed fields, FFmpeg encodes tutorial animations, and Pillow assembles rendered frames. They are not calibration dependencies.

The AdditiveFOAM environment adds the scripts to `PATH`; it does not create or activate the Python environment. Activate the intended Python environment in each new shell before invoking `calibrateHeatSource`.

## Verify the installation

```bash
command -v additiveFoam
command -v createScanPath
command -v primesToAdditiveFoam
command -v tabulatedProfileInfo
command -v calibrateHeatSource
echo "$ADDITIVEFOAM_VERSION"
```

Run `additiveFoam -help`. Solver startup also prints the AdditiveFOAM version, Git description, and commit SHA recorded at build time.

## Shell startup

Source both environments in this order in every new shell:

```bash
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
```

{: .custom }
Add those commands to the appropriate shell startup file only after confirming they work interactively.
