---
title: Installation
---

# Installation
[![OpenFOAM-10](https://img.shields.io/badge/OpenFOAM-10-blue.svg)](https://github.com/OpenFOAM/OpenFOAM-10)
AdditiveFOAM is built on source code released by the OpenFOAM Foundation [openfoam.org](https://openfoam.org/), which is available in public [OpenFOAM repositories](https://github.com/OpenFOAM). The current supported version is **OpenFOAM-10**.

## Spack install
[![Spack-Dev](https://img.shields.io/badge/Spack-Dev-blue.svg)](https://github.com/spack/spack-dev)

The easiest way to install AdditiveFOAM is using [spack](https://spack.readthedocs.io/en/latest/):  
```
spack install additivefoam
```
spack `develop` is currently required.

## Docker install
Alternatively, a Docker container with pre-built OpenFOAM-10 can be used:
```
docker pull openfoam/openfoam10-paraview510
docker run -it openfoam/openfoam10-paraview510
```

## Manual install
OpenFOAM-10 can be compiled from source code following the steps provided in the [OpenFOAM Foundation Documentation](https://openfoam.org/download/source/).

Once *OpenFOAM-10* is available on your system, perform the following steps:

1. Clone the AdditiveFOAM repository into the OpenFOAM project installation directory `WM_PROJECT_USER_DIR`:
   ```bash
   cd $WM_PROJECT_USER_DIR
   git clone https://github.com/ORNL/AdditiveFOAM.git
   ```

   If `git` is not available on your system (in the case of the OpenFOAM docker container) you can instead use:
   ```bash
   wget https://github.com/ORNL/AdditiveFOAM/archive/refs/heads/main.tar.gz
   mkdir AdditiveFOAM
   tar xzvf main.tar.gz -C AdditiveFOAM --strip-components=1
   ```
2. Build the `movingHeatSource` library and the `additiveFoam` executable:
   ```bash
   cd $WM_PROJECT_USER_DIR/AdditiveFOAM/applications/solvers/additiveFoam/movingHeatSource
   wmake libso
   cd $WM_PROJECT_USER_DIR/AdditiveFOAM/applications/solvers/additiveFoam
   wmake
   ```
