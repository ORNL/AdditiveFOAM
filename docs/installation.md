---
title: Installation
nav_order: 4
---

# Installation
[![OpenFOAM-10](https://img.shields.io/badge/OpenFOAM-10-blue.svg)](https://github.com/OpenFOAM/OpenFOAM-10)

AdditiveFOAM is built on source code released by the OpenFOAM Foundation [openfoam.org](https://openfoam.org/), which is available in public [OpenFOAM repositories](https://github.com/OpenFOAM). The current supported version is **OpenFOAM-10**.

---

## Spack install
[![Spack-Dev](https://img.shields.io/badge/Spack-Dev-blue.svg)](https://github.com/spack/spack)

The easiest way to install AdditiveFOAM is using [spack](https://spack.readthedocs.io/en/latest/):  
```
spack install additivefoam
```
spack `develop` is currently required.

### Spack installation with external modules

If you are installing OpenFOAM through spack, some configuration is needed to enable spack to use external modules, which is common on HPC systems. For OpenFOAM, this is largely needed for the C++ compiler and MPI dependencies, which are both generally available as preconfigured modules that can be loaded for functionality that has been configured for a specific HPC system. While module versions and names may vary from system to system, this example provides a template for configuring spack on your system.

First, the desired system compiler needs to be added to the spack compiler configuration.

```shell
module load gcc/8.1.0
spack compiler add $(dirname $(which gcc))
```

Second, you must get the location prefix of the OpenMPI installation by copying the output of the below command.

```shell
module load openmpi/3.1.5
echo $(dirname $(dirname $(which mpirun)))
```

> `/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/3.1.5/centos7.5_gnu8.1.0`

Next, open the package configuration file (~/.spack/packages.yaml)to add in system modules.

```bash
spack config edit packages
```

This will need to be edited to specify the external dependencies and specifications for the spack module. An example of `~/.spack/packages.yaml` is provided below. This file specifies that all packages should use the spack compiler (gcc@8.1.0) that was added in the first step. It also specifies that the spack module openmpi@3.1.5 will require the external `openmpi/3.1.5` and `gcc/8.1.0` modules from the system.

```yaml
packages:
  all:
    compiler: [gcc@8.1.0]
  openmpi:
    externals:
    - spec: openmpi@3.1.5
      prefix: /software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/3.1.5/centos7.5_gnu8.1.0
      modules:
      - openmpi/3.1.5
      - gcc/8.1.0
```

Once "packages.yaml" is configured, run `spack install openmpi@3.1.5` to install the package. You should get a message that the modules have external modules `openmpi/3.1.5` and `gcc/8.1.0`.

Finally, run the following to install OpenFOAM and additivefoam:

```bash
module load gcc/8.1.0 openmpi/3.1.5
spack load openmpi@3.1.5
spack install additivefoam
```

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
