---
title: Usage
nav_order: 5
usemathjax: true
---

# Usage

---

## Running a simulation
To run an AdditiveFOAM simulation, it is recommended to perform the following steps:
1. Prepare the case directory structure using a provided template:
   ```bash
   mkdir -p $FOAM_RUN/additivefoam
   cd $FOAM_RUN/additivefoam
   cp -r $WM_PROJECT_USER_DIR/AdditiveFOAM/tutorials/AMB2018-02-B userCase
   cd userCase
   ```

2. Modify the necessary input files according to your simulation requirements. These files are described in [Case Files](#case-files).

3. Run the simulation using shell scripting. Examples are provided in the `Allrun` scripts in each tutorial which creates a mesh, decomposes the mesh across multiple processors, and runs the AdditiveFOAM case in parallel using [MPI](https://www.mpi-forum.org/).

   ```bash
   #!/bin/sh
   cd ${0%/*} || exit 1    # Run from this directory
   
   # Source tutorial run functions
   . $WM_PROJECT_DIR/bin/tools/RunFunctions

   application=`getApplication`
   
   runApplication blockMesh
   runApplication decomposePar
   runParallel $application

   runApplication reconstructPar
   ```

4. Visualize and post-process the results using [ParaView](https://www.paraview.org/)
   ```bash
   touch case.foam
   paraview case.foam
   ```
   
## AdditiveFOAM File Structure

### Case Files

An AdditiveFOAM case directory has a typical organization like:
```
|-- case-directory/
|   |-- 0/
|   |   |-- alpha.powder
|   |   |-- alpha.solid
|   |   |-- p_rgh
|   |   |-- T
|   |   |-- U
|   |
|   |-- constant/
|   |   |-- heatSourceDict
|   |   |-- scanPath
|   |   |-- transportProperties
|   |   |-- thermoPath
|   |
|   |-- system/
|   |   |-- blockMeshDict
|   |   |-- controlDict
|   |   |-- fvScheme
|   |   |-- fvSolution
|   |
|   +-- ..
```

- `0`: Directory containing the initial conditions and boundary conditions for each field. The fields are specified by their file name:
  - `T`:            temperature             (required)
  - `U`:            velocity                (optional)
  - `p_rgh`:        reduced pressure        (optional)
  - `alpha.solid`:  solid volume fraction   (optional)
  - `alpha.powder`: powder volume fraction  (optional)

- `constant`: Directory containing the definitions for material conditions, including:

  - `transportProperties`: File that defines the transport properties of the material. An example of this file for IN625 is:
    ```cpp
    solid
    {
        kappa   (8.275  0.01472   0.0);
        Cp      (579.28 0.0       0.0);
    }

    liquid
    {
        kappa   (4.889    0.014743  0.0);
        Cp      (750.65   0.0       0.0);
    }

    powder
    {
        kappa   (-0.07707   0.00075   0.0);
        Cp      (747.568    0.0       0.0);
    }

    rho     [1 -3 0 0 0 0 0]    7569.92;   // reference density
    mu      [1 -1 -1  0 0 0 0]  0.003032;  // dynamic viscosity
    beta    [0 0 0 -1 0 0 0]    1.2e-4;    // coefficient of thermal expansion
    DAS     [0 1 0 0 0 0 0]     10e-6;     // dendrite arm spacing (for drag)
    Lf      [0  2 -2  0 0 0 0]  2.1754e5;  // latent heat of fusion
    ```
   The thermal conductivity `kappa` and specific heat `Cp` are given as temperature dependent second-order polynomials for the `solid`, `liquid`, and `powder` phases. The remaining transport properties are all assumed to be constant during the simulation.

  - `heatSourceDict`: File that defines the number and mathematical form of moving heat sources in the simulation. An example of this file is:
    ```cpp
    sources (beam);

    beam
    {
        pathName            scanPath;    

        absorptionModel     constant;
        
        constantCoeffs
        {
            eta             0.3;
        }
        
        heatSourceModel     superGaussian;
        
        superGaussianCoeffs
        {
            k               2.0;
            dimensions      (85.0e-6 85.0e-6 30e-6);
            nPoints         (10 10 10);
        }
    }
    ```

     The available entries for the absorption models are:
     - `constant`
     - [`Kelly`](https://opg.optica.org/ao/fulltext.cfm?uri=ao-5-6-925&id=14272): An effective absorption is calculated based of the geometry and apsect ratio of the heat source shape to simulate keyhole formation.

     The available entries for the heatSourceModel are:
     - `superGaussian`
     - `modifiedSuperGaussian`
       
     {: .custom }
     Each heat source model can dynamically update its depth to the a specified isotherm position by setting the `transient` flag to `True`. The recommended value for `isoValue` is the alloy liquidus temperature for simulating keyhole formation.

#### system/
Contains simulation configuration files.
  - `controlDict`: Set simulation time settings and numerical parameters.
  - `fvSchemes`: Set the discretization schemes used to solve the governing equations
  - `fvSolution`: Set solution algorithms and convergence criterias.
  
    {: .custom }
    Fluid flow can be turned off by setting `nOuterCorrectors` to `0` in the `PIMPLE` dictionary.


### Scan Path Files
AdditiveFOAM supports a scan path file format that decomposes the laser path into segments that are either a) line sources (mode = 0) or b) point sources (mode = 1).

An example scan path file for a bi-directional raster with a hatch-spacing of 100 $$\mu m$$ is:
```
Mode    X       Y       Z     Power   Parameter
1       0.000   0.000   0.0   0.0     0.0
0       0.002   0.000   0.0   195.0   0.8
0       0.002   1.e-4   0.0   0.0     0.8
0       0.000   1.e-4   0.0   195.0   0.8
```

The beam starts at (0, 0, 0) $$m$$ and moves horizontally to (0.002, 0, 0.0) $$m$$ with a power of 195 $$W$$ and speed of 0.8 $$m/s$$. Then, the laser turns off and moves vertically to (0.002, 0.0001, 0.0) $$m$$ at a speed of 0.8 $$m/s$$. Finally, the laser moves horizontally to (0.002, 0.0001, 0.0) $$m$$ with a power of 195 $$W$$ and speed of 0.8 $$m/s$$. At this point, the laser turns off for the remainder of the simulation.

A summary of the entries defining a path segment for the scan path is provided in the following table:

| Column   | Description                                                                                                                 |
|----------|-----------------------------------------------------------------------------------------------------------------------------|
| Column 1 | Mode = 0: line source, <br>Mode = 1: point source                                                                               |
| Columns 2-4 | Mode = 0: the final (x,y,z) position of the beam in meters. <br>Mode = 1: the current (x,y,z) position of the beam in meters |
| Column 5 | Value for laser power in watts                                                                                               |
| Column 6 | Mode = 0: the speed of the beam in meters/second. <br>Mode = 1: the time the beam remains at its current position in seconds                         |

## Exporting ExaCA Data
AdditiveFOAM is able to export thermal data to [ExaCA](https://github.com/LLNL/ExaCA), a cellular automata (CA) code for grain growth under additive manufacturing conditions. 

This feature is enabled using the `execute` flag in the `constant/foamToExaCADict` file. The solidification conditions at the specified `isotherm` is tracked in the represenative volume element defined by `box` and a resolution defined by `dx`. It is recommended to track the liquidus isotherm. Users should be warned that this interpolation may cause a significant load-balancing issues if the resolution of the ExaCA data is much finer than that of the AdditiveFOAM mesh, and therefore, this feature should be used selectively.
