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

3. Run the simulation:

   An example run script which creates a mesh, decomposes the mesh across multiple processors, and runs the AdditiveFOAM case in parallel using MPI is provided in tutorial `Allrun` script.

4. Visualize and post-process the results using **ParaView**
   ```bash
   touch case.foam
   paraview case.foam
   ```
   
## AdditiveFOAM File Structure

### Case Files

```
|-- case-directory
|   |-- 0
|   |   |-- alpha.powder
|   |   |-- alpha.solid
|   |   |-- p_rgh
|   |   |-- T
|   |   |-- U
|   |
|   |-- constant
|   |   |-- heatSourceDict
|   |   |-- scanPath
|   |   |-- transportProperties
|   |   |-- thermoPath
|   |
|   |-- system
|   |   |-- blockMeshDict
|   |   |-- controlDict
|   |   |-- fvScheme
|   |   |-- fvSolution
|   |
|   +-- ..
|
```

   - `0/`: Contains the initial fields. The available fields are provided in the files:
      - `T`:            temperature
      - `U`:            velocity
      - `p_rgh`:        reduced pressure
      - `alpha.solid`:  solid volume fraction
      - `alpha.powder`: powder volume fraction

  - `constant/`: Contains configuration and settings that define geometric and material conditions, including:

  - `transportProperties`: Sets the transport properties of the material. The thermal conductivity `kappa` and specific heat `Cp` are given as temperature dependent second-order polynomials for each phase in the material.
  
     The available phases are:
     - solid
     - liquid
     - powder
  
     The remaining properties are all assumed constant throughout the simulation.

  - `heatSourceDict`: Defines number, types, and paths of moving heat sources in the simulation.

     The available heat sources are:
     - superGaussian
     - modifiedSuperGaussian

     The available absorption models are:
     - constant
     - [Kelly](https://opg.optica.org/ao/fulltext.cfm?uri=ao-5-6-925&id=14272)
  
     Each heat source model has the ability to be update the depth of the heat source for keyhole modeling, by setting the `transient` flag to `True` and defining an `isoValue` to track the depth of an isotherm contained within the heat source radius. An example of this usage is provided in the [multiBeam](tutorials/multiBeam) tutorial.

   - `system/`: Contains simulation configuration files.
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
