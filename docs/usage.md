---
title: Usage
nav_order: 5
---

# Run AdditiveFOAM
To run an AdditiveFOAM simulation, it is recommended to perform the following steps:
1. Prepare the case directory structure using a provided template:
   ```bash
   mkdir -p $FOAM_RUN/additivefoam
   cd $FOAM_RUN/additivefoam
   cp -r $WM_PROJECT_USER_DIR/AdditiveFOAM/tutorials/AMB2018-02-B userCase
   cd userCase
   ```
2. Modify the necessary input files according to your simulation requirements. These files include:

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

   - `0/`: Contains the initial fields. The available fields are provided in the files:
      - `T`:            temperature
      - `U`:            velocity
      - `p_rgh`:        reduced pressure
      - `alpha.solid`:  solid volume fraction
      - `alpha.powder`: powder volume fraction

   - `system/`: Contains simulation configuration files.
      - `controlDict`: Set simulation time settings and numerical parameters.
      - `fvSchemes`: Set the discretization schemes used to solve the governing equations
      - `fvSolution`: Set solution algorithms and convergence criterias. Note: fluid flow can be turned off by setting `nOuterCorrectors` to `0` in the **PIMPLE** dictionary.
      
3. Run the simulation:
An example run script which creates a mesh, decomposes the mesh across multiple processors, and runs the AdditiveFOAM case in parallel using MPI is provided in tutorial `Allrun` script.

4. Visualize and post-process the results using **ParaView**
   ```bash
   touch case.foam
   paraview case.foam
   ```


## Creating Scan Path Files
AdditiveFOAM supports a scan path file format that decomposes the laser path into segments that are either a) line sources or b) point sources.

| Column   | Description                                                                                                                 |
|----------|-----------------------------------------------------------------------------------------------------------------------------|
| Column 1 | mode = 0 for line source, mode = 1 for point source                                                                               |
| Columns 2-4 | (x,y,z) coordinates in meters. <br>For a line (mode = 0), this is the final position of the line raster. <br>For a spot (mode = 1), this is the constant position of the laser |
| Column 5 | Value for laser power in watts                                                                                               |
| Column 6 | For a line (mode = 0), this is the velocity of the laser in meters/second. <br>For a point (mode = 1), this is the dwell time of the laser in seconds                         |

## Exporting ExaCA Data
One feature of AdditiveFOAM is its ability to export thermal information to [ExaCA](https://github.com/LLNL/ExaCA), a cellular automata (CA) code for grain growth under additive manufacturing conditions. This feature is enabled using the `execute` flag in the `constant/foamToExaCADict` file. The solidification conditions at the specified `isotherm` is tracked in the represenative volume element defined by `box` and a resolution defined by `dx`. It is recommended to track the liquidus isotherm. Users should be warned that this interpolation may cause a significant load-balancing issues if the resolution of the ExaCA data is much finer than that of the AdditiveFOAM mesh, and therefore, this feature should be used selectively.
