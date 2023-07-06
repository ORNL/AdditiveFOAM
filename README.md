# AdditiveFOAM
AdditiveFOAM is a free, open source heat and mass transfer software for simulations of Additive Manufacturing (AM) released by Oak Ridge National Laboratory. It is built upon OpenFOAM, a free, open source computational fluid dynamics (CFD) software package released by the OpenFOAM Foundation.

## Citing
If you use AdditiveFoam in your work, please cite the [source code](CITATION.bib). Also, please consider citing relevant AdditiveFOAM **Publications** listed below.

## Contributors
- [John Coleman](https://www.ornl.gov/staff-profile/john-s-coleman)
- [Kellis Kincaid](https://www.ornl.gov/staff-profile/kellis-c-kincaid)
- [Gerry L. Knapp](https://www.ornl.gov/staff-profile/gerald-l-knapp)
- [Benjamin Stump](https://www.ornl.gov/staff-profile/benjamin-c-stump)
- [Alex Plotkowski](https://www.ornl.gov/staff-profile/alex-j-plotkowski)

## Publications
Some select publications using AdditiveFOAM are provided:
1. [Coleman et al. "Sensitivity of Thermal Predictions to Uncertain Surface Tension Data in Laser Additive Manufacturing", J. Heat Transfer (2020) HT-19-1539](https://asmedigitalcollection.asme.org/heattransfer/article/doi/10.1115/1.4047916/1085538/Sensitivity-of-Thermal-Predictions-to-Uncertain)
2. [Knapp et al. "Calibrating uncertain parameters in melt pool simulations of additive manufacturing", Comp. Mat. Sci. (2023) 111904](https://www.sciencedirect.com/science/article/abs/pii/S0927025622006152)
3. [Rolchigo et al. "ExaCA: A performance portable exascale cellular automata application for alloy solidification modeling", Comp. Mat. Sci. (2022) 111692](https://www.sciencedirect.com/science/article/abs/pii/S0927025622004189)

## Repository Features
| Link                                                | Description                              |
|-----------------------------------------------------------|------------------------------------------|
| [solver](applications/solvers/additiveFoam)               | Development version of the solver        |
| [utilities](applications/utilities)                       | Utilities for post-processing and code wrappers |
| [tutorials](tutorials)                                     | Tutorial cases based on [NIST AMB2018](https://www.nist.gov/ambench/amb2018-02-description) single tracks |

## Installation
AdditiveFOAM is built on source code released by the OpenFOAM Foundation [openfoam.org](https://openfoam.org/), which is accessible to the public through the [OpenFOAM source code repositories at Github](https://github.com/OpenFOAM). The current supported version is **OpenFOAM-10**, which can be compiled from source code following the steps provided by the [OpenFOAM Foundation Documentation](https://openfoam.org/download/source/).

Once **OpenFOAM-10** is compiled, perform the following steps:

1. Clone the AdditiveFOAM repository into the OpenFOAM project installation directory `WM_PROJECT_USER_DIR`:
   ```bash
   cd $WM_PROJECT_USER_DIR
   git clone https://github.com/ORNL/AdditiveFOAM.git
   ```
2. Build the `movingHeatSource` library and the `additiveFoam` executable:
   ```bash
   cd $WM_PROJECT_USER_DIR/AdditiveFOAM/applications/solvers/additiveFoam/movingHeatSource
   wmake libso
   cd $WM_PROJECT_USER_DIR/AdditiveFOAM/applications/solvers/additiveFoam
   wmake
   ```
## Documentation
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

      - `transportProperties`: Sets the transport properties of the material. The thermal conductivity **kappa** and specific heat **Cp** are given as temperature dependent second-order polynomials for each phase in the material.
      
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
      
         Each heat source model has the ability to be update the depth of the heat source for keyhole modeling, by setting the **transient** flag to **True** and defining an **isoValue** to track the depth of an isotherm contained within the heat source radius. An example of this usage is provided in the [multiBeam](tutorials/multiBeam) tutorial.

   - `0/`: Contains the initial fields. The available fields are provided in the files:
      - `T`:            temperature
      - `U`:            velocity
      - `p_rgh`:        reduced pressure
      - `alpha.solid`:  solid volume fraction
      - `alpha.powder`: powder volume fraction

   - `system/`: Contains simulation configuration files.
      - `controlDict`: Set simulation time settings and numerical parameters.
      - `fvSchemes`: Set the discretization schemes used to solve the governing equations
      - `fvSolution`: Set solution algorithms and convergence criterias. Note: fluid flow can be turned off by setting **nOuterCorrectors** to **0** in the **PIMPLE** dictionary.
      
3. Run the simulation:
An example run script which creates a mesh, decomposes the mesh across multiple processors, and runs the AdditiveFOAM case in parallel using MPI is provided in tutorial `Allrun` script.

4. Visualize and post-process the results using **ParaView**
   ```bash
   touch case.foam
   paraview case.foam
   ```
### Creating Scan Path Files
AdditiveFOAM supports a scan path file format that decomposes the laser path into segments that are either a) line sources or b) point sources.

| Column   | Description                                                                                                                 |
|----------|-----------------------------------------------------------------------------------------------------------------------------|
| Column 1 | mode = 0 for line source, mode = 1 for point source                                                                               |
| Columns 2-4 | (x,y,z) coordinates in meters. <br>For a line (mode = 0), this is the final position of the line raster. <br>For a spot (mode = 1), this is the constant position of the laser |
| Column 5 | Value for laser power in watts                                                                                               |
| Column 6 | For a line (mode = 0), this is the velocity of the laser in meters/second. <br>For a point (mode = 1), this is the dwell time of the laser in seconds                         |

### Exporting ExaCA Data
One feature of AdditiveFOAM is its ability to export thermal information to [ExaCA](https://github.com/LLNL/ExaCA), a cellular automata (CA) code for grain growth under additive manufacturing conditions. This feature is enabled using the **execute** flag in the `constant/foamToExaCADict` file. The solidification conditions at the specified **isotherm** is tracked in the represenative volume element defined by **box** and a resolution defined by **dx**. It is recommended to track the liquidus isotherm. Users should be warned that this interpolation may cause a significant load-balancing issues if the resolution of the ExaCA data is much finer than that of the AdditiveFOAM mesh, and therefore, this feature should be used selectively.

## License
AdditiveFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file `LICENSE` in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

## Contact
For any questions, issues, or suggestions regarding AdditiveFOAM, you can reach out to the project maintainers through the GitHub repository's issue tracker or by contacting the development team in the **Contributors** links above.

We appreciate your interest in AdditiveFOAM and look forward to your contributions!
