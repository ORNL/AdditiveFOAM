---
title: 'AdditiveFOAM: A Continuum Heat and Mass Transfer Code for Laser-Based Additive Manufacturing'
tags:
  - Computational Fluid Dynamics
  - Heat Transfer
  - Mass Transfer
  - Laser-based Additive Manufacturing
  - OpenFOAM
authors:
  - name: [John Coleman]
    affiliation: "1"
    corresponding: true
  - name: [Kellis Kincaid]
    affiliation: "2"
  - name: [Alex Plotkowski]
    affiliation: "3"
  - name: [Gerry L. Knapp]
    affiliation: "3"
  - name: [Benjamin Stump]
    affiliation: "1"
affiliations:
 - name: [Computational Science and Engineering Division, Oak Ridge National Laboratory]
   index: 1
 - name: [Nuclear Energy and Fuel Cycle Division, Oak Ridge National Laboratory]
   index: 2
 - name: [Material Science and Technology Division, Oak Ridge National Laboratory]
   index: 3
date: [March 27, 2024]
bibliography: additiveFoam-joss-1.bib
---

# Summary

Additive manufacturing (AM) processes present novel methods to quickly develop and produce complex components. The physics involved in AM are often complex, particularly in the case of laser powder bed fusion (LBPF). In LPBF, a laser beam scans over a thin layer of powder to melt a desired pattern in two dimensions. More powder is added over this layer and more scans are conducted, resulting in a consolidated three-dimensional part after many layers. This process inherently involves physical phenomena across length and time scales, ranging from the micrometer size of the melt pool up to the order of meters for a large part [@Sing2020]. Numerical models of LPBF traditionally seek to represent physics on a single length scale, resulting in either highly resolved but computationally intensive melt pool models or lightweight but inaccurate simplified models of an entire part.

AdditiveFOAM is an open-source heat and momentum transfer model which serves to bridge the gap between these scales. It is built upon the OpenFOAM framework [@Weller1998], a package ubiquitous in computational fluid dynamics and multiphysics simulations. AdditiveFOAM relies on the meshing, discretization, and solution algorithms developed for OpenFOAM, a choice allowing for increased attention to be paid to the implementation of the thermal model. Two features in particular make AdditiveFOAM attractive for use in LPBF simulations: first, it will automatically switch between implicit and explicit solution methods, resulting in highly efficient computation; second, it contains a generalized heat source modeling framework capable of modeling any number of laser beams, each with unique shape, properties, and scan path. These features are designed to facilitate integration with the next generation of LPBF systems.

# Statement of need

Numerical models for LPBF are traditionally divided into one of two categories. The first category contains models which seek to resolve all relevant physics within the melt pool, including fluid flow, vapor cavity formation, and powder particles. Examples of such codes include ALE3D [@noble2017ale3d]. These models can produce excellent agreement with experiments and offer insights into the underlying causes of defect formation, but are generally quite expensive to run, limiting the length of scan that can be simulated to a few millimeters with current computational resources. This has limited the adoption of these models to specialized research instutions with large HPC infrastructure, rather than widespread industrial and academic application. At the other end of the spectrum are an array of simple heat transfer models, including analytical models (see e.g. @Stump2019) which can simulate an entire part in seconds, and slightly more complex counterparts which numerically solve the energy equation with a simple moving heat source. These models can provide insight on the thermal history of a part, but their approximate representation of the melt pool offers limited accuracy when this data is extended to predict secondary quantities including defect formation or microstructural evolution.

There is a clear need for an intermediate framework which is able to accurately predict melt pool shape, thermal gradients, and other meso-scale quantities across an entire part. These predictions can be provided to other modeling frameworks to predict grain structure, material properties, thermal stresses, and other important qualtities which affect the strength and durability of the component. While some commercial packages exist for this purpose, such as in FLOW-3D [@FLOW-3D], their implemenation is inherently opaque and the functionality and compatiblity with other code is typically limited. AdditiveFOAM was concieved to provide an open-source solution for numerical models of LPBF at this middle-of-the-road grade of accuracy and scale. A comprehensive, physics-based heat source library allows for accurate predictions of the melt pool shape and thermal gradients. An included optional fluid dynamics solver may increase accuracy in some scenarios, at the cost of computational expense. A dual implicit/explicit solver leverages numerical conditions to efficiently solve the equation set, and built-in functionalities allow interfacing with other physics-based modeling tools, notably ExaCA [@rolchigo2022exaca], which can be used to predict microstructural evolution from the thermal history generated during an AdditiveFOAM simulation. These qualities position AdditiveFOAM to be a useful tool for researchers to generate thermal data for LPBF processes in a range of conditions and applications, and extend this data to make useful predictions about the performance of components in service.

# Architecture and Methods

## Solver

The core off AdditiveFOAM is contained within a solver loop reminiscent of most OpenFOAM applications. The PIMPLE [@Moukalled2016] algorithm is used to advance the solution in time, and for pressure-velocity coupling when fluid flow is enabled. An energy equation, shown below, predicts the variation of temperature in space and time. The second and third terms on the right hand side account for the release of latent heat during phase change, and the power input from the scanning laser(s). The latter term is calculated using a generic heat source modeling library.

$$ \rho cp \left(\frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T\right) = \nabla \cdot (k \nabla T) + \rho L{f}\frac{\partial f_{s}}{\partial t} + Q $$

This equation is discretized and solved using the finite volume tools provided in OpenFOAM. The energy equation is solved iteratively until a certain tolerance is met, after which the physical properties are updated and the time loop continues. The user retains full control of the exact discretization and solution schemes, as is standard in the OpenFOAM framework.

## Heat Source Model

The main solver interfaces with the heat source models via the `movingHeatSource` library, which is perhaps the most notable feature of AdditiveFOAM. This library produces an effective volumetric heat source which is can represent an arbitrary number of independent, unique lasers. Each individual laser is represented by classes which prescribe the shape, absorption, and path of the beam in time. The `heatSourceModel` class determines the spatial distribution of power for each beam. Shape sub-models are implemented in derived classes which presently include a standard `superGaussian` distribution, as well as a `modifiedSuperGaussian` developed by John Coleman et al. [REF]. The fraction of prescribed laser power which is actually absorbed by the material is controlled by the `absorptionModel` class. Again, derived classes implement sub-models where the absorption can either be `constant` or be calculated dynamically using the `Kelly` model. In the latter model, melt pool geometry and internal reflections are considered to produce a more accurate prediction of absorption as detailed in [@kelly1966equation]. The transient effect of melt pool depth on the absorptivity can also be activitated by setting the `transient` flag in the heat source dictionary. The movement and power of each beam are controlled by the `scanPath` file, which uses a standardized file format for LPBF. This information is processed and relayed to the heat source model by the `movingBeam` library. Each model interfaces with the main `movingHeatSource` class, which returns the volumetric heat field to the solver with a simple `qDot()` function call.

## Boundary Conditions:
In addition to the standard boundary conditions provided in OpenFOAM, we provide two custom boundary conditions necessary for simulating AM processes.
The ```maragoni``` boundary condition is provided to calculate the shear stress induced by surface tensions gradients at along a boundary. The condition equates the viscous stress to the Marangoni stress:
$$
\mu \frac{\partial \mathbf{u}}{\partial \mathbf{n}} = \frac{\partial \gamma}{\partial T}\frac{\partial T}{\partial \mathbf{\tau}}
$$

Its implementation is derived from the ```partialSlip``` boundary conditions available in OpenFOAM-10, which uses two basic functions to update the boundary condition. The first function is ```Foam::marangoniFvPatchVectorField::snGrad()``` which sets the surface-normal gradient of the velocity field at the specified boundary faces. In OpenFOAM, ```snGrad()``` is notably used in the discretization of the laplacian operator in the momentum equation. Mathematically, the return of the function call is

$$
[\mathbf{I} - \mathbf{n}\cdot\mathbf{n}] \cdot\mathbf{u}_{i} - \mathbf{u}_{i} + c[\mathbf{I} - \mathbf{n}\cdot\mathbf{n}] \cdot \nabla T_{i}
$$

and is represented in the source code as:
```
146    // calculate the surface normal gradient on the patch
147    return
148    (
149        transform(I - sqr(nHat), pif) - pif
150     + coeff_*transform(I - sqr(nHat), tGrad) / this->patch().deltaCoeffs()
151   )*this->patch().deltaCoeffs();
```
The first two terms uses vector transormations to subtract out the component of the velocoity field normal to the surface which may be numerically nonzero, and the last term equations the surface-normal gradient to the product of the maragoni coefficient $c=\frac{1}{\mu}\frac{\partial \gamma}{\partial T}$ and the component of the temperature gradient that is tangential to the boundary.

The second function equated in the boundary condition is ```Foam::marangoniFvPatchVectorField::evaluate()```

$$
\mathbf{u}_{b} = [\mathbf{I} - \mathbf{n}\cdot\mathbf{n}] \cdot\mathbf{u}_{i}
$$

```
vectorField::operator=(transform(I - sqr(nHat), pif));
```

The ```mixedTemperature``` boundary condition is provided to calculate the mixed convective and radiative heat transport across a boundary and has a standard mathematical form.


## Implicit / Explicit Solution:

At each time step, the solver calculates a diffusion number, a non-dimensionalized quantity which relates the variation in temperature to the local mesh and time step size. In the case that the global maximum of the diffusion number is less than `1.0`, the solver will automatically perform an explicit integration of the energy equation, drastically reducing the computation time associated with this operation. This switching can be disabled by the user via the `explicitSolve` flag in the `fvSolution` file.


# Tutorials and Support
AdditiveFOAM's documentation includes tutorials based on standard test cases designed to illustrate the usage of various features within the code. These tutorials include a basic line scan with and without solution of the momemtum equation, a multi-beam example to illustrate the `movingHeatSourceModel` framework, and a multi-layer build demonstrating the capacity of the code to perform on large scales. A basic continuous integration pipeline ensures that the code continues to compile and run after changes are made in the repository.

# Installation
AdditiveFOAM can be installed using spack, simplifying the deployment process. For environments where spack is not suitable, a Docker container with pre-built OpenFOAM-10 is available, demonstrating the software’s commitment to accessibility. Detailed installation instructions in the README facilitate easy setup and configuration.

# Acknowledgements

This manuscript has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy (DOE). The publisher, by accepting the article for publication, acknowledges that the United States Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this manuscript, or allow others to do so, for United States Government purposes. The DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan.

This research used resources of the Oak Ridge Leadership Computing Facility (OLCF), supported by DOE under contract DE-AC05-00OR22725.

# References
