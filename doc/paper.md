---
title: 'AdditiveFOAM 1.0: A Continuum Physics Code for Additive Manufacturing'
tags:
  - Computational Fluid Dynamics
  - Heat Transfer
  - Mass Transfer
  - Laser-based Additive Manufacturing
  - OpenFOAM
authors:
  - name: John Coleman
    affiliation: "1"
    corresponding: true
  - name: Kellis Kincaid
    affiliation: "2"
  - name: Gerry L. Knapp
    affiliation: "3"
  - name: Benjamin Stump
    affiliation: "1"
  - name: Alex Plotkowski
    affiliation: "1"
affiliations:
 - name: Computational Science and Engineering Division, Oak Ridge National Laboratory
   index: 1
 - name: Nuclear Energy and Fuel Cycle Division, Oak Ridge National Laboratory
   index: 2
 - name: Material Science and Technology Division, Oak Ridge National Laboratory
   index: 3
date: 12 June 2024
bibliography: paper.bib
---

# Summary

AdditiveFOAM is a computational framework that simulates transport phenomena in Additive Manufacturing (AM) processes, built on OpenFOAM [@openfoam], the leading free, open-source software package for computational fluid dynamics (CFD). OpenFOAM offers an extensible platform for solving complex multiphysics problems using state-of-the-art finite volume methods. AdditiveFOAM utilizes the capabilities of OpenFOAM to develop specialized tools aimed at addressing challenges in AM processing.

Metal additive manufacturing, also known as metal 3D printing, is an advanced manufacturing technique that creates physical parts from a three-dimensional (3D) digital model by melting metal powder or wire feedstock. A significant area of research in metal AM is focused on process planning to mitigate anomalous features during printing that are deleterious to part performance (e.g. porosity and cracking), as well as controlling localized microstructure and material properties. Given the high costs and substantial time requirements associated with experimental methods for qualifying new materials and processes, there is a compelling incentive for researchers to utilize advanced numerical models. In this context, AdditiveFOAM offers modeling capabilities to better understand undesirable features in printing, thereby enhancing process planning and reducing the reliance on labor-intensive experimental campaigns.

# Statement of need

Numerical models for AM traditionally seek to represent physics at a single length scale depending on the target question. Two main categories of models with available software solutions are distinguished: 1) high-fidelity melt pool-scale models, and 2) part-scale heat transfer models. 

The first category contains models which seek to resolve all relevant physics within the melt pool. These models can produce excellent agreement with experiments and offer insights into the underlying causes of anomalous features during processing, but are generally expensive to run, limiting the scan length that can be simulated to a few millimeters with current computational resources. Due to their computational expense, these models are well-suited for use at specialized research institutions with large HPC infrastructure to answer target questions about melt pool scale physical phenomena. An exhaustive review of high-fidelity melt pool models is beyond the current scope of this article; however, two primary software packages are highlighted: 

- ALE3D [@ALE3D] is a versatile multiphysics simulation tool that uses the Arbitrary Lagrangian-Eulerian approach and has been used for powder-resolved simulations of laser-material interactions in AM. However, ALE3D is a limited access code for use by DoD and DOE laboratories and their contractors and ALE3D4I is available for U.S. companies and academics through individual use agreements. 

- FLOW-3D [@FLOW-3D] is a commercial CFD software known for its capabilities in simulating complex free-surface problems, and has several specialized models for AM. However, the Flow-3D is proprietary, meaning its source code is not available for public inspection or modification and users must purchase a license to use the software. 

The second category contains models that simplify melt pool physics to rapidly simulate heat transport, residual stress, and distortion across an entire part. Commercial finite element thermomechanics software solutions built on Abaqus [@abaqus] and Ansys [@ansys] are available; however, these software packages are not open and free to use and develop upon.   

There is an established need for intermediate frameworks that are capable of accurately predict melt pool shape, thermal gradients, and other meso-scale quantities across an entire part. Some open-source mesoscale heat transfer software tools have recently emerged to address these challenges: 

- Adamantine [@adamantine] is a thermomechanics simulation tool built on an open-source software stack designed for high performance computing across computing architectures, including the deal.II finite element software package [@deal.ii].

- GO-MELT [@go-melt] is a gpu-optimized multidomain thermal simulation tool implemented as a finite element solver using jax [@jax].  

These software tools are promising platforms for further development of GPU-capable mesoscale simulations of AM processes, however at present, they remain limited to single gaussian heat sources and do not consider fluid flow or general thermodynamic pathways which may be important for accurate predictions of the conditions that lead to anomalous features and microstructure evolution during printing. 

# Software Features
AdditiveFOAM features a volumetric heat source library that supports any number of independently moving sources with distinct energy profiles that have been validated for a number of AM processes. The volumetric source term is used in the solution of the unsteady, nonlinear heat equation which considers three distinct phases (solid, liquid, and powder) as a continuum. A novel thermodynamic algorithm enables variable-order time integration for any thermodynamic pathway supported through tabulated solid fraction-temperature lookup tables. The available time integration schemes are explicit forward Euler, implicit backward Euler, Backward differentiation formula (BDF-2), and Crank-Nicolson. An adaptive time integration and time stepping method automatically switches between implicit and explicit schemes to balances the computational cost and solution accuracy depending on the problem state and user defined tolerances. Additionally, AdditiveFOAM features boundary conditions for Marangoni driven fluid flow in the melt pool as well as convective and radiative heat transfer. Finally, AdditiveFOAM leverages existing OpenFOAM libraries to sample thermal data needed for microstructure predictions with a specific library for coupling to the grain structure prediction software ExaCA [@exaca]. Future releases of AdditiveFOAM will focus on resource-optimized adaptive mesh refinement (AMR) as well as dynamic load-balancing using the Zoltan [@zoltan] library.

AdditiveFOAM was designed to be used openly by researchers, industry, and academics. It has already been used in a number of scientific publications primarily focused on studying how AM process planning influences the development of local microstructural features in a part. Select publications that demonstrates AdditiveFOAM’s usage towards understanding the connection between process planning and microstructure evolution in metal AM include: 

- Knapp et al. [@knapp] developed an automated framework for statistical calibration of the gaussian heat source parameters available in AdditiveFOAM to accurately predict the melt pool shape, solidification conditions, and solidification grain structure for single track welds on IN625 plate. 

- Rolchigo and co-workers [@exaca; @texture-regimes] used AdditiveFOAM to generate thermal data for solidification grain structure predictions is laser powder bed fusion, validated against the NIST AM Benchmark dataset [@amb2018]. 

- Haines et al. [@haines-1] used AdditiveFOAM to confirm experimentally observed phase transformation pathways in 17–4 PH stainless steel during the laser powder bed fusion (L-PBF) process, and [@haines-2] 

- Plotkowski et al. [@plotkowski] used AdditiveFOAM to confirm experimentally observed grain growth in Fe-Si steel during the Laser Engineered Net Shaping (LENS) process. 

The free, open-source nature of AdditiveFOAM will continue to enable scientific explorations in AM processing. The source code for AdditiveFOAM has been archived to Zenodo with the link DOI: [@AdditiveFOAM_1.0.0]. Documentation for AdditiveFOAM is hosted on GitHub pages with the link: https://ornl.github.io/AdditiveFOAM/. 

# Acknowledgements

This manuscript has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy (DOE). The publisher, by accepting the article for publication, acknowledges that the United States Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this manuscript, or allow others to do so, for United States Government purposes. The DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan. 

This research used resources of the Oak Ridge Leadership Computing Facility (OLCF), supported by DOE under contract DE-AC05-00OR22725. 

This research used resources of the Compute and Data Environment for Science (CADES) at the Oak Ridge National Laboratory, supported by the Office of Science of the U.S. Department of Energy under Contract No. DE-AC05-00OR22725. 

# References
