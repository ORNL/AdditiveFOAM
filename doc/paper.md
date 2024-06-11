---
title: 'AdditiveFOAM: A Continuum Heat and Mass Transfer Code for Laser-Based Additive Manufacturing'
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
  - name: Alex Plotkowski
    affiliation: "3"
  - name: Gerry L. Knapp
    affiliation: "3"
  - name: Benjamin Stump
    affiliation: "1"
affiliations:
 - name: Computational Science and Engineering Division, Oak Ridge National Laboratory
   index: 1
 - name: Nuclear Energy and Fuel Cycle Division, Oak Ridge National Laboratory
   index: 2
 - name: Material Science and Technology Division, Oak Ridge National Laboratory
   index: 3
date: 06 June 2024
bibliography: paper.bib
---

# Summary

Additive manufacturing (AM) processes present novel methods to quickly develop and produce complex components. The physics involved in AM are often complex, particularly in the case of laser powder bed fusion (LBPF). In LPBF, a laser beam scans over a thin layer of powder to melt a desired pattern in two dimensions. More powder is added over this layer and more scans are conducted, resulting in a consolidated three-dimensional part after many layers. This process inherently involves physical phenomena across length and time scales, ranging from phenomena within the melt pool that occur at the micrometer and microsecond length and time scales to bulk phenomena that occur over meter length scales for a large part [@Sing2020]. Due to the wide range of process parameters that can be adjusted, characterizing new materials, process variations, and part geometries through experimental work alone can be expensive and time-consuming. Computational models can be used to provide mechanistic understanding of the process to guide experiments and obtain more optimal results from the LPBF process with less expense.

# Statement of need

Numerical models for LPBF traditionally seek to represent physics on a single length scale. This results in two categories of models: 1) highly resolved but computationally intensive melt pool models which are intended to simulate only a few millimeters of scan path, or 2) coarse, approximate models meant to simulate length scales of an entire part. The first category contains models which seek to resolve all relevant physics within the melt pool, including fluid flow, vapor cavity formation, and powder particles. Examples of such codes include ALE3D [@noble2017ale3d]. These models can produce excellent agreement with experiments and offer insights into the underlying causes of anomalous features during processing, but are generally quite expensive to run, limiting the length of scan that can be simulated to a few millimeters with current computational resources. Due to the computational expense, these models are well-suited for use at specialized research instutions with large HPC infrastructure to answer target questions about melt pool scale physical phenomena. The second category is largely made up of continuum heat transfer and thermomechanical models. This includes analytical models solving the linear heat transfer equation (see e.g. @Stump2019) which can rapidly simulate the bulk heat transport during the LPBF process, and slightly more complex heat transfer models that numerically solve the energy equation with a moving heat source. These models can provide insight on the thermal history of a part, but their approximate representation of the melt pool offers limited accuracy when attempting to link this data the melt pool features needed to predict secondary quantities such as microstructural evolution.

There is a clear need for an intermediate framework which is able to accurately predict melt pool shape, thermal gradients, and other meso-scale quantities across an entire part. These predictions can be provided to other modeling frameworks to predict grain structure, material properties, thermal stresses, and other important qualtities which affect the strength and durability of the component. Some multi-physics commercial packages exist, such as FLOW-3D [@FLOW-3D] which focuses on the melt pool scale phenomena and the ANSYS Additive Simulation tools that target either the melt pool scale or part scale simulations. These tools typically focus on one end of spectrum of length and time scales and are also not open for custom modification. AdditiveFOAM was concieved to provide an open-source solution for numerical models of LPBF at this middle-of-the-road grade of accuracy and scale. A comprehensive, physics-based heat source library allows for accurate predictions of the melt pool shape and thermal gradients. An included optional fluid dynamics solver may increase accuracy in some scenarios, at the cost of computational expense. A dual implicit/explicit solver leverages numerical conditions to efficiently solve the equation set, and built-in functionalities allow interfacing with other physics-based modeling tools, notably ExaCA [@rolchigo2022exaca], which can be used to predict microstructural evolution from the thermal history generated during an AdditiveFOAM simulation. These qualities position AdditiveFOAM to be a useful tool for researchers to generate thermal data for LPBF processes in a range of conditions and applications and extend this data to make useful predictions about the performance of components in service. Because AdditiveFOAM is built on the OpenFOAM framework, it can also provide a test bed for developing and testing the effects of novel numerical methods or physics.

# Acknowledgements

This manuscript has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy (DOE). The publisher, by accepting the article for publication, acknowledges that the United States Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this manuscript, or allow others to do so, for United States Government purposes. The DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan.

This research used resources of the Oak Ridge Leadership Computing Facility (OLCF), supported by DOE under contract DE-AC05-00OR22725.

# References
