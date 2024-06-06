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
bibliography: paper.bib
---

# Summary

Additive manufacturing (AM) processes present novel methods to quickly develop and produce complex components. The physics involved in AM are often complex, particularly in the case of laser powder bed fusion (LBPF). In LPBF, a laser beam scans over a thin layer of powder to melt a desired pattern in two dimensions. More powder is added over this layer and more scans are conducted, resulting in a consolidated three-dimensional part after many layers. This process inherently involves physical phenomena across length and time scales, ranging from the micrometer size of the melt pool up to the order of meters for a large part [@Sing2020]. Characterizing new materials, process variations, and part geometries through experimental work alone is expensive and time-consuming, leading to an interest in the usage of computational models to predict various quantities associated with the LPBF process.

# Statement of need

Numerical models for LPBF traditionally seek to represent physics on a single length scale, resulting in either highly resolved but computationally intensive melt pool models which can only model a few millimeters of scan path, or lightweight but inaccurate simplified models of an entire part. The first category contains models which seek to resolve all relevant physics within the melt pool, including fluid flow, vapor cavity formation, and powder particles. Examples of such codes include ALE3D [@noble2017ale3d]. These models can produce excellent agreement with experiments and offer insights into the underlying causes of defect formation, but are generally quite expensive to run, limiting the length of scan that can be simulated to a few millimeters with current computational resources. This has limited the adoption of these models to specialized research instutions with large HPC infrastructure, rather than widespread industrial and academic application. At the other end of the spectrum are an array of simple heat transfer models, including analytical models (see e.g. @Stump2019) which can simulate an entire part in seconds, and slightly more complex counterparts which numerically solve the energy equation with a simple moving heat source. These models can provide insight on the thermal history of a part, but their approximate representation of the melt pool offers limited accuracy when this data is extended to predict secondary quantities including defect formation or microstructural evolution.

There is a clear need for an intermediate framework which is able to accurately predict melt pool shape, thermal gradients, and other meso-scale quantities across an entire part. These predictions can be provided to other modeling frameworks to predict grain structure, material properties, thermal stresses, and other important qualtities which affect the strength and durability of the component. While some commercial packages exist for this purpose, such as in FLOW-3D [@FLOW-3D], their implemenation is inherently opaque and the functionality and compatiblity with other code is typically limited. AdditiveFOAM was concieved to provide an open-source solution for numerical models of LPBF at this middle-of-the-road grade of accuracy and scale. A comprehensive, physics-based heat source library allows for accurate predictions of the melt pool shape and thermal gradients. An included optional fluid dynamics solver may increase accuracy in some scenarios, at the cost of computational expense. A dual implicit/explicit solver leverages numerical conditions to efficiently solve the equation set, and built-in functionalities allow interfacing with other physics-based modeling tools, notably ExaCA [@rolchigo2022exaca], which can be used to predict microstructural evolution from the thermal history generated during an AdditiveFOAM simulation. These qualities position AdditiveFOAM to be a useful tool for researchers to generate thermal data for LPBF processes in a range of conditions and applications, and extend this data to make useful predictions about the performance of components in service.

# Acknowledgements

This manuscript has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy (DOE). The publisher, by accepting the article for publication, acknowledges that the United States Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this manuscript, or allow others to do so, for United States Government purposes. The DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan.

This research used resources of the Oak Ridge Leadership Computing Facility (OLCF), supported by DOE under contract DE-AC05-00OR22725.

# References
