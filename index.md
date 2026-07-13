---
layout: versioned
title: Home
nav_order: 1
permalink: /
description: An open-source CFD code for additive manufacturing built on OpenFOAM.
---

<h1 class="additivefoam-wordmark">
  <img src="{{ '/assets/images/AdditiveFOAM-wordmark.svg' | relative_url }}" alt="AdditiveFOAM">
</h1>

An open-source CFD code for additive manufacturing built on OpenFOAM.
{: .fs-6 .fw-300 }

[Install AdditiveFOAM]({{ '/docs/installation/' | relative_url }}){: .btn .btn-primary .fs-5 }
[Run your first case]({{ '/docs/quick-start/' | relative_url }}){: .btn .fs-5 }
[Source code](https://github.com/ORNL/AdditiveFOAM){: .btn .fs-5 }

## Capabilities

- Transient conduction, convection, melting, solidification, latent heat, buoyancy, mushy-zone drag, and Marangoni flow.
- Multiple simultaneous moving beams with configurable spatial profiles.
- Temperature-dependent solid, liquid, and powder properties and tabulated solid-fraction paths.
- Dynamic load-balanced adaptive mesh refinement (AMR) around thermal regions and moving heat sources.
- Single-track and multi-layer powder-bed-fusion workflows.
- Runtime melt-pool measurements, solidification-condition extraction, and thermal histories for [ExaCA](https://github.com/LLNL/ExaCA).

## Requirements

AdditiveFOAM 2.0 requires the OpenFOAM Foundation release **OpenFOAM 13**. Start with the [installation guide]({{ '/docs/installation/' | relative_url }}) and then copy a validated [tutorial]({{ '/tutorials/' | relative_url }}) instead of constructing a case from an empty directory.

## Contributors

- [John Coleman](https://www.ornl.gov/staff-profile/john-s-coleman)
- [Kellis Kincaid](https://www.ornl.gov/staff-profile/kellis-c-kincaid)
- [Gerry L. Knapp](https://www.ornl.gov/staff-profile/gerald-l-knapp)
- [Benjamin Stump](https://www.ornl.gov/staff-profile/benjamin-c-stump)
- [Alex Plotkowski](https://www.ornl.gov/staff-profile/alex-j-plotkowski)
- [Matt Rolchigo](https://www.ornl.gov/staff-profile/matt-rolchigo)
- [Sam T. Reeve](https://www.ornl.gov/staff-profile/samuel-t-reeve)

## Citation

If AdditiveFOAM contributes to published work, cite the [JOSS article](https://doi.org/10.21105/joss.07770) and the version-specific software record. The [Zenodo concept DOI](https://doi.org/10.5281/zenodo.8034097) resolves to the latest archived release.

<div class="citation-badges">
  <a href="https://doi.org/10.21105/joss.07770"><img src="{{ '/assets/images/joss-doi.svg' | relative_url }}" alt="JOSS article DOI: 10.21105/joss.07770"></a>
  <a href="https://doi.org/10.5281/zenodo.8034097"><img src="{{ '/assets/images/zenodo-doi.svg' | relative_url }}" alt="Zenodo DOI: 10.5281/zenodo.8034097"></a>
</div>

## License and support

AdditiveFOAM is distributed under the GNU General Public License v3 or later. Report defects and request features through the [GitHub issue tracker](https://github.com/ORNL/AdditiveFOAM/issues).
