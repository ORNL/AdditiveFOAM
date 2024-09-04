<img src="https://raw.githubusercontent.com/ORNL-MDF/additivefoam-website-assets/main/images/logo.png" alt="image" width="800">

---
`AdditiveFOAM` is a free, open source heat and mass transfer software for simulations of Additive Manufacturing (AM) released by Oak Ridge National Laboratory. It is built upon OpenFOAM, a free, open source computational fluid dynamics (CFD) software package released by the OpenFOAM Foundation.

### Repository Features
| Link                                                | Description                              |
|-----------------------------------------------------------|------------------------------------------|
| [solver](applications/solvers/additiveFoam)               | Development version of the solver        |
| [utilities](applications/utilities)                       | Utilities for post-processing and code wrappers |
| [tutorials](tutorials)                                     | Tutorial cases |

## Documentation
[![Documentation Status][docs-badge]][docs-url]

The documentation for `AdditiveFOAM` is hosted on [GitHub Pages](https://ornl.github.io/AdditiveFOAM/).

## Installation and Dependencies
[![OpenFOAM-10](https://img.shields.io/badge/OpenFOAM-10-blue.svg)](https://github.com/OpenFOAM/OpenFOAM-10)

AdditiveFOAM is built on source code released by the OpenFOAM Foundation [openfoam.org](https://openfoam.org/), which is available in public [OpenFOAM repositories](https://github.com/OpenFOAM).

[![Spack-Dev](https://img.shields.io/badge/Spack-Dev-blue.svg)](https://github.com/spack/spack)

The easiest way to install AdditiveFOAM is using [spack](https://spack.readthedocs.io/en/latest/):  
```
spack install additivefoam
```

See the [installation instructions](https://ornl.github.io/AdditiveFOAM/docs/installation/#installation) in the [documentation](https://ornl.github.io/AdditiveFOAM/) for other options for building `AdditiveFOAM`.

## Citing
If you use AdditiveFOAM in your work, please cite the Zenodo DOI [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8034098.svg)](https://doi.org/10.5281/zenodo.8034098) of the version you used as a software citation:
```bibtex
@software{AdditiveFOAM_1.0.0,
  author       = {John Coleman and
                  Kellis Kincaid and
                  Gerald L. Knapp and
                  Benjamin Stump and
                  Alexander J. Plotkowski},
  title        = {AdditiveFOAM: Release 1.0},
  month        = jun,
  year         = 2023,
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.8034098},
  url          = {https://doi.org/10.5281/zenodo.8034098}
}
```

## License
[![GPL](https://img.shields.io/badge/GPL-3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

AdditiveFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file `LICENSE` in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

## Contact
For any questions, issues, or suggestions regarding AdditiveFOAM, you can reach out to the project maintainers through the GitHub repository's issue tracker or by contacting the development team directly.

## Contributing

We encourage you to contribute to AdditiveFOAM! Please check the
[guidelines](CONTRIBUTING.md) on how to do so.

We appreciate your interest in AdditiveFOAM and look forward to your contributions!

#### Contributors
- [John Coleman](https://www.ornl.gov/staff-profile/john-s-coleman)
- [Kellis Kincaid](https://www.ornl.gov/staff-profile/kellis-c-kincaid)
- [Gerry L. Knapp](https://www.ornl.gov/staff-profile/gerald-l-knapp)
- [Benjamin Stump](https://www.ornl.gov/staff-profile/benjamin-c-stump)
- [Alex Plotkowski](https://www.ornl.gov/staff-profile/alex-j-plotkowski)
- [Sam T. Reeve](https://www.ornl.gov/staff-profile/samuel-t-reeve)


[docs-badge]: https://img.shields.io/badge/docs-latest-brightgreen.svg
[docs-url]: https://ornl.github.io/AdditiveFOAM/
