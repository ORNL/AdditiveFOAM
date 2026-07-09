<img src="https://raw.githubusercontent.com/ORNL-MDF/additivefoam-website-assets/main/images/AdditiveFOAM-wordmark.svg" alt="image">

---
`AdditiveFOAM` is a free, open source heat and mass transfer software for simulations of Additive Manufacturing (AM) released by Oak Ridge National Laboratory. It is built upon OpenFOAM, a free, open source computational fluid dynamics (CFD) software package released by the OpenFOAM Foundation.

## Documentation
[![Documentation Status][docs-badge]][docs-url]

The documentation for `AdditiveFOAM` is hosted on [GitHub Pages](https://ornl.github.io/AdditiveFOAM/).

### Repository Features
| Link                                                | Description                              |
|-----------------------------------------------------------|------------------------------------------|
| [solver](applications/solvers/additiveFoam)               | Development version of the solver        |
| [tutorials](tutorials)                                     | Tutorial cases |

## Installation
[![OpenFOAM-13](https://img.shields.io/badge/OpenFOAM-13-blue.svg)](https://github.com/OpenFOAM/OpenFOAM-13)

AdditiveFOAM requires OpenFOAM-13 from the OpenFOAM Foundation. Install or
build OpenFOAM-13 first, then source the OpenFOAM environment:

```sh
source /path/to/OpenFOAM-13/etc/bashrc
```

Clone AdditiveFOAM, enter the repository, and source the AdditiveFOAM
environment:

```sh
git clone https://github.com/ORNL/AdditiveFOAM.git
cd AdditiveFOAM
source etc/bashrc
```

The AdditiveFOAM `etc/bashrc` checks that OpenFOAM-13 is active and sets the required paths.

Build all AdditiveFOAM libraries, solvers, and utilities with the master build
script from the repository root:

```sh
./Allwmake
```

For regular use, source both environments in each new shell or add them to your
shell startup file:

```sh
source /path/to/OpenFOAM-13/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
```

## Citing
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07770/status.svg)](https://doi.org/10.21105/joss.07770)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8034097.svg)](https://doi.org/10.5281/zenodo.8034097)

If you use AdditiveFOAM in your work, please cite the JOSS article and
consider citing the Zenodo DOI for the version used.

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
- [Matt Rolchigo](https://www.ornl.gov/staff-profile/matt-rolchigo)

[docs-badge]: https://img.shields.io/badge/docs-latest-brightgreen.svg
[docs-url]: https://ornl.github.io/AdditiveFOAM/
