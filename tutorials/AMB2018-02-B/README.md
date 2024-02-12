## AMB2018-02-B Tutorial Suite
Reference: https://www.nist.gov/ambench/amb2018-01-description
2.0 mm single track scan

### Standard Tutorial
The standard tutorial uses the `superGaussian` heat source model and can be run without modifications to any files.

### AMR Tutorial
The AMR tutorial incorporates adaptive mesh refinement coupled with the `refinementController` class implemented in AdditiveFOAM. For effective AMR, the [Zoltan](https://sandialabs.github.io/Zoltan) must be compiled with OpenFOAM. This can be accomplished by rerunning the ThirdParty Allwmake script with the [Zoltan 3.90 tarball](https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.90.tar.gz) present in the ThirdParty directory.

In addition to compiling Zoltan, two files must be changed from the standard tutorial. First, the included `dynamicMeshDictAMR` file should be renamed to `dynamicMeshDict`:

`$ mv constant/dynamicMeshDictAMR constant/dynamicMeshDict`

This will activate OpenFOAM's AMR capability. In additon, the `system/decomposeParDict` should be modified to use the `hierachical` decomposer with `zoltan` as the distributor. The `scotch` decomposer entry should be commented or removed. Be sure to include the `libzoltanDecomp.so` file and define appropriate entries for the decomposer and distributor. The `refinementHistory` constrain will ensure that all refined cells originating from a single parent cell remain on the same processor throughout the load distribution process.
