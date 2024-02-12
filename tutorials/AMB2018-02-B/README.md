## AMB2018-02-B Tutorial Suite
A 2 mm single track case based on [AM-BENCH](https://www.nist.gov/ambench/amb2018-02-description)

### Standard Tutorial
The standard tutorial uses the `superGaussian` heat source model and can be run without modifications to any files.

### AMR Tutorial
The AMR tutorial incorporates adaptive mesh refinement via with the `refinementController` class implemented in AdditiveFOAM. For CPU load-balanced AMR, the [Zoltan](https://sandialabs.github.io/Zoltan) must be compiled with OpenFOAM. This can be accomplished by rerunning the ThirdParty Allwmake script with the [Zoltan 3.90 tarball](https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.90.tar.gz) present in the ThirdParty directory.

In addition to compiling Zoltan, the following changes must be made from the standard tutorial:
1) The `constant/dynamicMeshDictAMR` file should be renamed to `constant/dynamicMeshDict`:

`$ cp -r constant/dynamicMeshDictAMR constant/dynamicMeshDict`

3) In the `constant/heatSourceModels` file, change the `refinementController` entry from `none` to `uniformIntervals`.

4) In the `system/decomposeParDict` file, use the `hierachical` decomposer combined with the `zoltan` distributor. The `scotch` decomposer entry should be commented or removed. Be sure to include the `libzoltanDecomp.so` file and define appropriate entries for the decomposer and distributor. The `refinementHistory` constrain will ensure that all refined cells originating from a single parent cell remain on the same processor throughout the load distribution process.
