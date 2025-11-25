# Based on this build for OpenFOAM 12
#   - https://github.com/NixOS/nixpkgs/blob/07ba4b68bad1931618851055077ed20d677b81ed/pkgs/by-name/op/openfoam-org/package.nix
# Useful as well
#   - https://git.computecanada.ca/nix/ccpkgs/-/blob/cc-20.09/pkgs/openfoam.nix

{
  stdenv,
  bash,
  m4,
  flex,
  bison,
  fftw,
  scotch,
  boost,
  openmpi,
  cgal,
  zlib,
  fetchFromGitHub,
  lib,
  trilinos-mpi
}:
let
  openfoamDrv = stdenv.mkDerivation rec {
    pname = "openfoam-org";
    version = "10";
    src = fetchFromGitHub {
      owner = "OpenFOAM";
      repo = "OpenFOAM-12";
      rev = "refs/tags/version-${version}";
      hash = "sha256-1vcBZELsThlfSJeW3iFm8sTh+uOgKKEYU0g+XlxczdA=";

    };
    
    nativeBuildInputs = [
      bash
      flex
      bison
    ];
    
    
    buildInputs = [
      m4
      fftw
      boost
      cgal
      zlib
      trilinos-mpi
      openmpi
      scotch
    ];

    propagatedBuildInputs = [
      openmpi
    ];

    sourceRoot = ".";

    patchPhase = ''
    runHook prePatch

    mkdir -p builduser/OpenFOAM
    mkdir -p builduser/.OpenFOAM
    export HOME=$(pwd)/builduser
    export FOAM_DIR=$HOME/OpenFOAM/OpenFOAM-${version}
    mv source $FOAM_DIR
    mkdir -p $FOAM_DIR/paraviewout
    
    set +e
    for f in \
        $FOAM_DIR/wmake/scripts/* \
        $FOAM_DIR/wmake/*
    do
      substituteInPlace $f --replace-quiet /bin/bash ${bash}/bin/bash
    done
    set -e

    rm $FOAM_DIR/etc/config.sh/bash_completion
    touch $FOAM_DIR/etc/config.sh/bash_completion

    echo "set +e" | cat $FOAM_DIR/etc/bashrc > tmp
    rm $FOAM_DIR/etc/bashrc
    mv tmp $FOAM_DIR/etc/bashrc

    echo "set +e" | cat $FOAM_DIR/Allwmake > tmp
    rm $FOAM_DIR/Allwmake
    mv tmp $FOAM_DIR/Allwmake

    alias wmRefresh="placeholder"
    chmod +x $FOAM_DIR/Allwmake
    chmod +x $FOAM_DIR/applications/utilities/postProcessing/graphics/PVReaders/Allwmake
    touch $HOME/.OpenFOAM/prefs.sh

    # only if version 10
    substituteInPlace $FOAM_DIR/etc/bashrc --replace-fail "export WM_PROJECT_VERSION=dev" "export WM_PROJECT_VERSION=10"
    sed -i '47 i libDir=${openmpi.dev}/lib' $FOAM_DIR/etc/config.sh/mpi
    sed -ie 's|SCOTCH_ARCH_PATH=.*$|SCOTCH_ARCH_PATH=${scotch.dev}|' $FOAM_DIR/etc/config.sh/scotch

    set +e
    for f in \
        $FOAM_DIR/src/parallel/decompose/*/Make/options
    do
      substituteInPlace $f --replace-quiet /usr/include/scotch ${scotch.dev}/include
      substituteInPlace $f --replace-quiet "\$(SCOTCH_ARCH_PATH)/lib" ${scotch.out}/lib
    done
    set -e

    substituteInPlace $FOAM_DIR/wmake/rules/General/mplibOPENMPI --replace-fail "\$(MPI_ARCH_PATH)/include" "${openmpi.dev}/include"
    substituteInPlace $FOAM_DIR/wmake/rules/General/mplibOPENMPI --replace-fail "\$(MPI_ARCH_PATH)/lib" "${openmpi.out}/lib"

    # This has a broken link
    rm $FOAM_DIR/tutorials/mesh/snappyHexMesh/iglooWithFridges

    runHook postPatch

  '';

    buildPhase = ''
    runHook preBuild

    cd $FOAM_DIR
    source ./etc/bashrc

    # Can't get it working without these paths.
    # Shouldn't be required
    export LD_LIBRARY_PATH="${flex}/lib''${LD_LIBRARY_PATH}"
    export C_INCLUDE_PATH="${flex}/include''${C_INCLUDE_PATH}"
    export CPLUS_INCLUDE_PATH="${flex}/include''${CPLUS_INCLUDE_PATH}"

    ./Allwmake -j $NIX_BUILD_CORES -q

    runHook postBuild
  '';

    installPhase = ''
    runHook preInstall

    FOAM_OUT=$out/opt/OpenFOAM-${version}
    mkdir -p $FOAM_OUT

    cp -r bin $FOAM_OUT
    cp -r platforms $FOAM_OUT
    cp -r etc $FOAM_OUT
    cp -r applications $FOAM_OUT
    cp -r src $FOAM_OUT
    cp -r doc $FOAM_OUT
    cp -r tutorials $FOAM_OUT
    cp -r wmake $FOAM_OUT

    runHook postInstall
  '';

    passthru = {
      BASHRC = "${openfoamDrv}/opt/OpenFOAM-${version}/etc/bashrc";
    };
  };
in
  openfoamDrv
