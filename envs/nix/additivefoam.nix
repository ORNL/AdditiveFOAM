{
  stdenv,
  makeWrapper,
  openfoam,
  openmpi,
  exaca,
  src,
  version,
  glibc,
  zlib
}:
stdenv.mkDerivation rec {
  pname = "additivefoam";
  inherit version;
  inherit src;

  nativeBuildInputs = [
    makeWrapper
  ];
  
  buildInputs = [
    openfoam
    openmpi
    zlib
  ];

  propagatedBuildInputs = [
    exaca
    openmpi
  ];
  
  buildPhase = ''
    runHook preBuild

    mkdir -p builduser/.OpenFOAM
    mkdir -p builduser/OpenFOAM
    export HOME=$(pwd)/builduser
    export USER=builduser

    source ${openfoam.BASHRC} || true

    APPS_DIR=$(pwd)/applications/solvers/additiveFoam

    cd $APPS_DIR/movingHeatSource
    wmake libso

    cd $APPS_DIR/functionObjects/ExaCA
    wmake libso

    cd $APPS_DIR
    wmake

    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall

    mkdir -p $out

    cp -r $FOAM_USER_APPBIN $out
    cp -r $FOAM_USER_LIBBIN $out

    cp -r ${src}/applications $out/
    cp -r ${src}/tutorials $out/

    runHook postInstall
  '';

  postInstall = ''
    wrapProgram $out/bin/additiveFoam \
      --suffix LD_LIBRARY_PATH : "$out/lib:${openmpi}/lib:${stdenv.cc.cc.lib}/lib:${zlib}/lib"

    LINKER=${glibc}/lib/ld-linux-x86-64.so.2;
    patchelf --set-interpreter $LINKER $out/bin/.additiveFoam-wrapped
  '';
  
  doCheck = false;
  doInstallCheck = true;

  installCheckPhase = ''

    cd $HOME
    mkdir -p app
    cp -r ${src}/tutorials/AMB2018-02-B/* app/
    chmod u+w -R app
    cd app

    # Ensure ExaCA is found
    substituteInPlace $HOME/app/Allrun --replace-fail "~/install/exaca/bin/ExaCA" "ExaCA"
    
    ./Allrun -withExaCA

    if [ -f "ExaCA/Output.vtk" ]; then
        echo "Check PASSED: Output.vtk generated."
    else
        echo "Check FAILED: Output.vtk not found."
        exit 1 
    fi

    cd ..
    rm -rf app

  '';
  
}
