if (adjustTimeStep)
{
    // Interface Courant number
    const volScalarField::Internal alphaCoNumvf
    (
        0.5*fvc::surfaceSum(mag(phi))()
       *pos0(alpha1() - 0.01)*pos0(0.99 - alpha1())
       /mesh.V()*runTime.deltaT()
    );
    const scalar alphaCoNum = gMax(alphaCoNumvf);

    Info<< "Interface Courant Number max: " << alphaCoNum << endl;

    const scalar maxAlphaCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", 1.0);

    const scalar maxDeltaTFact =
        min(maxCo/(CoNum + small), maxAlphaCo/(alphaCoNum + small));

    const scalar deltaTFact =
        min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    scalar deltaT = min(deltaTFact*runTime.deltaTValue(), maxDeltaT);

    #include "capillaryDeltaT.H"

    sources.adjustDeltaT(deltaT);

    runTime.setDeltaT(deltaT);
}

Info<< "deltaT = " << runTime.deltaTValue() << endl;
