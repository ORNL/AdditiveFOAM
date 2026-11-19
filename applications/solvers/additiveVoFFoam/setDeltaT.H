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

    // The write-time adjustment of setDeltaT divides the time to the next
    // write evenly and replaces the time step, so the scan-path adjustment
    // is applied to its result
    runTime.setDeltaT(deltaT);

    deltaT = runTime.deltaTValue();
    sources.adjustDeltaT(deltaT);
    runTime.setDeltaTNoAdjust(deltaT);
}

Info<< "deltaT = " << runTime.deltaTValue() << endl;
