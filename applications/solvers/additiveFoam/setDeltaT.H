if (adjustTimeStep)
{
    // temporal resolution w.r.t freezing range
    scalar maxAlphaCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", 0.1);

    const scalar dtBydT = mag(runTime.deltaT() / (Tliq - Tsol)).value();

    alphaCoNum = gMax(mag(fvc::ddt(T))())*dtBydT;

    Info << "Courant Number thermo: " << alphaCoNum << endl;

    FoNum = gMin(Foam::pow(mesh.V().field(), 2.0/3.0) / (kappa/rho/Cp)());

    scalar maxFo = runTime.controlDict().lookupOrDefault<scalar>("maxFo", 1);

    scalar deltaT = min
    (
        min(maxCo/(CoNum + small), maxFo/(FoNum + small)),
        maxAlphaCo/(alphaCoNum + small)
    )*runTime.deltaTValue();

    deltaT = min
    (
        min(deltaT, runTime.deltaTValue() + 0.1*deltaT),
        1.2*runTime.deltaTValue()
    );

    deltaT = min(deltaT, maxDeltaT);

    // adjust time step to desired hit scan path resolution
    sources.adjustDeltaT(deltaT);

    runTime.setDeltaT(deltaT);
}

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
