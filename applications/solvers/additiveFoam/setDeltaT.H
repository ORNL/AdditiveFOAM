if (adjustTimeStep)
{
    // temporal resolution w.r.t freezing range
    scalar maxAlphaCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", 0.1);

    alphaCoNum = gMax(mag(T - T.oldTime())()) / (Tliq - Tsol).value();

    Info << "alphaCoNum: " << alphaCoNum << endl;

    scalar maxFo = runTime.controlDict().lookupOrDefault<scalar>("maxFo", 1);

    const scalar diffusionTime
    (
        gMin(Foam::pow(mesh.V().field(), 2.0/3.0) / (kappa/rho/Cp)())
    );

    FoNum = runTime.deltaTValue() / diffusionTime;

    Info << "FoNum: " << FoNum << endl;

    scalar deltaT = min
    (
        min(maxCo/(CoNum + small), maxAlphaCo/(alphaCoNum + small)),
        maxFo/(FoNum + small)
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
