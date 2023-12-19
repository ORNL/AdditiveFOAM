if (adjustTimeStep)
{
    scalar maxAlphaCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", 1.0);

    const volScalarField::Internal alphaNum
    (
        mag(fvc::ddt(T))()()
       /(Tliq - Tsol)
      *runTime.deltaT()
    );

    alphaCoNum = gMax(alphaNum);

    Info << "Thermo Number max: " << alphaCoNum << endl;

    scalar maxDi = runTime.controlDict().lookupOrDefault<scalar>("maxDi", 10.0);

    const volScalarField::Internal DiNumvf
    (
        fvc::surfaceSum
        (
            mesh.magSf()
           *fvc::interpolate(kappa)
           *mesh.surfaceInterpolation::deltaCoeffs()
        )()()
       /(mesh.V()*rho.value()*Cp())
       *runTime.deltaT()
    );

    DiNum = gMax(DiNumvf);

    Info<< "Diffusion Number max: " << DiNum << endl;

    scalar maxDeltaTFact = min
    (
        min(maxCo/(CoNum + small), maxAlphaCo/(alphaCoNum + small)),
        maxDi/(DiNum + small)
    );

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    scalar deltaT = min(deltaTFact*runTime.deltaTValue(), maxDeltaT);

    sources.adjustDeltaT(deltaT);

    runTime.setDeltaT(deltaT);
}

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
