if (adjustTimeStep)
{
    // Phase-change Courant number: the change of the solid fraction a cell
    // can make in one step on a linear path through the melting range
    const scalar maxPhaseChangeCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxPhaseChangeCo", 1.0);

    const volScalarField::Internal phaseChangeNum
    (
        mag(fvc::ddt(T)) / (liquidus - solidus)
       *runTime.deltaT()
    );

    phaseChangeCoNum = gMax(phaseChangeNum);

    Info<< "Phase-change Courant Number max: " << phaseChangeCoNum << endl;

    scalar maxDi = runTime.controlDict().lookupOrDefault<scalar>("maxDi", 10.0);

    const volScalarField::Internal DiNumvf
    (
        fvc::surfaceSum
        (
            mesh.magSf()
           *fvc::interpolate(kappa)
           *mesh.surfaceInterpolation::deltaCoeffs()
        )
       /(mesh.V()*rhoCp)
       *runTime.deltaT()
    );

    DiNum = gMax(DiNumvf);

    Info<< "Diffusion Number max: " << DiNum << endl;

    scalar maxDeltaTFact = min
    (
        min(maxCo/(CoNum + small), maxPhaseChangeCo/(phaseChangeCoNum + small)),
        maxDi/(DiNum + small)
    );

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    scalar deltaT = min(deltaTFact*runTime.deltaTValue(), maxDeltaT);

    // The write-time adjustment of setDeltaT divides the time to the next
    // write evenly and replaces the time step, so the scan-path adjustment
    // is applied to its result
    runTime.setDeltaT(deltaT);

    deltaT = runTime.deltaTValue();
    sources.adjustDeltaT(deltaT);
    runTime.setDeltaTNoAdjust(deltaT);
}

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
