//- Set solution controls for thermodynamic algorithm
scalar nThermoCorr = 
    pimple.dict().lookupOrDefault<scalar>("nThermoCorrectors", 20);

scalar thermoTol = 
    pimple.dict().lookupOrDefault<scalar>("thermoTolerance", 1e-6);

// Integrate the total power input to domain
const scalar totalPower = fvc::domainIntegrate(sources.qDot()).value();

Info<< "absorbed power: " << totalPower << endl;

//- prevent iteration of PIMPLE and THERMO where needed
bool beamOn(totalPower > small);

bool fluidInDomain((1 - gMin(alpha1)) > small);

if (!(beamOn || fluidInDomain))
{
    nThermoCorr = 1;
}

//- solve energy equation explicitly whenever possible
Switch explicitSolve =
    pimple.dict().lookupOrDefault<Switch>("explicitSolve", false);

if (explicitSolve)
{
    scalar maxExplicitTimeStep
    (
        Foam::pow(0.5, mesh.nSolutionD())
      * gMin(Foam::pow(mesh.V().field(), 2.0/3.0) / (kappa/rho/Cp)())
    );

    if (runTime.deltaTValue() < maxExplicitTimeStep)
    {
        explicitSolve = true;
    }
    else
    {
        explicitSolve = false;
    }
}
