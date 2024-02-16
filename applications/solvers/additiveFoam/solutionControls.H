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

if (explicitSolve && !(sources.refinementControl().update()))
{
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

    explicitSolve = (gMax(DiNumvf) < 1.0) ? true : false;
}
