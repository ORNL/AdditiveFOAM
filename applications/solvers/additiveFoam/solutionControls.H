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
