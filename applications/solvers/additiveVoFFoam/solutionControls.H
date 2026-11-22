//- Set solution controls for the phase-change iteration
const label nThermoCorr =
    pimple.dict().lookupOrDefault<label>("nThermoCorrectors", 20);

const scalar thermoTol =
    pimple.dict().lookupOrDefault<scalar>("thermoTolerance", 1e-6);
