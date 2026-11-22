IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- Phase-change path and solid-liquid thermodynamics of the material
thermoPath path(transportProperties, runTime);

const solidLiquidThermo thermo
(
    path,
    transportProperties.subDict("solid"),
    transportProperties.subDict("liquid")
);

const dimensionedScalar solidus("solidus", dimTemperature, thermo.solidus());
const dimensionedScalar liquidus("liquidus", dimTemperature, thermo.liquidus());

//- Powder, below the solidus
const dictionary& powderDict(transportProperties.subDict("powder"));

const Polynomial<3> CpPowder(powderDict.lookup("Cp"));
const Polynomial<3> kappaPowder(powderDict.lookup("kappa"));

// Reference density [kg/m^3]
const dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);

// Thermal expansion coefficient [1/K]
const dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    transportProperties.lookup("beta")
);

//- Mushy-zone constant per unit density; a large value keeps the solid rigid
const dimensionedScalar Cu
(
    "Cu",
    dimless/dimTime,
    transportProperties.lookupOrDefault<scalar>("Cu", 1e10)
);

//- Latent heat of fusion [J/kg]
const dimensionedScalar Lf
(
    "Lf",
    dimEnergy/dimMass,
    transportProperties.lookup("Lf")
);

//- Liquid dynamic viscosity [kg/(m s)]
const dimensionedScalar mu
(
    "mu",
    dimDynamicViscosity,
    transportProperties.lookup("mu")
);

//- Liquid kinematic viscosity [m^2/s]
const dimensionedScalar nu
(
    "nu",
    dimDynamicViscosity/dimDensity,
    mu.value() / rho.value()
);
