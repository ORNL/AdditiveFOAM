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

const int PolySize(3);
const dictionary& solidDict(transportProperties.subDict("solid"));
const dictionary& liquidDict(transportProperties.subDict("liquid"));
const dictionary& powderDict(transportProperties.subDict("powder"));


Polynomial<PolySize> kappa1(solidDict.lookup("kappa"));
Polynomial<PolySize> kappa2(liquidDict.lookup("kappa"));
Polynomial<PolySize> kappa3(powderDict.lookup("kappa"));

Polynomial<PolySize> Cp1(solidDict.lookup("Cp"));
Polynomial<PolySize> Cp2(liquidDict.lookup("Cp"));
Polynomial<PolySize> Cp3(powderDict.lookup("Cp"));

// Reference density [kg]
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

//- Length scale of mushy zone [m]
const dimensionedScalar DAS
(
    "DAS",
    dimLength,
    transportProperties.lookup("DAS")
);

//- Latent heat of fusion [J/kg]
const dimensionedScalar Lf
(
    "Lf", 
    dimEnergy/dimMass, 
    transportProperties.lookup("Lf")
);

//- Liquid dynamic viscosity [m^2/(kg s)]
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
