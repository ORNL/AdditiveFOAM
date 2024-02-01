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

Polynomial<PolySize> kappa1(solidDict.get<Polynomial<PolySize>>("kappa"));
Polynomial<PolySize> kappa2(liquidDict.get<Polynomial<PolySize>>("kappa"));
Polynomial<PolySize> kappa3(powderDict.get<Polynomial<PolySize>>("kappa"));

Polynomial<PolySize> Cp1(solidDict.get<Polynomial<PolySize>>("Cp"));
Polynomial<PolySize> Cp2(liquidDict.get<Polynomial<PolySize>>("Cp"));
Polynomial<PolySize> Cp3(powderDict.get<Polynomial<PolySize>>("Cp"));

// Reference density [kg]
const dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties.get<dimensionedScalar>("rho").value()
);

// Thermal expansion coefficient [1/K]
const dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    transportProperties.get<dimensionedScalar>("beta").value()
);

//- Length scale of mushy zone [m]
const dimensionedScalar DAS
(
    "DAS",
    dimLength,
    transportProperties.get<dimensionedScalar>("DAS").value()
);

//- Latent heat of fusion [J/kg]
const dimensionedScalar Lf
(
    "Lf", 
    dimEnergy/dimMass, 
    transportProperties.get<dimensionedScalar>("Lf").value()
);

//- Liquid dynamic viscosity [m^2/(kg s)]
const dimensionedScalar mu
(
    "mu",
    dimDynamicViscosity,
    transportProperties.get<dimensionedScalar>("mu").value()
);

//- Liquid kinematic viscosity [m^2/s]
const dimensionedScalar nu
(
    "nu",
    dimDynamicViscosity/dimDensity,
    mu.value() / rho.value()
);
