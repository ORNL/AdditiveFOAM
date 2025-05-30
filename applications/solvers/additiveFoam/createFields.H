Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Setting auxiliary fields for solidification model\n" << endl;
volScalarField alpha3
(
    IOobject
    (
        IOobject::groupName("alpha", "powder"),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha.powder", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alpha1
(
    IOobject
    (
        IOobject::groupName("alpha", "solid"),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha.solid", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

IFstream thermoFile(runTime.constant()/"thermoPath");

graph thermo
(
    "thermo",
    "T",
    "alpha1",
    thermoFile
);

dimensionedScalar Tliq
(
    "Tliq",
    dimTemperature,
    interpolateXY(0.0, thermo.y(), thermo.x())
);

dimensionedScalar Tsol
(
    "Tsol",
    dimTemperature,
    interpolateXY(1.0, thermo.y(), thermo.x())
);

// set solid fraction field consistent with temperature
forAll(mesh.cells(), cellI)
{
    scalar alpha1_ = interpolateXY(T[cellI], thermo.x(), thermo.y());
    alpha1[cellI] = min(max(alpha1_, 0.0), 1.0);
}

alpha1.correctBoundaryConditions();

#include "readTransportProperties.H"

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cp", dimensionSet(0, 2, -2, -1, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("kappa", dimensionSet(1, 1, -3, -1, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

#include "updateProperties.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p_rgh", dimensionSet(0, 2, -2, 0, 0), 0)
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, Zero)
);

#include "createPhi.H"

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - Tliq)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.schemes().setFluxRequired(p_rgh.name());

#include "createMRF.H"
