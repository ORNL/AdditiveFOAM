Info<< "Reading transportProperties\n" << endl;
#include "readTransportProperties.H"

Info<< "Reading thermoPath\n" << endl;
thermoPath thermo(transportProperties, runTime);

dimensionedScalar liquidus
(
    "liquidus",
    dimTemperature,
    thermo.liquidus()
);

dimensionedScalar solidus
(
    "solidus",
    dimTemperature,
    thermo.solidus()
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.name(),
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
        runTime.name(),
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
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha.solid", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// set solid fraction field consistent with temperature
forAll(mesh.cells(), celli)
{
    scalar alpha1_ = interpolateXY(T[celli], thermo.x(), thermo.y());
    alpha1[celli] = min(max(alpha1_, 0.0), 1.0);
}

alpha1.correctBoundaryConditions();

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.name(),
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
        runTime.name(),
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
        runTime.name(),
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
        runTime.name(),
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
        runTime.name(),
        mesh
    ),
    1.0 - beta*(T - liquidus)
);

#include "readGravitationalAcceleration.H"

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, 0)
);

dimensionedScalar ghRef(-mag(g)*hRef);
volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rgh + rhok*gh
);

// Create pressure reference
Foam::pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict()
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );

    p_rgh = p - rhok*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

#include "createMRF.H"

autoPtr<surfaceVectorField> Uf;
if (mesh.dynamic() || MRF.size())
{
    U.correctBoundaryConditions();

    Uf = new surfaceVectorField
    (
        IOobject
        (
            "Uf",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U)
    );
}
