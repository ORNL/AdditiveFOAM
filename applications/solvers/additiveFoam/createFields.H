Info<< "Reading transportProperties\n" << endl;
#include "readTransportProperties.H"

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
volScalarField alphaPowder
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

volScalarField alphaSolid
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

//- Solid fraction on the phase-change path at the initial temperature
forAll(alphaSolid, celli)
{
    alphaSolid[celli] = thermo.fs(T[celli]);
}
alphaSolid.correctBoundaryConditions();

//- Phase-change iteration of the energy equation on the thermoPath
phaseChangeSolver phaseChange(T, alphaSolid, path, rho*Lf);

//- Heat capacity, sensible enthalpy and conductivity of the cells
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTemperature, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rhoH
(
    IOobject
    (
        "rhoH",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimLength/dimTemperature, 0),
    zeroGradientFvPatchScalarField::typeName
);

//- Sensible enthalpy law of a cell: the solid-liquid mixture and, below the
//  solidus, the powder it may still be
const sensibleEnergy sensible
{
    alphaPowder.primitiveField(),
    rhoCp.primitiveField(),
    rhoH.primitiveField(),
    thermo,
    CpPowder,
    rho.value()
};

#include "thermo/updateProperties.H"
rhoCp.oldTime();
rhoH.oldTime();

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
