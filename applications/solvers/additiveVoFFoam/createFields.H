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

Info<< "Reading field alpha.metal\n" << endl;
volScalarField alpha1
(
    IOobject
    (
        IOobject::groupName("alpha", "metal"),
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField alpha2("alpha.gas", 1.0 - alpha1);

const dimensionedScalar deltaN
(
    "deltaN",
    1e-8/pow(average(mesh.V()), 1.0/3.0)
);

surfaceScalarField nHatf
(
    IOobject("nHatf", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimArea, 0)
);

#include "interfaceNormal.H"

//- Clipped metal fraction
volScalarField alpha1c("alpha1c", min(max(alpha1, scalar(0)), scalar(1)));

//- Solid volume fraction of the cell, on the phase-change path
volScalarField alphaSolid
(
    IOobject
    (
        IOobject::groupName("alpha", "solid"),
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1c,
    zeroGradientFvPatchScalarField::typeName
);

forAll(alphaSolid, celli)
{
    alphaSolid[celli] = alpha1c[celli]*thermo.fs(T[celli]);
}
alphaSolid.correctBoundaryConditions();

//- Phase-change iteration of the energy equation on the thermoPath, with the
//  latent heat weighted by the metal fraction of the cell
phaseChangeSolver phaseChange(T, alphaSolid, path, rho1*Lf, &alpha1c);

//- Volume fraction of the cell occupied by metal that has been liquid
volScalarField alphaMelted
(
    IOobject
    (
        IOobject::groupName("alpha", "melted"),
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

mesh.schemes().setFluxRequired(alpha1.name());

surfaceScalarField alphaPhi1
(
    IOobject
    (
        IOobject::groupName("alphaPhi", "metal"),
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.name(),
        mesh
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.name(),
        mesh
    ),
    alphaPhi1*(rho1 - rho2) + phi*rho2
);

volScalarField mu("mu", alpha1c*mu1 + (1.0 - alpha1c)*mu2);

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
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

Foam::pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

#include "createMRF.H"

autoPtr<surfaceVectorField> Uf;

tmp<volScalarField> rAU;

//- Phase heat capacities and sensible enthalpies, mixture heat capacity,
//  enthalpy, conductivity and convected enthalpy flux
volScalarField Cp1
(
    IOobject("Cp1", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimSpecificHeatCapacity, 0)
);

volScalarField Cp2
(
    IOobject("Cp2", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimSpecificHeatCapacity, 0)
);

volScalarField h1
(
    IOobject("h1", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimEnergy/dimMass, 0)
);

volScalarField h2
(
    IOobject("h2", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimEnergy/dimMass, 0)
);

volScalarField rhoCp
(
    IOobject("rhoCp", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTemperature, 0)
);

volScalarField rhoH
(
    IOobject("rhoH", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume, 0)
);

volScalarField kappa
(
    IOobject("kappa", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimPower/dimLength/dimTemperature, 0)
);

surfaceScalarField rhoHPhi
(
    IOobject("rhoHPhi", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimPower, 0)
);

//- Sensible enthalpy law of a cell: the metal on its path and the gas, on
//  the clipped metal fraction like the other thermal properties
const sensibleEnergy sensible
{
    alpha1c.primitiveField(),
    rhoCp.primitiveField(),
    rhoH.primitiveField(),
    thermo,
    CpGas,
    rho1.value(),
    rho2.value()
};

//- Surface cooling (evaporation, radiation): implicit coefficient and
//  explicit remainder of the linearized sink, and the step's powers
volScalarField surfaceCoolingSp
(
    IOobject("surfaceCoolingSp", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimPower/dimVolume/dimTemperature, 0)
);
volScalarField surfaceCoolingSu
(
    IOobject("surfaceCoolingSu", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimPower/dimVolume, 0)
);
scalar evaporationPower = 0, radiationPower = 0, surfaceArea = 0;

#include "thermo/updateProperties.H"
rhoCp.oldTime();
rhoH.oldTime();

mesh.schemes().setFluxRequired(T.name());
