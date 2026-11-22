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

//- Metal: phase-change path and solid-liquid thermodynamics
thermoPath path(transportProperties, runTime);

const solidLiquidThermo thermo
(
    path,
    transportProperties.subDict("solid"),
    transportProperties.subDict("liquid")
);

const dimensionedScalar solidus("solidus", dimTemperature, thermo.solidus());
const dimensionedScalar liquidus("liquidus", dimTemperature, thermo.liquidus());

const dimensionedScalar rho1
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);

const dimensionedScalar mu1
(
    "mu",
    dimDynamicViscosity,
    transportProperties.lookup("mu")
);

const dimensionedScalar beta1
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

const dimensionedScalar Lf
(
    "Lf",
    dimEnergy/dimMass,
    transportProperties.lookup("Lf")
);

//- Surface tension at the liquidus and its temperature coefficient
const dimensionedScalar sigma0
(
    "sigma",
    dimForce/dimLength,
    transportProperties.lookup("sigma")
);

const dimensionedScalar dSigmadT
(
    "dSigmadT",
    dimForce/dimLength/dimTemperature,
    transportProperties.lookup<scalar>("dSigmadT")
);

//- Optional recoil pressure, 0.54*p0*exp(Lv*M/R*(1/Tvap - 1/T))
const bool recoil(transportProperties.found("vaporization"));
const dictionary vaporizationDict
(
    transportProperties.subOrEmptyDict("vaporization")
);
const scalar Tvap(recoil ? vaporizationDict.lookup<scalar>("Tvap") : 1);
const scalar TvapTheta
(
    recoil
  ? vaporizationDict.lookup<scalar>("Lv")
   *vaporizationDict.lookup<scalar>("M")
   /constant::physicoChemical::R.value()
  : 0
);
const scalar pRecoil0
(
    0.54*vaporizationDict.lookupOrDefault<scalar>("p0", 101325)
);

//- Gas
const dictionary& gasDict(transportProperties.subDict("gas"));

const dimensionedScalar rho2("rho", dimDensity, gasDict.lookup("rho"));
const dimensionedScalar mu2("mu", dimDynamicViscosity, gasDict.lookup("mu"));
const dimensionedScalar beta2
(
    "beta",
    dimless/dimTemperature,
    gasDict.lookupOrDefault<scalar>("beta", 0)
);
const Polynomial<3> CpGas(gasDict.lookup("Cp"));
const Polynomial<3> kappaGas(gasDict.lookup("kappa"));
