#include "pU/surfaceTensionForce.H"

// Mushy-zone drag of the solid volume fraction (Carman-Kozeny)
const volScalarField porosity(min(max(1.0 - alphaS, scalar(0)), scalar(1)));

const volScalarField D
(
    "D",
    rho
   *min
    (
        Cu*sqr(1.0 - porosity)/(pow3(porosity) + 0.001),
        dimensionedScalar(dimless/dimTime, 1e12)
    )
);

const surfaceScalarField Df("Df", fvc::interpolate(D));

// Buoyancy density: the fluid weight with thermal expansion; the solid
// carries its own weight
const volScalarField alphaL(alpha1 - min(alphaS, alpha1));

const volScalarField rhoBuoyant
(
    "rhoBuoyant",
    rho - min(alphaS, alpha1)*(rho1 - rho2)
  - alphaL*rho1*beta1*(T - liquidus)
  - (1.0 - alphaL)*rho2*beta2*(T - liquidus)
);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  - fvm::laplacian(mu, U)
  - fvc::div(mu*dev2(Foam::T(fvc::grad(U))))
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn + fvm::Sp(D, U)
     ==
        UEqn.A()*fvc::reconstruct
        (
            (
                (
                    surfaceTensionForce
                  - ghf*fvc::snGrad(rhoBuoyant)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
              - Df*phi
            )/fvc::interpolate(UEqn.A())
        )
      + D*fvc::reconstruct(phi)
    );
}
