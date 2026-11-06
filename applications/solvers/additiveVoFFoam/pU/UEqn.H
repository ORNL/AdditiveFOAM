#include "pU/surfaceTensionForce.H"

// Mushy-zone drag of the solid fraction (Kozeny-Carman)
const volScalarField alphaSc(min(alphaS, scalar(0.95)));

const volScalarField D
(
    "D",
    mu1*(180.0/sqr(DAS))*sqr(alphaSc)/pow3(1.0 - alphaSc)
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
