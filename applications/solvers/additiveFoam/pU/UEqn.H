// Solve the Momentum equation

// Mushy-zone drag of the solid volume fraction (Carman-Kozeny), the same
// coefficient as additiveVoFFoam per unit density
const volScalarField porosity(min(max(1.0 - alphaSolid, scalar(0)), scalar(1)));

const volScalarField D
(
    "D",
    min
    (
        Cu*sqr(1.0 - porosity)/(pow3(porosity) + 0.001),
        dimensionedScalar(dimless/dimTime, 1e12)
    )
);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  - fvm::laplacian(nu, U)
 ==
  - fvm::Sp(D, U)
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );
}
