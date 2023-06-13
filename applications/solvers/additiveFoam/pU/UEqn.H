// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

const volScalarField D
(
    "D",
    nu*(180.0/sqr(DAS))*pow(min(alpha1, 0.95), 2.0)
  / pow(1.0 - min(alpha1, 0.95), 3.0)
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
