{
    rAU = 1.0/UEqn.A();

    // Face mobility of the interpolated momentum diagonal
    const surfaceScalarField rAUf("rAUf", 1.0/fvc::interpolate(UEqn.A()));
    const surfaceScalarField rAUDf("rAUDf", rAUf/(1.0 + rAUf*Df));

    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));

    surfaceScalarField phiH("phiH", rAUf*fvc::flux(UEqn.H()));
    {
        const surfaceScalarField phiHbyAc(fvc::flux(HbyA));

        forAll(phiH.boundaryField(), patchi)
        {
            if (!phiH.boundaryField()[patchi].coupled())
            {
                phiH.boundaryFieldRef()[patchi] =
                    phiHbyAc.boundaryField()[patchi];
            }
        }
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            phiH
          + rAUf*fvc::interpolate(rho)*fvc::ddtCorr(U, phi, Uf)
        )/(1.0 + rAUf*Df)
    );

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    const surfaceScalarField phig
    (
        (surfaceTensionForce - ghf*fvc::snGrad(rhoBuoyant))*rAUDf*mesh.magSf()
    );

    phiHbyA += phig;

    constrainPressure(p_rgh, U, phiHbyA, rAUDf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUDf, p_rgh)
        );

        p_rghEqn.setReference
        (
            pressureReference.refCell(),
            getRefCellValue(p_rgh, pressureReference.refCell())
        );

        p_rghEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            // Velocity from the same face-weighted acceleration as phi
            U =
            (
                HbyA
              + rAU()*D*fvc::reconstruct(phi)
              + fvc::reconstruct(rAUf*((phig + p_rghEqn.flux())/rAUDf - Df*phi))
            )/(1 + rAU()*D);

            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rhoBuoyant*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureReference.refValue()
          - getRefCellValue(p, pressureReference.refCell())
        );
        p_rgh = p - rhoBuoyant*gh;
    }
}
