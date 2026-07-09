while (pimple.correct())
{
    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf)
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    rhok = 1.0 - beta*(T - liquidus);

    surfaceScalarField phig
    (
        "phig",
        - ghf*fvc::snGrad(rhok)*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
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

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rhok*gh;

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
}
