{
    const dictionary& alphaControls =
        mesh.solution().solverDict(alpha1.name());

    const label nAlphaCorr =
        alphaControls.lookupOrDefault<label>("nCorrectors", 1);

    const Switch MULESCorr =
        alphaControls.lookupOrDefault<Switch>("MULESCorr", false);

    const MULES::control MULEScontrols(alphaControls);

    if (MULESCorr)
    {
        // Bounded upwind predictor
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
        );

        alpha1Eqn.solve();

        alphaPhi1 = alpha1Eqn.flux();
    }

    for (label aCorr = 0; aCorr < nAlphaCorr; aCorr++)
    {
        const surfaceScalarField alphaPhi1Un
        (
            fvc::flux(phi, alpha1, "div(phi,alpha)")
        );

        if (MULESCorr)
        {
            surfaceScalarField alphaPhi1Corr(alphaPhi1Un - alphaPhi1);
            const volScalarField alpha10("alpha10", alpha1);

            MULES::correct
            (
                MULEScontrols,
                geometricOneField(),
                alpha1,
                alphaPhi1Un,
                alphaPhi1Corr,
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi1 += alphaPhi1Corr;
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi1 += 0.5*alphaPhi1Corr;
            }
        }
        else
        {
            alphaPhi1 = alphaPhi1Un;

            MULES::explicitSolve
            (
                MULEScontrols,
                geometricOneField(),
                alpha1,
                phi,
                alphaPhi1,
                oneField(),
                zeroField()
            );
        }
    }

    alpha2 = 1.0 - alpha1;

    #include "interfaceNormal.H"

    alpha1c = min(max(alpha1, scalar(0)), scalar(1));

    rho == alpha1*rho1 + alpha2*rho2;
    rhoPhi = alphaPhi1*(rho1 - rho2) + phi*rho2;
    mu = alpha1c*mu1 + (1.0 - alpha1c)*mu2;

    Info<< "Metal volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
