// Capillary wave limit (Brackbill et al. 1992) while a liquid-gas interface
// exists or is about to form
{
    const scalar maxCapillaryCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCapillaryCo", 1.0);

    // Surface tension is largest at the solidus for a negative dSigmadT
    const scalar sigmaMax =
        max
        (
            sigma0.value(),
            sigma0.value() + dSigmadT.value()*(solidus - liquidus).value()
        );

    const scalar dxMin = 1.0/gMax(mesh.deltaCoeffs().primitiveField());

    const scalar deltaTSigma =
        Foam::sqrt
        (
            (rho1.value() + rho2.value())*pow3(dxMin)
           /(4*constant::mathematical::pi*max(sigmaMax, small))
        );

    bool capillaryActive = false;
    {
        const surfaceScalarField alpha1f(fvc::interpolate(alpha1));
        const surfaceScalarField alphaLf
        (
            fvc::interpolate(max(alpha1 - alphaS, scalar(0)))
        );
        const surfaceScalarField Tf(fvc::interpolate(T));
        const surfaceScalarField dTf
        (
            fvc::interpolate
            (
                max(T - T.oldTime(), dimensionedScalar(dimTemperature, 0))
            )
        );
        const scalar steps = 2.0*runTime.deltaTValue()/runTime.deltaT0Value();

        forAll(alpha1f, facei)
        {
            const scalar interface = min(alpha1f[facei], 1 - alpha1f[facei]);

            if
            (
                min(alphaLf[facei], 1 - alpha1f[facei]) > 1e-3
             || (
                    interface > 1e-3
                 && Tf[facei] + steps*dTf[facei] >= solidus.value()
                )
            )
            {
                capillaryActive = true;
                break;
            }
        }
        reduce(capillaryActive, orOp<bool>());
    }

    if (capillaryActive)
    {
        deltaT = min(deltaT, maxCapillaryCo*deltaTSigma);
    }
}
