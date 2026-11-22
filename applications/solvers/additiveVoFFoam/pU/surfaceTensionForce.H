surfaceScalarField surfaceTensionForce
(
    IOobject("surfaceTensionForce", runTime.name(), mesh),
    mesh,
    dimensionedScalar(dimForce/dimVolume, 0)
);

{
    // Smoothed metal fraction for the interface normal
    volScalarField alpha1a("alpha1a", alpha1);
    for (label i = 0; i < 3; i++)
    {
        alpha1a = 0.5*fvc::average(alpha1a) + 0.5*alpha1a;
        alpha1a.correctBoundaryConditions();
    }
    const volScalarField alpha2a(1 - alpha1a);

    const surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2a)*fvc::interpolate(fvc::grad(alpha1a, "nHat"))
      - fvc::interpolate(alpha1a)*fvc::interpolate(fvc::grad(alpha2a, "nHat"))
    );

    const surfaceVectorField nHatf(gradAlphaf/(mag(gradAlphaf) + deltaN));

    // Curvature, smoothed toward the interface center
    const volScalarField K(-fvc::div(nHatf & mesh.Sf()));

    const volScalarField w(sqrt(alpha1c*(1.0 - alpha1c) + 0.001));
    const volScalarField wInterface(2.0*sqrt(alpha1c*(1.0 - alpha1c)));

    auto smooth = [&](const volScalarField& K)
    {
        volScalarField Ks
        (
            wInterface*K
          + (1.0 - wInterface)*fvc::average(w*K)/fvc::average(w)
        );
        Ks.correctBoundaryConditions();
        return Ks;
    };

    const surfaceScalarField Kf
    (
        fvc::interpolate(w*smooth(smooth(K)))/fvc::interpolate(w)
    );

    const volScalarField sigma
    (
        "sigma",
        sigma0 + dSigmadT*(max(T, solidus) - liquidus)
    );

    const volScalarField alphaLiquid(max(alpha1c - alphaSolid, scalar(0)));
    const volScalarField alpha2c(1.0 - alpha1c);

    // A property of the liquid surface is evaluated on the liquid side of the
    // face. In a smeared thermal field the gas cell next to the surface is
    // colder than the surface, and a plain face interpolation would average
    // the surface tension, its gradient and the recoil pressure with values
    // that belong to no surface: the gas side carries sigma(solidus) and a
    // zero gradient below the solidus, and a recoil pressure of nothing
    const surfaceScalarField alphaLiquidf(fvc::interpolate(alphaLiquid));

    auto liquidFace = [&](const auto& field)
    {
        return fvc::interpolate(alphaLiquid*field)/max(alphaLiquidf, small);
    };

    // Normal pressure jump of the liquid-gas pair: capillary pressure and
    // recoil pressure of the evaporating liquid
    surfaceScalarField pJumpf(liquidFace(sigma)*Kf);

    if (recoil)
    {
        // Saturation pressure at the temperature of the liquid surface
        const volScalarField pRecoil
        (
            dimensionedScalar(dimPressure, pRecoil0)
           *exp
            (
                dimensionedScalar(dimTemperature, TvapTheta)
               *(
                    1.0/dimensionedScalar(dimTemperature, Tvap)
                  - 1.0/max(T, dimensionedScalar(dimTemperature, 1))
                )
            )
        );

        pJumpf += liquidFace(pRecoil);
    }

    surfaceTensionForce +=
        pJumpf
       *(
            fvc::interpolate(alpha2c)*fvc::snGrad(alphaLiquid)
          - fvc::interpolate(alphaLiquid)*fvc::snGrad(alpha2c)
        );

    // Marangoni stress, from the surface-tension gradient of the liquid side
    const surfaceVectorField gradSigmaf(liquidFace(fvc::grad(sigma)));

    surfaceTensionForce +=
        (
            (gradSigmaf - (gradSigmaf & nHatf)*nHatf) & mesh.Sf()/mesh.magSf()
        )
       *fvc::interpolate(alpha1c)
       *(
            mag(fvc::interpolate(fvc::grad(alpha1c)))
          + mag(fvc::interpolate(fvc::grad(alpha2c)))
        );
}
