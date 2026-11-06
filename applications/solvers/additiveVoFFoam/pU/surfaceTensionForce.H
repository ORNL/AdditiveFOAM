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

    // Curvature, smoothed towards the interface centre
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

    // Normal pressure jump of the liquid-gas pair: capillary pressure and
    // recoil pressure of the evaporating liquid
    surfaceScalarField pJumpf(fvc::interpolate(sigma)*Kf);

    if (recoil)
    {
        pJumpf += fvc::interpolate
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
    }

    const volScalarField alphaL(max(alpha1c - alphaS, scalar(0)));
    const volScalarField alpha2c(1.0 - alpha1c);

    surfaceTensionForce +=
        pJumpf
       *(
            fvc::interpolate(alpha2c)*fvc::snGrad(alphaL)
          - fvc::interpolate(alphaL)*fvc::snGrad(alpha2c)
        );

    // Marangoni stress
    const surfaceVectorField gradSigmaf(fvc::interpolate(fvc::grad(sigma)));

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
