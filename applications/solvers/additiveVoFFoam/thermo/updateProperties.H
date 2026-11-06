//- Phase properties and sensible enthalpy at the current temperature
{
    const scalar r1 = rho1.value();
    const scalar r2 = rho2.value();
    const scalar Ts = solidus.value();

    auto set = [&]
    (
        const scalar a,
        const scalar t,
        scalar& cp1,
        scalar& cp2,
        scalar& hh1,
        scalar& hh2,
        scalar& rcp,
        scalar& k
    )
    {
        cp1 = metal.Cp(t);
        cp2 = CpGas.value(t);
        hh1 = metal.h(t);
        hh2 = CpGas.integral(Ts, t);
        rcp = a*r1*cp1 + (1 - a)*r2*cp2;
        k = a*metal.kappa(t) + (1 - a)*kappaGas.value(t);
    };

    forAll(T, celli)
    {
        set
        (
            alpha1[celli], T[celli],
            Cp1[celli], Cp2[celli], h1[celli], h2[celli],
            rhoCp[celli], kappa[celli]
        );
    }

    forAll(T.boundaryField(), patchi)
    {
        const scalarField& Tp = T.boundaryField()[patchi];
        const scalarField& alpha1p = alpha1.boundaryField()[patchi];

        forAll(Tp, facei)
        {
            set
            (
                alpha1p[facei], Tp[facei],
                Cp1.boundaryFieldRef()[patchi][facei],
                Cp2.boundaryFieldRef()[patchi][facei],
                h1.boundaryFieldRef()[patchi][facei],
                h2.boundaryFieldRef()[patchi][facei],
                rhoCp.boundaryFieldRef()[patchi][facei],
                kappa.boundaryFieldRef()[patchi][facei]
            );
        }
    }

    rhoH = alpha1*rho1*h1 + alpha2*rho2*h2;
    rhoHPhi =
        rho1*fvc::flux(alphaPhi1, h1, "div(rhoPhi,T)")
      + rho2*fvc::flux(phi - alphaPhi1, h2, "div(rhoPhi,T)");
}
