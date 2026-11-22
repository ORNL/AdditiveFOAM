//- Phase properties and sensible enthalpy at the current temperature
{
    const scalar r1 = rho1.value();
    const scalar r2 = rho2.value();
    const scalar Ts = solidus.value();

    //- Solid fraction of the metal from the state of the cell: on an
    //  isothermal step of the path it is not a function of the temperature
    auto set = [&]
    (
        const scalar a,
        const scalar aSolid,
        const scalar t,
        scalar& cp1,
        scalar& cp2,
        scalar& hh1,
        scalar& hh2,
        scalar& rcp,
        scalar& k
    )
    {
        // The thermal properties are evaluated with the clipped metal
        // fraction: a transport undershoot must not give a cell a negative
        // heat capacity or conductivity (the mass flux keeps the raw fraction)
        const scalar ac = min(max(a, scalar(0)), scalar(1));
        const scalar f =
            ac > small
          ? min(max(aSolid/ac, scalar(0)), scalar(1))
          : thermo.fs(t);

        cp1 = thermo.Cp(t, f);
        cp2 = CpGas.value(t);
        hh1 = thermo.h(t);
        hh2 = CpGas.integral(Ts, t);
        rcp = ac*r1*cp1 + (1 - ac)*r2*cp2;
        k = ac*thermo.kappa(t, f) + (1 - ac)*kappaGas.value(t);
    };

    // primitiveFieldRef() stores the old-time fields before they are
    // overwritten
    scalarField& Cp1Cells = Cp1.primitiveFieldRef();
    scalarField& Cp2Cells = Cp2.primitiveFieldRef();
    scalarField& h1Cells = h1.primitiveFieldRef();
    scalarField& h2Cells = h2.primitiveFieldRef();
    scalarField& rhoCpCells = rhoCp.primitiveFieldRef();
    scalarField& kappaCells = kappa.primitiveFieldRef();

    forAll(T, celli)
    {
        set
        (
            alpha1[celli], alphaSolid[celli], T[celli],
            Cp1Cells[celli], Cp2Cells[celli], h1Cells[celli], h2Cells[celli],
            rhoCpCells[celli], kappaCells[celli]
        );
    }

    forAll(T.boundaryField(), patchi)
    {
        const scalarField& Tp = T.boundaryField()[patchi];
        const scalarField& alpha1p = alpha1.boundaryField()[patchi];
        const scalarField& alphaSolidp = alphaSolid.boundaryField()[patchi];

        forAll(Tp, facei)
        {
            set
            (
                alpha1p[facei], alphaSolidp[facei], Tp[facei],
                Cp1.boundaryFieldRef()[patchi][facei],
                Cp2.boundaryFieldRef()[patchi][facei],
                h1.boundaryFieldRef()[patchi][facei],
                h2.boundaryFieldRef()[patchi][facei],
                rhoCp.boundaryFieldRef()[patchi][facei],
                kappa.boundaryFieldRef()[patchi][facei]
            );
        }
    }

    rhoH = alpha1c*rho1*h1 + (1 - alpha1c)*rho2*h2;
    rhoHPhi =
        rho1*fvc::flux(alphaPhi1, h1, "div(rhoPhi,T)")
      + rho2*fvc::flux(phi - alphaPhi1, h2, "div(rhoPhi,T)");
}
