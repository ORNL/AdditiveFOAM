{
    #include "thermo/thermoScheme.H"

    //- solve thermodynamic coupling with the optimum approach: a Newton
    //  linearization of the solid fraction, then projection of each cell onto
    //  the thermodynamic path at the energy of the linearized solution
    volScalarField dFdT
    (
        IOobject
        (
            "dFdT",
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedScalar("dFdT", dimless/dimTemperature, 0.0)
    );

    volScalarField T0("T0", T);

    //- store the old-time fields before the cell updates below modify them
    T.oldTime();
    alpha1.oldTime();

    const scalarField& V = mesh.V();
    const scalar L = Lf.value();

    for (int tCorr=0; tCorr < nThermoCorr; tCorr++)
    {
        forAll(T, celli)
        {
            T0[celli] = T[celli];
            dFdT[celli] = thermo.dAlphadT(T[celli]);
        }
        T0.correctBoundaryConditions();
        dFdT.correctBoundaryConditions();

        const volScalarField alpha10("alpha10", alpha1);

        //- optional implicit limiting of temperature field
        volScalarField A
        (
            "A",
            dimensionedScalar(rho.dimensions()*Cp.dimensions(), 1e15)
           *pos(T - Tmax)
        );

        //- solve the energy equation
        solve
        (
            TEqn + rDeltaT*(fvm::Sp(A, T) - A*Tmax)
         ==
            (
                (ddtScheme == "Euler")
              ? rho*Lf*fv::EulerDdtScheme<scalar>(mesh).fvcDdt(alpha1)
              : rho*Lf*fvc::ddt(alpha1)
            )
          + rDeltaT*rho*Lf*(fvm::Sp(dFdT, T) - dFdT*T0)
        );

        //- put each cell on the path at the energy of the linearized solution,
        //  Cp*(t - T*) = Lf*(alpha(t) - alpha*), whose left side minus right
        //  side increases with t
        scalar energyCorrection = 0;

        forAll(T, celli)
        {
            const scalar Tstar = T[celli];
            const scalar c = Cp[celli]/L;
            const scalar alphaStar =
                alpha10[celli] + dFdT[celli]*(Tstar - T0[celli]);

            auto f = [&](const scalar t)
            {
                return c*(t - Tstar) - (thermo.alpha(t) - alphaStar);
            };

            scalar lo = min(T0[celli], Tstar);
            scalar hi = max(T0[celli], Tstar);
            scalar span = max(hi - lo, scalar(1));
            while (f(lo) > 0) { lo -= span; span *= 2; }
            span = max(hi - lo, scalar(1));
            while (f(hi) < 0) { hi += span; span *= 2; }

            scalar t = Tstar;
            for (label it = 0; it < 60; it++)
            {
                const scalar g = f(t);
                if (g > 0) { hi = t; } else { lo = t; }
                if (hi - lo < 1e-9*max(t, scalar(1))) { break; }
                const scalar tn = t - g/(c - thermo.dAlphadT(t));
                t = (tn > lo && tn < hi) ? tn : 0.5*(lo + hi);
            }

            //- on an isothermal step of the path the solid fraction follows
            //  from the energy, within the step
            const scalar eps = 1e-6*max(t, scalar(1));
            T[celli] = t;
            alpha1[celli] = min
            (
                max(alphaStar + c*(t - Tstar), thermo.alpha(t + eps)),
                thermo.alpha(t - eps)
            );

            energyCorrection +=
                V[celli]
               *(
                    Cp[celli]*mag(t - T0[celli])
                  + L*mag(alpha1[celli] - alpha10[celli])
                );
        }

        T.correctBoundaryConditions();
        alpha1.correctBoundaryConditions();

        //- check convergence: energy moved by this iteration relative to the
        //  energy exchanged in the step
        const scalar energyStep = gSum
        (
            V
           *(
                Cp.primitiveField()
               *mag(T.primitiveField() - T.oldTime().primitiveField())
              + L
               *mag
                (
                    alpha1.primitiveField() - alpha1.oldTime().primitiveField()
                )
              + mag(sources.qDot().primitiveField())
               /(rho.value()*rDeltaT.value())
            )
        );

        const scalar residual =
            returnReduce(energyCorrection, sumOp<scalar>())
           /max(energyStep, vSmall);

        Info<< "Thermo: iteration " << tCorr
            << " residual: " << residual << endl;

        if ((residual < thermoTol) && (tCorr > 0))
        {
            break;
        }
    }
}
