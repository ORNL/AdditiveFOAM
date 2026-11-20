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

    const scalar L = Lf.value();

    //- solid fraction before each solve, for the linearization and the
    //  residual; allocated once
    volScalarField alpha10("alpha10", alpha1);

    //- the implicit temperature limit is assembled only when one is set
    const bool limitT = Tmax.value() < vGreat;

    for (int tCorr=0; tCorr < nThermoCorr; tCorr++)
    {
        forAll(T, celli)
        {
            T0[celli] = T[celli];
            dFdT[celli] = thermo.dAlphadT(T[celli]);
        }
        T0.correctBoundaryConditions();
        dFdT.correctBoundaryConditions();

        alpha10 = alpha1;

        //- the energy equation with the linearized latent heat
        tmp<fvScalarMatrix> tTEqnLatent
        (
            TEqn
         ==
            (
                (ddtScheme == "Euler")
              ? rho*Lf*fv::EulerDdtScheme<scalar>(mesh).fvcDdt(alpha1)
              : rho*Lf*fvc::ddt(alpha1)
            )
          + rDeltaT*rho*Lf*(fvm::Sp(dFdT, T) - dFdT*T0)
        );

        //- optional implicit limiting of temperature field
        if (limitT)
        {
            const volScalarField A
            (
                "A",
                dimensionedScalar(rho.dimensions()*Cp.dimensions(), 1e15)
               *pos(T - Tmax)
            );

            tTEqnLatent.ref() += rDeltaT*(fvm::Sp(A, T) - A*Tmax);
        }

        tTEqnLatent.ref().solve();

        //- put each cell on the path at the energy of the linearized solution,
        //  Cp*(t - T*) = Lf*(alpha(t) - alpha*), whose left side minus right
        //  side increases with t

        scalar residual = 0;

        forAll(T, celli)
        {
            const scalar Tstar = T[celli];
            const scalar c = Cp[celli]/L;
            const scalar alphaStar =
                alpha10[celli] + dFdT[celli]*(Tstar - T0[celli]);

            //- a cell whose linearized state is on the path is its own
            //  projection: every cell outside the melting range that stays
            //  there, which is most of the mesh on every solve
            if (thermo.alpha(Tstar) == alphaStar)
            {
                alpha1[celli] = alphaStar;
                residual = max(residual, mag(alphaStar - alpha10[celli]));
                continue;
            }

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

            //- safeguarded Newton iteration: a Newton step inside the bracket,
            //  otherwise bisection. It ends when the update no longer moves
            //  the temperature or the bracket has closed
            scalar t = Tstar;
            const scalar tolT = 1e-12*max(mag(t), scalar(1));
            for (label it = 0; it < 60; it++)
            {
                const scalar g = f(t);
                if (g > 0) { hi = t; } else { lo = t; }
                const scalar tn = t - g/(c - thermo.dAlphadT(t));
                const scalar tNew = (tn > lo && tn < hi) ? tn : 0.5*(lo + hi);
                const bool done = mag(tNew - t) <= tolT || hi - lo <= tolT;
                t = tNew;
                if (done) { break; }
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
            residual = max(residual, mag(alpha1[celli] - alpha10[celli]));
        }

        T.correctBoundaryConditions();
        alpha1.correctBoundaryConditions();

        //- check convergence: largest change of the solid fraction in a cell
        reduce(residual, maxOp<scalar>());

        Info<< "Thermo: iteration " << tCorr
            << " residual: " << residual << endl;

        if ((residual < thermoTol) && (tCorr > 0))
        {
            break;
        }
    }
}
