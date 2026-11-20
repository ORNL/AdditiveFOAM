{
    const label nThermoCorr =
        pimple.dict().lookupOrDefault<label>("nThermoCorrectors", 20);
    const scalar thermoTol =
        pimple.dict().lookupOrDefault<scalar>("thermoTolerance", 1e-6);

    const dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
    const surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    const scalar r1 = rho1.value();
    const scalar r2 = rho2.value();
    const scalar L = Lf.value();
    const scalar Ts = solidus.value();

    volScalarField dFdT
    (
        IOobject("dFdT", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless/dimTemperature, 0)
    );

    volScalarField T0("T0", T);

    alphaSolid.oldTime();

    forAll(alphaSolid, celli)
    {
        alphaSolid[celli] = alpha1c[celli]*metal.fs(T[celli]);
    }
    alphaSolid.correctBoundaryConditions();

    #include "thermo/updateProperties.H"

    bool converged = false;

    for (label tCorr = 0; tCorr < nThermoCorr && !converged; tCorr++)
    {
        // Newton linearization of the solid fraction about the current T
        forAll(T, celli)
        {
            T0[celli] = T[celli];
            dFdT[celli] = alpha1c[celli]*metal.dfsdT(T[celli]);
        }
        T0.correctBoundaryConditions();
        dFdT.correctBoundaryConditions();

        const volScalarField Tprev("Tprev", T);
        const volScalarField alphaSolidPrev("alphaSolidPrev", alphaSolid);

        const surfaceScalarField rhoCpPhi
        (
            "rhoCpPhi",
            rho1*alphaPhi1*linearInterpolate(Cp1)
          + rho2*alphaPhi2*linearInterpolate(Cp2)
        );

        // Upwind in both convection terms: their high-order parts cancel
        const fv::gaussConvectionScheme<scalar> upwindCp
        (
            mesh,
            rhoCpPhi,
            tmp<surfaceInterpolationScheme<scalar>>
            (
                new upwind<scalar>(mesh, rhoCpPhi)
            )
        );

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
          - fvc::ddt(rhoCp, T) + fvc::ddt(rhoH)
          + upwindCp.fvmDiv(rhoCpPhi, T)
          - upwindCp.fvcDiv(rhoCpPhi, T)
          + fvc::div(rhoHPhi)
          - fvm::laplacian(kappa, T)
         ==
            sources.qDot()
          + rho1*Lf*fvc::ddt(alphaSolid)
          + rDeltaT*rho1*Lf*(fvm::Sp(dFdT, T) - dFdT*T0)
        );

        TEqn.relax();
        TEqn.solve();

        // Optimum approach: the cell energy follows from the solved
        // linearization and T is put back on the energy curve
        scalarField& Tc = T.primitiveFieldRef();

        forAll(Tc, i)
        {
            const scalar a = alpha1[i];
            const scalar ac = alpha1c[i];
            const scalar Tk = Tprev[i];

            const scalar Estar =
                rhoH[i] + rhoCp[i]*(Tc[i] - Tk)
              + r1*L*(ac - alphaSolidPrev[i] - dFdT[i]*(Tc[i] - Tk));

            auto energy = [&](const scalar t)
            {
                return
                    a*r1*metal.h(t) + (1 - a)*r2*CpGas.integral(Ts, t)
                  + r1*L*ac*(1 - metal.fs(t));
            };

            auto capacity = [&](const scalar t)
            {
                return
                    a*r1*metal.Cp(t) + (1 - a)*r2*CpGas.value(t)
                  - r1*L*ac*metal.dfsdT(t);
            };

            scalar lo = min(Tk, Tc[i]);
            scalar hi = max(Tk, Tc[i]);
            scalar span = max(hi - lo, scalar(1));
            while (energy(lo) > Estar) { lo -= span; span *= 2; }
            span = max(hi - lo, scalar(1));
            while (energy(hi) < Estar) { hi += span; span *= 2; }

            // Safeguarded Newton iteration: a Newton step inside the bracket,
            // otherwise bisection. It ends when the update no longer moves
            // the temperature or the bracket has closed.
            scalar t = Tc[i];
            const scalar tolT = 1e-12*max(mag(t), scalar(1));
            for (label it = 0; it < 60; it++)
            {
                const scalar f = energy(t) - Estar;
                if (f > 0) { hi = t; } else { lo = t; }
                const scalar tn = t - f/max(capacity(t), small);
                const scalar tNew = (tn > lo && tn < hi) ? tn : 0.5*(lo + hi);
                const bool done = mag(tNew - t) <= tolT || hi - lo <= tolT;
                t = tNew;
                if (done) { break; }
            }
            Tc[i] = t;
        }
        T.correctBoundaryConditions();

        forAll(alphaSolid, celli)
        {
            alphaSolid[celli] = alpha1c[celli]*metal.fs(T[celli]);
        }
        alphaSolid.correctBoundaryConditions();

        #include "thermo/updateProperties.H"

        // Largest change of the solid fraction in a cell by this iteration
        const scalar residual = gMax
        (
            mag
            (
                alphaSolid.primitiveField()
              - alphaSolidPrev.primitiveField()
            )
        );

        Info<< "Thermo: iteration " << tCorr
            << " residual: " << residual << endl;

        converged = residual < thermoTol && tCorr > 0;
    }

    if (!converged)
    {
        WarningInFunction
            << "Energy equation not converged in " << nThermoCorr
            << " iterations" << endl;
    }

    // Solid metal does not move, so the melted metal of a cell is its largest
    // liquid fraction, limited by the metal that remains
    alphaMelted = min(alpha1c, max(alphaMelted, alpha1c - alphaSolid));
    alphaMelted.correctBoundaryConditions();
}
