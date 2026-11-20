{
    #include "thermo/thermoScheme.H"

    //- Phase-change iteration with a convergence guarantee. The Newton
    //  iteration of the optimum approach (thermo/TEqn.H: slope of the path
    //  at the iterate, each solve followed by the energy-conserving
    //  projection of every cell onto the path) is run first. If its active
    //  path segments repeat with period two (a limit cycle: the same linear
    //  system recurs, so it cannot progress) or nNewtonSolves is reached, the
    //  step falls back to the nested iteration of Brugnano and Casulli (SIAM
    //  J. Sci. Comput. 31, 2009, Algorithm 1) on the piecewise-linear path:
    //  the liquid fraction lambda = 1 - alpha is written as a sum of ramps
    //  over the breakpoints, lambda(T) = sum_j d_j (T - T_j)+, with slope
    //  increments d_j. Ramps with d_j > 0 (set P) switch in an inner cycle
    //  that starts with all of them active, ramps with d_j < 0 (set Q) in an
    //  outer cycle that starts with none active. With the ramps frozen the
    //  energy equation is linear with a Stieltjes matrix, and the nested
    //  cycles reach the exact solution of the piecewise-linear problem in
    //  finitely many solves (their Theorems 1 and 2). Isothermal steps are
    //  given a small width for the ramps and are made exact by the final
    //  projection onto the path.

    const label nNewtonSolves =
        pimple.dict().lookupOrDefault<label>("nNewtonSolves", label(nThermoCorr));
    const label nNestedCorr =
        pimple.dict().lookupOrDefault<label>("nNestedCorrectors", 50);
    const bool debugPC =
        pimple.dict().lookupOrDefault<Switch>("phaseChangeDebug", false);

    scalarField kinkT;
    scalarField kinkD;
    {
        scalarField Tp(thermo.x());
        const scalarField& ap = thermo.y();
        const scalar epsIso = 1e-3;

        for (label i = 1; i < Tp.size(); i++)
        {
            if (Tp[i] <= Tp[i - 1])
            {
                Tp[i] = Tp[i - 1] + epsIso;
            }
        }

        kinkT.setSize(Tp.size());
        kinkD.setSize(Tp.size());

        scalar prev = 0;

        forAll(Tp, i)
        {
            const scalar slope =
                (i < Tp.size() - 1)
              ? (ap[i] - ap[i + 1])/(Tp[i + 1] - Tp[i])
              : 0;

            kinkT[i] = Tp[i];
            kinkD[i] = slope - prev;
            prev = slope;
        }
    }

    DynamicList<label> pKinks;
    DynamicList<label> qKinks;

    forAll(kinkD, j)
    {
        if (kinkD[j] > small)
        {
            pKinks.append(j);
        }
        else if (kinkD[j] < -small)
        {
            qKinks.append(j);
        }
    }

    if (debugPC)
    {
        Info<< "Thermo: kinks T " << kinkT << nl
            << "Thermo: kinks d " << kinkD << nl
            << "Thermo: P " << pKinks << " Q " << qKinks << endl;
    }

    volScalarField dFdT
    (
        IOobject("dFdT", runTime.name(), mesh),
        mesh,
        dimensionedScalar("dFdT", dimless/dimTemperature, 0.0)
    );

    volScalarField alphaOffset
    (
        IOobject("alphaOffset", runTime.name(), mesh),
        mesh,
        dimensionedScalar("alphaOffset", dimless, 0.0)
    );

    //- store the old-time fields before the cell updates below modify them
    T.oldTime();
    alpha1.oldTime();

    const scalar L = Lf.value();

    //- optional implicit limiting of temperature field
    volScalarField A
    (
        "A",
        dimensionedScalar(rho.dimensions()*Cp.dimensions(), 1e15)
       *pos(T - Tmax)
    );

    List<boolList> Pon(pKinks.size(), boolList(mesh.nCells(), true));
    List<boolList> Qon(qKinks.size(), boolList(mesh.nCells(), false));

    const scalar tolT = 1e-9;
    label nSolves = 0;
    label nClamped = 0;
    scalar residual = great;

    //- the latent fraction of the current ramp sets, lambda = a + b T, as
    //  an offset from the current alpha1 so that together with ddt(alpha1)
    //  the source is rho Lf (alphaLin(T) - alpha1.oldTime())/deltaT
    auto setCoeffs = [&]()
    {
        forAll(T, celli)
        {
            scalar a = 0;
            scalar b = 0;

            forAll(pKinks, j)
            {
                if (Pon[j][celli])
                {
                    const label i = pKinks[j];
                    b += kinkD[i];
                    a -= kinkD[i]*kinkT[i];
                }
            }

            scalar aq = 0;
            scalar bq = 0;

            forAll(qKinks, j)
            {
                if (Qon[j][celli])
                {
                    const label i = qKinks[j];
                    bq += kinkD[i];
                    aq -= kinkD[i]*kinkT[i];
                }
            }

            // the total slope of consistent sets is non-negative; it is zero
            // for cells above the liquidus, so round-off must not count
            if (b + bq < -1e-9*(mag(b) + mag(bq)))
            {
                nClamped++;
            }
            a += aq;
            b = max(b + bq, scalar(0));

            dFdT[celli] = -b;
            alphaOffset[celli] = 1 - a - alpha1[celli];
        }

        dFdT.correctBoundaryConditions();
        alphaOffset.correctBoundaryConditions();
    };

    //- solve the energy equation with the frozen latent heat
    auto solveEnergy = [&]()
    {
        solve
        (
            TEqn + rDeltaT*(fvm::Sp(A, T) - A*Tmax)
         ==
            (
                (ddtScheme == "Euler")
              ? rho*Lf*fv::EulerDdtScheme<scalar>(mesh).fvcDdt(alpha1)
              : rho*Lf*fvc::ddt(alpha1)
            )
          + rDeltaT*rho*Lf*(fvm::Sp(dFdT, T) + alphaOffset)
        );

        nSolves++;
    };

    //- put each cell on the path at the energy of the linear solution; the
    //  bracket starts between the previous iterate and the solution, as in
    //  TEqn.H, so a cell that did not move is done in one evaluation
    auto project = [&](const volScalarField& Tref)
    {
        forAll(T, celli)
        {
            const scalar Tstar = T[celli];
            const scalar c = Cp[celli]/L;
            const scalar alphaStar =
                alpha1[celli] + alphaOffset[celli] + dFdT[celli]*Tstar;

            //- a cell whose linearized state is on the path is its own
            //  projection (as in TEqn.H)
            if (thermo.alpha(Tstar) == alphaStar)
            {
                alpha1[celli] = alphaStar;
                continue;
            }

            auto f = [&](const scalar t)
            {
                return c*(t - Tstar) - (thermo.alpha(t) - alphaStar);
            };

            scalar lo = min(Tref[celli], Tstar);
            scalar hi = max(Tref[celli], Tstar);
            scalar span = max(hi - lo, scalar(1));
            while (f(lo) > 0) { lo -= span; span *= 2; }
            span = max(hi - lo, scalar(1));
            while (f(hi) < 0) { hi += span; span *= 2; }

            scalar t = Tstar;
            const scalar tolP = 1e-12*max(mag(t), scalar(1));
            for (label it = 0; it < 60; it++)
            {
                const scalar g = f(t);
                if (g > 0) { hi = t; } else { lo = t; }
                const scalar tn = t - g/(c - thermo.dAlphadT(t));
                const scalar tNew = (tn >= lo && tn <= hi) ? tn : 0.5*(lo + hi);
                const bool done = mag(tNew - t) <= tolP || hi - lo <= tolP;
                t = tNew;
                if (done) { break; }
            }

            const scalar eps = 1e-6*max(t, scalar(1));
            T[celli] = t;
            alpha1[celli] = min
            (
                max(alphaStar + c*(t - Tstar), thermo.alpha(t + eps)),
                thermo.alpha(t - eps)
            );
        }

        T.correctBoundaryConditions();
        alpha1.correctBoundaryConditions();
    };

    //- switch one family of ramps at the current temperature; true if any
    //  cell changed
    auto switchRamps = [&](const DynamicList<label>& kinks, List<boolList>& on)
    {
        bool changed = false;

        forAll(kinks, j)
        {
            const scalar Tj = kinkT[kinks[j]];

            forAll(T, celli)
            {
                if (mag(T[celli] - Tj) > tolT*max(mag(Tj), scalar(1)))
                {
                    const bool active = T[celli] > Tj;

                    if (active != on[j][celli])
                    {
                        on[j][celli] = active;
                        changed = true;
                    }
                }
            }
        }

        reduce(changed, orOp<bool>());

        return changed;
    };

    //- signatures of the ramp sets of every cell after the last solves, to
    //  detect repeated sets (limit cycles of period up to maxPeriod)
    const label maxPeriod = 6;
    List<uint64_t> sig0(mesh.nCells(), uint64_t(0));
    List<List<uint64_t>> sigHistory;
    DynamicList<scalar> resHistory;

    auto signature = [&](List<uint64_t>& sig)
    {
        forAll(sig, celli)
        {
            uint64_t h = 1;
            forAll(pKinks, j)
            {
                h = h*31u + (Pon[j][celli] ? 2u : 1u);
            }
            forAll(qKinks, j)
            {
                h = h*31u + (Qon[j][celli] ? 2u : 1u);
            }
            sig[celli] = h;
        }
    };

    //- true if the sets after the current solve equal those after an
    //  earlier solve of this step, two or more solves ago, with a residual
    //  that has not halved since: the same linear system recurs, so the
    //  iteration is periodic and cannot progress
    auto repeatedSets = [&](const scalar res)
    {
        signature(sig0);
        label period = 0;
        for (label q = sigHistory.size() - 2; q >= 0 && !period; q--)
        {
            if (res < 0.5*resHistory[q]) { continue; }
            bool same = true;
            const List<uint64_t>& old = sigHistory[q];
            forAll(sig0, celli)
            {
                if (sig0[celli] != old[celli]) { same = false; break; }
            }
            if (same) { period = sigHistory.size() - q; }
        }
        reduce(period, maxOp<label>());
        sigHistory.append(sig0);
        resHistory.append(res);
        if (sigHistory.size() > maxPeriod + 1)
        {
            for (label q = 0; q < sigHistory.size() - 1; q++)
            {
                sigHistory[q] = sigHistory[q + 1];
                resHistory[q] = resHistory[q + 1];
            }
            sigHistory.setSize(maxPeriod + 1);
            resHistory.setSize(maxPeriod + 1);
        }
        return period > 0;
    };

    auto report = [&]()
    {
        Info<< "Thermo: iteration " << nSolves - 1
            << " residual: " << residual << endl;

        if (debugPC && mesh.nCells() > 0)
        {
            Info<< "Thermo:   cell0 T " << T[0] << " alpha " << alpha1[0]
                << " dFdT " << dFdT[0] << " P";
            forAll(pKinks, j) { Info<< (Pon[j][0] ? 1 : 0); }
            Info<< " Q";
            forAll(qKinks, j) { Info<< (Qon[j][0] ? 1 : 0); }
            Info<< endl;
        }
    };

    //- Newton phase: the optimum approach with cycle detection
    switchRamps(pKinks, Pon);
    switchRamps(qKinks, Qon);

    bool converged = false;
    bool stalled = false;

    for (label m = 0; m < nNewtonSolves && !converged && !stalled; m++)
    {
        const volScalarField alphaPrev("alphaPrev", alpha1);
        const volScalarField Tprev("Tprev", T);

        setCoeffs();
        solveEnergy();
        project(Tprev);

        residual = gMax(mag(alpha1.primitiveField() - alphaPrev.primitiveField()));
        report();

        const bool changedP = switchRamps(pKinks, Pon);
        const bool changedQ = switchRamps(qKinks, Qon);
        converged = (residual < thermoTol && m > 0) || nThermoCorr <= 1;

        // sets that recur while still changing are a limit cycle
        const bool repeated = repeatedSets(residual);
        stalled = !converged && (changedP || changedQ) && repeated;
    }

    //- fallback: Algorithm 1 from the canonical starting sets, on the
    //  piecewise-linear path without intermediate projection
    if (!converged)
    {
        Info<< "Thermo: Newton " << (stalled ? "stalled" : "did not settle")
            << " in " << nSolves << " solves, nested iteration" << endl;

        forAll(qKinks, j)
        {
            Qon[j] = false;
        }

        bool outerDone = false;
        bool innerDone = false;

        for (label k = 0; k < nNestedCorr && !outerDone; k++)
        {
            forAll(pKinks, j)
            {
                Pon[j] = true;
            }

            innerDone = false;

            for (label l = 0; l < nNestedCorr && !innerDone; l++)
            {
                setCoeffs();
                solveEnergy();

                // mismatch between the path and the frozen ramps
                residual = 0;
                forAll(T, celli)
                {
                    residual = max
                    (
                        residual,
                        mag
                        (
                            thermo.alpha(T[celli])
                          - (alpha1[celli] + alphaOffset[celli] + dFdT[celli]*T[celli])
                        )
                    );
                }
                reduce(residual, maxOp<scalar>());
                report();

                innerDone = !switchRamps(pKinks, Pon);
            }

            outerDone = !switchRamps(qKinks, Qon);
        }

        converged = outerDone && innerDone;

        project(T.oldTime());
    }

    reduce(nClamped, sumOp<label>());

    if (!converged)
    {
        WarningInFunction
            << "Phase-change iteration not converged in " << nSolves
            << " solves" << endl;
    }

    if (nClamped)
    {
        WarningInFunction
            << nClamped << " cell evaluations had a negative latent slope"
            << endl;
    }
}
