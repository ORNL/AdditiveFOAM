{
    //- Phase-change iteration with the latent heat treated as a sum of
    //  piecewise-linear ramps over the breakpoints of the thermo path, as in
    //  additiveFoam/thermo/TEqnNested.H. The sensible enthalpy of the metal
    //  and gas keeps its Newton linearization, so the iteration continues
    //  until the ramp sets are stable and the change of the solid fraction is
    //  below thermoTolerance. If the ramp sets repeat with period two (a
    //  limit cycle) or nNewtonSolves is reached, the ramp sets are iterated
    //  with the nested cycles of Brugnano and Casulli (SIAM J. Sci. Comput.
    //  31, 2009, Algorithm 1).

    const label nThermoCorr =
        pimple.dict().lookupOrDefault<label>("nThermoCorrectors", 20);
    const scalar thermoTol =
        pimple.dict().lookupOrDefault<scalar>("thermoTolerance", 1e-6);
    const label nNewtonSolves =
        pimple.dict().lookupOrDefault<label>("nNewtonSolves", nThermoCorr);
    const label nNestedCorr =
        pimple.dict().lookupOrDefault<label>("nNestedCorrectors", 50);

    const dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
    const surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    const scalar r1 = rho1.value();
    const scalar r2 = rho2.value();
    const scalar L = Lf.value();
    const scalar Ts = solidus.value();

    //- breakpoints and slope increments of the liquid fraction of the metal
    scalarField kinkT;
    scalarField kinkD;
    {
        scalarField Tp(path.x());
        const scalarField& ap = path.y();
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

    volScalarField dFdT
    (
        IOobject("dFdT", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless/dimTemperature, 0)
    );

    volScalarField alphaOffset
    (
        IOobject("alphaOffset", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    alphaSolid.oldTime();

    forAll(alphaSolid, celli)
    {
        alphaSolid[celli] = alpha1c[celli]*metal.fs(T[celli]);
    }
    alphaSolid.correctBoundaryConditions();

    #include "thermo/updateProperties.H"

    List<boolList> Pon(pKinks.size(), boolList(mesh.nCells(), true));
    List<boolList> Qon(qKinks.size(), boolList(mesh.nCells(), false));

    const scalar tolT = 1e-9;
    label nSolves = 0;
    label nClamped = 0;
    scalar residual = great;

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

    //- the solid fraction of the current ramp sets, alphaSolid = a1c (1 - a - b T)
    auto setCoeffs = [&](const volScalarField& alphaSolidPrev)
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

            const scalar ac = alpha1c[celli];
            dFdT[celli] = -ac*b;
            alphaOffset[celli] = ac*(1 - a) - alphaSolidPrev[celli];
        }

        dFdT.correctBoundaryConditions();
        alphaOffset.correctBoundaryConditions();
    };

    //- switch one family of ramps at the current temperature
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

    //- one linearized solve, projection of every cell onto the energy curve
    //  and property update; sets the residual
    auto iterate = [&]()
    {
        const volScalarField Tprev("Tprev", T);
        const volScalarField alphaSolidPrev("alphaSolidPrev", alphaSolid);

        setCoeffs(alphaSolidPrev);

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
          + rDeltaT*rho1*Lf*(fvm::Sp(dFdT, T) + alphaOffset)
        );

        TEqn.relax();
        TEqn.solve();
        nSolves++;

        // the cell energy follows from the solved linearization and T is put
        // back on the energy curve
        scalarField& Tc = T.primitiveFieldRef();

        forAll(Tc, i)
        {
            const scalar a = alpha1[i];
            const scalar ac = alpha1c[i];
            const scalar Tk = Tprev[i];

            const scalar alphaStar =
                alphaSolidPrev[i] + alphaOffset[i] + dFdT[i]*Tc[i];

            const scalar Estar =
                rhoH[i] + rhoCp[i]*(Tc[i] - Tk) + r1*L*(ac - alphaStar);

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

            scalar t = Tc[i];
            const scalar tolP = 1e-12*max(mag(t), scalar(1));
            for (label it = 0; it < 60; it++)
            {
                const scalar f = energy(t) - Estar;
                if (f > 0) { hi = t; } else { lo = t; }
                const scalar tn = t - f/max(capacity(t), small);
                const scalar tNew = (tn > lo && tn < hi) ? tn : 0.5*(lo + hi);
                const bool done = mag(tNew - t) <= tolP || hi - lo <= tolP;
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

        residual = gMax
        (
            mag
            (
                alphaSolid.primitiveField()
              - alphaSolidPrev.primitiveField()
            )
        );

        Info<< "Thermo: iteration " << nSolves - 1
            << " residual: " << residual << endl;
    };

    //- Newton: ramps switched at the current temperature
    switchRamps(pKinks, Pon);
    switchRamps(qKinks, Qon);

    bool converged = false;
    bool stalled = false;

    for (label m = 0; m < nNewtonSolves && !converged && !stalled; m++)
    {
        iterate();
        const bool changedP = switchRamps(pKinks, Pon);
        const bool changedQ = switchRamps(qKinks, Qon);
        converged = !(changedP || changedQ) && residual < thermoTol && m > 0;

        // sets that recur while still changing are a limit cycle
        const bool repeated = repeatedSets(residual);
        stalled = !converged && (changedP || changedQ) && repeated;
    }

    //- fallback: Algorithm 1 from the canonical starting sets
    if (!converged)
    {
        Info<< "Thermo: Newton " << (stalled ? "stalled" : "did not settle")
            << " in " << nSolves << " solves, nested iteration" << endl;

        forAll(qKinks, j)
        {
            Qon[j] = false;
        }

        bool outerDone = false;

        for (label k = 0; k < nNestedCorr && !outerDone; k++)
        {
            forAll(pKinks, j)
            {
                Pon[j] = true;
            }

            bool innerDone = false;

            for (label l = 0; l < nNestedCorr && !innerDone; l++)
            {
                iterate();
                innerDone = !switchRamps(pKinks, Pon) && residual < thermoTol;
            }

            outerDone = !switchRamps(qKinks, Qon) && innerDone;
        }

        converged = outerDone;
    }

    reduce(nClamped, sumOp<label>());

    if (!converged)
    {
        WarningInFunction
            << "Energy equation not converged in " << nSolves
            << " solves" << endl;
    }

    if (nClamped)
    {
        WarningInFunction
            << nClamped << " cell evaluations had a negative latent slope "
            << "and dropped their Q ramps" << endl;
    }

    // Solid metal does not move, so the melted metal of a cell is its largest
    // liquid fraction, limited by the metal that remains
    alphaMelted = min(alpha1c, max(alphaMelted, alpha1c - alphaSolid));
    alphaMelted.correctBoundaryConditions();
}
