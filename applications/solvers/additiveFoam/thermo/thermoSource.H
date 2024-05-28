{
    scalarField xvals(thermo.x());
    const scalarField& yvals(thermo.y());
    label n = xvals.size();
    const scalar slope = 1e10;

    forAll(mesh.cells(), celli)
    {
        // update slope and temperature estimates
        const scalar& x = T[celli];
        const scalar& y = alpha1[celli];

        if ((x < Tliq.value()) && (x > Tsol.value()))
        {
            label lo = 0;
            for (lo=0; lo<n && xvals[lo]>x; ++lo)
            {}

            label low = lo;
            if (low < n)
            {
                for (label i=low; i<n; ++i)
                {
                    if (xvals[i] > xvals[lo] && xvals[i] <= x)
                    {
                        lo = i;
                    }
                }
            }

            label hi = 0;
            for (hi=0; hi<n && xvals[hi]<x; ++hi)
            {}

            label high = hi;
            if (high < n)
            {
                for (label i=high; i<n; ++i)
                {
                    if (xvals[i] < xvals[hi] && xvals[i] >= x)
                    {
                        hi = i;
                    }
                }
            }

            dFdT[celli] = (yvals[hi] - yvals[lo]) / (xvals[hi] - xvals[lo]);

            T0[celli] = xvals[lo] + (y - yvals[lo])/dFdT[celli];
        }
        else if ((x >= Tliq.value()) && (y > thermoTol))
        {
            dFdT[celli] = slope;
            T0[celli] = Tliq.value() - thermoTol;
        }
        else if ((x <= Tsol.value()) && (y < 1 - thermoTol))
        {
            dFdT[celli] = slope;
            T0[celli] = Tsol.value() + thermoTol;
        }
        else
        {
            dFdT[celli] = 0.0;
            T0[celli] = x;
        }
    }

    dFdT.correctBoundaryConditions();
    T0.correctBoundaryConditions();
}
