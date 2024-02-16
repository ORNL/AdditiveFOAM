{
    const scalarField& xvals(thermo.x());
    const scalarField& yvals(thermo.y());
    label n = xvals.size();

    forAll(mesh.cells(), celli)
    {
        // update solid fraction estimate
        alpha10[celli] = alpha1[celli];

        // update slope and temperature estimates
        const scalar& x(T[celli]);

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

            T0[celli] = xvals[lo] + (alpha10[celli] - yvals[lo])/dFdT[celli];
        }
        else
        {
            dFdT[celli] = 0.0;
            T0[celli] = T[celli];
        }
    }

    // update boundary conditions
    alpha10.correctBoundaryConditions();
    dFdT.correctBoundaryConditions();
    T0.correctBoundaryConditions();
}
