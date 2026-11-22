//- Evaporative and radiative cooling of the free surface, linearized about
//  the current temperature for the phase-change iteration. The sink acts in
//  the surface cells of the beam (the first metal cell of each column), on
//  the true area of the surface in the cell: its footprint divided by the
//  cosine of the surface inclination, from the normal of the metal fraction.
//  It is evaluated at the temperature of that metal cell, the temperature of
//  the liquid surface, as the recoil pressure is.
{
    surfaceCoolingSp.primitiveFieldRef() = 0;
    surfaceCoolingSu.primitiveFieldRef() = 0;
    evaporationPower = 0;
    radiationPower = 0;
    surfaceArea = 0;

    if (evaporationCoeff > 0 || emissivity > 0)
    {
        const scalar R = constant::physicoChemical::R.value();
        const scalar sigmaSB = constant::physicoChemical::sigma.value();
        const pointField& points = mesh.points();
        const volVectorField gradAlpha(fvc::grad(alpha1c));

        scalarField& Sp = surfaceCoolingSp.primitiveFieldRef();
        scalarField& Su = surfaceCoolingSu.primitiveFieldRef();

        forAll(surfaceCells, i)
        {
            const label celli = surfaceCells[i];
            const scalar t = T[celli];

            scalar zMin = great, zMax = -great;
            const labelList& vertices = mesh.cellPoints()[celli];
            forAll(vertices, j)
            {
                zMin = min(zMin, points[vertices[j]].z());
                zMax = max(zMax, points[vertices[j]].z());
            }
            const scalar cosTheta =
                max(mag(gradAlpha[celli].z())/(mag(gradAlpha[celli]) + vSmall), 0.2);
            const scalar area = mesh.V()[celli]/(zMax - zMin)/cosTheta;

            scalar q = 0, dqdT = 0;
            if (evaporationCoeff > 0)
            {
                const scalar qe =
                    evaporationCoeff*Foam::exp(TvapTheta*(1.0/Tvap - 1.0/t))/Foam::sqrt(t);
                q += qe;
                dqdT += qe*(TvapTheta/sqr(t) - 0.5/t);
                evaporationPower += qe*area;
            }
            if (emissivity > 0)
            {
                const scalar qr = emissivity*sigmaSB*(pow4(t) - pow4(Tambient));
                q += qr;
                dqdT += 4*emissivity*sigmaSB*pow3(t);
                radiationPower += qr*area;
            }
            surfaceArea += area;

            // q(T) = q(t) + dqdT*(T - t): the implicit part and the remainder
            Sp[celli] = dqdT*area/mesh.V()[celli];
            Su[celli] = (q - dqdT*t)*area/mesh.V()[celli];
        }

        reduce(evaporationPower, sumOp<scalar>());
        reduce(radiationPower, sumOp<scalar>());
        reduce(surfaceArea, sumOp<scalar>());
    }
}
