/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseChangeSolver.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SensibleEnergy>
Foam::label Foam::phaseChangeSolver::solveCell
(
    const label celli,
    const scalar rhoHstar,
    const scalar TLo,
    const scalar THi,
    const SensibleEnergy& sensible,
    scalar& TCell,
    scalar& alphaSolid
) const
{
    const scalar rhoLf = rhoLf_.value();
    const scalar alphai = alpha(celli);

    // Enthalpy of the cell on the path at a trial temperature, less the
    // target: it increases with temperature at least as fast as the sensible
    // enthalpy does
    auto enthalpyResidual = [&](const scalar TTrial)
    {
        return
            sensible.rhoH(celli, TTrial)
          + rhoLf*alphai*(1 - path_.alpha(TTrial))
          - rhoHstar;
    };

    const scalar tolT = 1e-12*max(mag(TCell), scalar(1));

    scalar residual = enthalpyResidual(TCell);
    label it = 0;

    // A cell already at its enthalpy on the path is its own projection: every
    // cell outside the melting range that stays there, which is most of the
    // mesh on every solve
    if (mag(residual) > sensible.rhoCp(celli, TCell)*tolT)
    {
        // Bracket the root
        scalar TLow = min(TLo, TCell);
        scalar THigh = max(THi, TCell);
        scalar span = max(THigh - TLow, scalar(1));
        while (enthalpyResidual(TLow) > 0) { TLow -= span; span *= 2; }
        span = max(THigh - TLow, scalar(1));
        while (enthalpyResidual(THigh) < 0) { THigh += span; span *= 2; }

        // Safeguarded Newton iteration: a Newton step inside the bracket,
        // otherwise bisection. It ends when the enthalpy is matched to within
        // the tolerance, which bounds the temperature error by tolT because
        // the slope of the residual is at least the sensible rho*Cp, or when
        // the bracket has closed on an isothermal step of the path
        for (; it < maxCellIter; it++)
        {
            if (it > 0)
            {
                residual = enthalpyResidual(TCell);
            }

            const scalar rhoCpCell = sensible.rhoCp(celli, TCell);

            if (mag(residual) <= rhoCpCell*tolT)
            {
                break;
            }

            if (residual > 0) { THigh = TCell; } else { TLow = TCell; }

            const scalar TNewton =
                TCell
              - residual
               /max(rhoCpCell - rhoLf*alphai*path_.dAlphadT(TCell), small);

            const scalar TNext =
                (TNewton > TLow && TNewton < THigh)
              ? TNewton
              : 0.5*(TLow + THigh);
            const bool done =
                mag(TNext - TCell) <= tolT || THigh - TLow <= tolT;
            TCell = TNext;
            if (done) { break; }
        }
    }

    // The solid fraction that conserves the enthalpy; on an isothermal step of
    // the path it lies between the fractions on either side
    const scalar width = 1e-6*max(TCell, scalar(1));

    alphaSolid = min
    (
        max
        (
            alphai - (rhoHstar - sensible.rhoH(celli, TCell))/rhoLf,
            alphai*path_.alpha(TCell + width)
        ),
        alphai*path_.alpha(TCell - width)
    );

    return it;
}


template<class SensibleEnergy>
void Foam::phaseChangeSolver::project(const SensibleEnergy& sensible)
{
    const scalar rhoLf = rhoLf_.value();

    scalarField& T = T_.primitiveFieldRef();
    scalarField& alphaSolid = alphaSolid_.primitiveFieldRef();

    scalar residual = 0;
    label nCapped = 0;

    forAll(T, celli)
    {
        const scalar T0 = T0_[celli];
        const scalar Tstar = T[celli];

        // Solid fraction and enthalpy of the linearized solution
        const scalar alphaStar =
            alphaSolid0_[celli] + dFdT_[celli]*(Tstar - T0);

        const scalar rhoHstar =
            sensible.rhoHstar(celli, T0, Tstar)
          + rhoLf*(alpha(celli) - alphaStar);

        scalar TCell = Tstar;

        if
        (
            solveCell
            (
                celli, rhoHstar, T0, Tstar, sensible, TCell, alphaSolid[celli]
            )
         == maxCellIter
        )
        {
            nCapped++;
        }

        T[celli] = TCell;

        residual = max(residual, mag(alphaSolid[celli] - alphaSolid0_[celli]));
    }

    T_.correctBoundaryConditions();
    alphaSolid_.correctBoundaryConditions();

    reduce(residual, maxOp<scalar>());
    reduce(nCapped, sumOp<label>());

    residual_ = residual;
    converged_ = (residual_ < tol_) && (iter_ > 1);

    Info<< "Thermo: iteration " << iter_ - 1
        << " residual: " << residual_ << endl;

    if (nCapped)
    {
        WarningInFunction
            << nCapped << " cells were not put back on the thermoPath in "
            << maxCellIter << " iterations" << endl;
    }
}


template<class SensibleEnergy>
void Foam::phaseChangeSolver::setFromEnthalpy
(
    const volScalarField& rhoH,
    const SensibleEnergy& sensible
)
{
    scalarField& T = T_.primitiveFieldRef();
    scalarField& alphaSolid = alphaSolid_.primitiveFieldRef();

    label nCapped = 0;

    forAll(T, celli)
    {
        scalar TCell = T[celli];

        if
        (
            solveCell
            (
                celli, rhoH[celli], TCell, TCell, sensible, TCell,
                alphaSolid[celli]
            )
         == maxCellIter
        )
        {
            nCapped++;
        }

        T[celli] = TCell;
    }

    T_.correctBoundaryConditions();
    alphaSolid_.correctBoundaryConditions();

    reduce(nCapped, sumOp<label>());

    if (nCapped)
    {
        WarningInFunction
            << nCapped << " cells were not put back on the thermoPath in "
            << maxCellIter << " iterations" << endl;
    }
}


// ************************************************************************* //
