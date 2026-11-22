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
#include "fvcDdt.H"
#include "fvmDdt.H"
#include "EulerDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::phaseChangeSolver::maxCellIter;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::phaseChangeSolver::linearize()
{
    scalarField& T0 = T0_.primitiveFieldRef();
    scalarField& dFdT = dFdT_.primitiveFieldRef();

    forAll(T_, celli)
    {
        const scalar T = T_[celli];
        const scalar alphai = alpha(celli);

        T0[celli] = T;

        // Slope of the path at T. On an isothermal step of the path, with the
        // solid fraction strictly inside the step, the segment is vertical:
        // its slope over the width of the step as the projection resolves it
        // keeps the temperature on the step while the enthalpy goes into the
        // solid fraction
        const scalar width = 1e-6*max(T, scalar(1));

        dFdT[celli] = alphai*path_.dAlphadT(T);

        if (alphai*path_.step(T, width) > small)
        {
            const scalar alphaHigh = alphai*path_.alpha(T - width);
            const scalar alphaLow = alphai*path_.alpha(T + width);
            const scalar a = alphaSolid_[celli];

            if (a > alphaLow + small && a < alphaHigh - small)
            {
                dFdT[celli] = -(alphaHigh - alphaLow)/(2*width);
            }
        }
    }

    T0_.correctBoundaryConditions();
    dFdT_.correctBoundaryConditions();

    alphaSolid0_ = alphaSolid_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeSolver::phaseChangeSolver
(
    volScalarField& T,
    volScalarField& alphaSolid,
    const thermoPath& path,
    const dimensionedScalar& rhoLf,
    const volScalarField* alpha
)
:
    mesh_(T.mesh()),
    T_(T),
    alphaSolid_(alphaSolid),
    path_(path),
    rhoLf_(rhoLf),
    alphaPtr_(alpha),
    T0_
    (
        IOobject("T0", T.time().name(), mesh_),
        T
    ),
    dFdT_
    (
        IOobject("dFdT", T.time().name(), mesh_),
        mesh_,
        dimensionedScalar(dimless/dimTemperature, 0)
    ),
    alphaSolid0_
    (
        IOobject("alphaSolid0", T.time().name(), mesh_),
        alphaSolid
    ),
    nCorr_(1),
    tol_(0),
    iter_(0),
    residual_(great),
    converged_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseChangeSolver::begin(const label nCorr, const scalar tol)
{
    // Store the old-time fields before the cell updates modify them
    T_.oldTime();
    alphaSolid_.oldTime();

    nCorr_ = nCorr;
    tol_ = tol;
    iter_ = 0;
    residual_ = great;
    converged_ = false;
}


bool Foam::phaseChangeSolver::loop()
{
    if (converged_)
    {
        return false;
    }

    if (iter_ >= nCorr_)
    {
        if (nCorr_ > 1)
        {
            WarningInFunction
                << "Phase change not converged in " << nCorr_
                << " iterations, residual " << residual_ << endl;
        }

        return false;
    }

    linearize();
    iter_++;

    return true;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::phaseChangeSolver::latentSource
(
    const bool eulerDdt
) const
{
    // The implicit part is the derivative of the latent heat with temperature
    // under the same time-derivative operator as the sensible heat, so the
    // linearization is consistent with the scheme; the old-time terms of the
    // operator cancel between its implicit and explicit forms
    if (eulerDdt)
    {
        fv::EulerDdtScheme<scalar> euler(mesh_);

        return
            rhoLf_*euler.fvcDdt(alphaSolid_)
          + rhoLf_*(euler.fvmDdt(dFdT_, T_) - euler.fvcDdt(dFdT_, T_));
    }

    return
        rhoLf_*fvc::ddt(alphaSolid_)
      + rhoLf_*(fvm::ddt(dFdT_, T_) - fvc::ddt(dFdT_, T_));
}


// ************************************************************************* //
