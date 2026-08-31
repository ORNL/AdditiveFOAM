/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Oak Ridge National Laboratory
     \\/     M anipulation  |
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

#include "KellyAbsorption.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace absorptionModels
{
    defineTypeNameAndDebug(Kelly, 0);
    addToRunTimeSelectionTable(absorptionModel, Kelly, dictionary);
}
}

Foam::absorptionModels::Kelly::Kelly
(
    const dictionary& dict,
    const fvMesh&
)
:
    geometry_(dict.lookup<word>("geometry")),
    eta0_("eta0", dimless, dict),
    etaMin_("etaMin", dimless, dict),
    aspectRatioSwitch_
    (
        dict.lookupOrDefault<scalar>("aspectRatioSwitch", 1.0)
    )
{
    if (geometry_ != "cone" && geometry_ != "cylinder")
    {
        FatalIOErrorInFunction(dict)
            << "Kelly geometry must be cone or cylinder, found " << geometry_
            << exit(FatalIOError);
    }
}

Foam::scalar Foam::absorptionModels::Kelly::eta
(
    const scalar aspectRatio
) const
{
    if (aspectRatio <= aspectRatioSwitch_)
    {
        return etaMin_.value();
    }

    const scalar theta = Foam::atan(1.0/aspectRatio);

    scalar F = 0;
    scalar G = 0;

    if (geometry_ == "cone")
    {
        F = 0.25*(3.0*Foam::sin(theta) - Foam::sin(3.0*theta));
        G = 1.0/(1.0 + Foam::sqrt(1.0 + sqr(aspectRatio)));
    }

    if (geometry_ == "cylinder")
    {
        F = 0.5*(1.0 - Foam::cos(2.0*theta));
        G = 0.5/(1.0 + aspectRatio);
    }

    return
        max
        (
            etaMin_.value(),
            (
                eta0_*(1.0 + (1.0 - eta0_)*(G - F))
              / (1.0 - (1.0 - eta0_)*(1.0 - G))
            ).value()
        );
}

// ************************************************************************* //
