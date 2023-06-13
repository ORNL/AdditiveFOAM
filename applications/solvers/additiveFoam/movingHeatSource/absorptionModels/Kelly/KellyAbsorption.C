/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace absorptionModels
{
    defineTypeNameAndDebug(Kelly, 0);
    addToRunTimeSelectionTable(absorptionModel, Kelly, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absorptionModels::Kelly::Kelly
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionModel(typeName, sourceName, dict, mesh),
    geometry_(absorptionModelCoeffs_.lookup<word>("geometry")),
    eta0_("eta0", dimless, absorptionModelCoeffs_),
    etaMin_("etaMin", dimless, absorptionModelCoeffs_)
{
    if ((geometry_ != "cone") && (geometry_ != "cylinder"))
    {
        FatalErrorInFunction
        << "Kelly absorption model geometry type " << geometry_
        << " not recognized. Choose either cylinder or cone."
        << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar 
Foam::absorptionModels::Kelly::eta
(
    const scalar& aspectRatio
) const
{
    if (aspectRatio > 1.0)
    {
        const scalar theta = Foam::atan(1.0 / aspectRatio);
        
        scalar F = 0.0;
        scalar G = 0.0;
        
        if (geometry_ == "cone")
        {
            F = 0.25 * (3.0 * Foam::sin(theta) - Foam::sin(3.0 * theta));
            G = 1.0 / (1.0 + Foam::sqrt(1.0 + pow(aspectRatio, 2)));
        }
        
        else if (geometry_ == "cylinder")
        {
            F = 0.5 * (1.0 - Foam::cos(2.0 * theta));
            G = 0.5 / (1.0 + aspectRatio);
        }
        
        return (eta0_ * (1.0 + (1.0 - eta0_)*(G - F)) 
                / (1.0 - (1.0 - eta0_)*(1.0 - G))).value();
    }
    else
    {
        return etaMin_.value();
    }
};

bool Foam::absorptionModels::Kelly::read()
{
    if (absorptionModel::read())
    {
        absorptionModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        absorptionModelCoeffs_.lookup("geometry") >> geometry_;
        absorptionModelCoeffs_.lookup("eta0") >> eta0_;
        absorptionModelCoeffs_.lookup("etaMin") >> etaMin_;

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
