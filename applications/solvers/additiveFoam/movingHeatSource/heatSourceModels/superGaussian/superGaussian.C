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

#include "superGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(superGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, superGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::superGaussian::superGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh)
{
    k_ = heatSourceModelCoeffs_.lookup<scalar>("k");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::superGaussian::weight(const vector& d)
{
    vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

    scalar x = Foam::pow(magSqr(cmptDivide(d, s)), k_/2.0);

    return Foam::exp(-x);
}


inline Foam::dimensionedScalar
Foam::heatSourceModels::superGaussian::V0()
{
    vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

    const dimensionedScalar V0
    (
        "V0",
        dimVolume,
        (2.0 / 3.0)*s.x()*s.y()*s.z()*pi*Foam::tgamma(1.0 + 3.0/k_)
    );

    return V0;
}


bool Foam::heatSourceModels::superGaussian::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("k") >> k_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
