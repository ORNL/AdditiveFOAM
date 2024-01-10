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

#include "modifiedSuperGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(modifiedSuperGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, modifiedSuperGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::modifiedSuperGaussian::modifiedSuperGaussian
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
    m_ = heatSourceModelCoeffs_.lookup<scalar>("m");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::modifiedSuperGaussian::weight(const vector& d)
{
    scalar a = Foam::pow(2.0, 1.0/k_);

    vector s = cmptDivide(dimensions_, vector(a, a, 1.0));

    if (d.z() < s.z())
    {
        s *= Foam::pow(1.0 - Foam::pow(d.z() / s.z(), m_), 1.0/m_);

        vector di = vector(d.x(), d.y(), 0.0);

        scalar x = Foam::pow(magSqr(cmptDivide(di, s)), k_/2.0);

        return Foam::exp(-x);
    }
    else
    {
        return 0.0;
    }
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::modifiedSuperGaussian::V0()
{
    const scalar a = Foam::pow(2.0, 1.0/k_);

    const vector s = cmptDivide(dimensions_, vector(a, a, 1.0));

    const dimensionedScalar V0
    (
        "V0",
        dimVolume,
        s.x()*s.y()*s.z()*pi*Foam::tgamma(1.0 + 2.0/k_)
      * Foam::tgamma(1.0 + 1.0/m_)*Foam::tgamma(1.0 + 2.0/m_)
      / Foam::tgamma(1.0 + 3.0/m_)
    );

    return V0;
}

bool Foam::heatSourceModels::modifiedSuperGaussian::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("k") >> k_;
        heatSourceModelCoeffs_.lookup("m") >> m_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
