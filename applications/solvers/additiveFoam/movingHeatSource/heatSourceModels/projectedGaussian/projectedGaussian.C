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

#include "projectedGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(projectedGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, projectedGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::projectedGaussian::projectedGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh)
{
    A_ = heatSourceModelCoeffs_.lookup<scalar>("A");
    B_ = heatSourceModelCoeffs_.lookup<scalar>("B");
    
    // set initial shape function
    const scalar x_ = max(dimensions_.z() / staticDimensions_.x(), 1.0);
    const scalar n_ = min(max(A_*std::log2(x_) + B_, 0.0), 9.0);
    k_ = std::pow(2.0, n_);            
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::projectedGaussian::weight(const vector& d)
{    
    const scalar f_ =
        std::exp
        (
            -2.0
          * (
                Foam::sqr(d.x() / dimensions_.x())
              + Foam::sqr(d.y() / dimensions_.y())
            )
        );
    
    const scalar s_ =
        std::exp(-3.0 * std::pow(mag(mag(d.z()) / dimensions_.z()), k_));
           
    return f_ * s_;
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::projectedGaussian::V0()
{
    const scalar x_ = max(dimensions_.z() / staticDimensions_.x(), 1.0);
    
    const scalar n_ = min(max(A_*std::log2(x_) + B_, 0.0), 9.0);
    
    k_ = std::pow(2.0, n_);
    
    const dimensionedScalar V0
    (
        "V0",
        dimVolume,
        0.5 * pi * dimensions_.x() * dimensions_.y() * dimensions_.z()
      * Foam::tgamma(1.0 / k_)
      / ( k_ * std::pow(3.0, 1.0 / k_) )
    );

    return V0;
}

bool Foam::heatSourceModels::projectedGaussian::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("A") >> A_;
        heatSourceModelCoeffs_.lookup("B") >> B_;
        
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
