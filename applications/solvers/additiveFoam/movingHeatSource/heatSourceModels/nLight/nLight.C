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

#include "nLight.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(nLight, 0);
    addToRunTimeSelectionTable(heatSourceModel, nLight, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::nLight::nLight
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh),
    spotCoeffs_(heatSourceModelCoeffs_.optionalSubDict("spotCoeffs")),
    ringCoeffs_(heatSourceModelCoeffs_.optionalSubDict("ringCoeffs"))
{
    alpha_ = heatSourceModelCoeffs_.lookup<scalar>("alpha");
    
    //- Spot parameters
    ks_ = spotCoeffs_.lookup<scalar>("k");
    ms_ = spotCoeffs_.lookup<scalar>("m");
    ds_ = heatSourceModelCoeffs_.lookup<vector>("dimensions");
    
    //- Ring parameters
    kr_ = ringCoeffs_.lookup<scalar>("k");
    mr_ = ringCoeffs_.lookup<scalar>("m");
    R_ = ringCoeffs_.lookup<scalar>("R");
    r_ = ringCoeffs_.lookup<scalar>("r");
    
    //- Overwrite the dimensions used by the base class
    dimensions_ = vector(R_ + r_, R_ + r_, ds_.z());
    staticDimensions_ = dimensions_;
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::nLight::weight(const vector& d)
{
    // Spot parameters
    const scalar as = Foam::pow(2.0, 1.0/ks_);
    vector ss =
        cmptDivide
        (
            vector(ds_.x(), ds_.y(), dimensions_.z()),
            vector(as, as, 1.0)
        );
    
    // Ring parameters
    const scalar ar = Foam::pow(2.0, 1.0/kr_);
    vector sr =
        cmptDivide
        (
            vector(r_, r_, dimensions_.z()),
            vector(ar, ar, 1.0)
        );

    // Check if cell is within heat source depth
    if (d.z() < max(ss.z(), sr.z()))
    {
        // x-y distance from beam center
        vector di = vector(d.x(), d.y(), 0.0);
        
        // Spot calcs
        ss *= Foam::pow(1.0 - Foam::pow(d.z() / ss.z(), ms_), 1.0/ms_);

        scalar xs = Foam::pow(magSqr(cmptDivide(di, ss)), ks_/2.0);
        
        // Ring calcs
        const scalar rsz =
            sr.x() * Foam::pow
                     (
                         1.0 - Foam::pow(d.z() / sr.z(), mr_),
                         1.0 / mr_
                     );
        
        scalar xr =
            Foam::pow
            (
                Foam::sqr((mag(di) - R_) / rsz), 
                kr_ / 2.0
            );

        // Return weights with power split and volume factored in
        return alpha_ / Vs().value() * Foam::exp(-xs) 
            + (1.0 - alpha_) / Vr().value() * Foam::exp(-xr);
    }
    else
    {
        return 0.0;
    }
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::nLight::V0()
{
    //- Because of linear combination of shapes in nLight,
    //  volume is integrated into the weight calculation
    return dimensionedScalar("V0", dimVolume, 1.0);
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::nLight::Vs()
{
    const scalar a = Foam::pow(2.0, 1.0/ks_);

    const vector s =
        cmptDivide
        (
            vector(ds_.x(), ds_.y(), dimensions_.z()),
            vector(a, a, 1.0)
        );

    const dimensionedScalar Vs
    (
        "Vs",
        dimVolume,
        s.x()*s.y()*s.z()*pi*Foam::tgamma(1.0 + 2.0/ks_)
      * Foam::tgamma(1.0 + 1.0/ms_)*Foam::tgamma(1.0 + 2.0/ms_)
      / Foam::tgamma(1.0 + 3.0/ms_)
    );

    return Vs;
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::nLight::Vr()
{
    const scalar a = Foam::pow(2.0, 1.0/kr_);
    
    const vector s =
        cmptDivide
        (
            vector(r_, r_, dimensions_.z()),
            vector(a, a, 1.0)
        );
    
    const dimensionedScalar Vr
    (
        "Vr",
        dimVolume,
        4.0 * pi * R_ * s.x() * s.z() * pi / kr_ // this kr is in original nLight impl but not sure where it comes from
      * Foam::tgamma(1.0 + 1.0 / mr_) * Foam::tgamma(1.0 + 1.0 / mr_)
      * Foam::tgamma(1.0 + 1.0 / kr_) / Foam::tgamma(1.0 + 2.0 / mr_)
    );
    
    return Vr;
}

bool Foam::heatSourceModels::nLight::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- General parameters
        heatSourceModelCoeffs_.lookup("alpha") >> alpha_;
        
        //- Spot parameters
        spotCoeffs_.lookup("k") >> ks_;
        spotCoeffs_.lookup("m") >> ms_;
        
        //- Ring parameters
        ringCoeffs_.lookup("k") >> kr_;
        ringCoeffs_.lookup("m") >> kr_;
        ringCoeffs_.lookup("R") >> R_;
        ringCoeffs_.lookup("r") >> r_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
