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
    // Select parameters from specified mode
    mode_ = heatSourceModelCoeffs_.lookup<label>("mode");

    List<scalar> alphas_ =
        {0.99, 0.83, 0.7, 0.54, 0.33, 0.21, 0.29};

    List<scalar> spot4Sigma_ =
        {156.99e-6, 158.23e-6, 165.31e-6, 182.28e-6, 181.15e-6, 214.63e-6, 402.04e-6};

    List<scalar> ring4Sigma_ =
        {48.49e-6, 112.15e-6, 116.07e-6, 135.32e-6, 116.74e-6, 112.46e-6, 116.62e-6};

    dimensions_ = heatSourceModelCoeffs_.lookup<vector>("dimensions");

    alpha_ = alphas_[mode_];

    Info<< "Parameters for mode: " << mode_ << endl
        << "\talpha: " << alpha_ << endl
        << "\tD4s (spot): " << spot4Sigma_[mode_] << endl
        << "\tD4s (ring): " << ring4Sigma_[mode_] << endl;
    
    //- Spot parameters
    ks_ = spotCoeffs_.lookup<scalar>("k");
    ms_ = spotCoeffs_.lookup<scalar>("m");
    ds_ = vector(0.5*spot4Sigma_[mode_], 0.5*spot4Sigma_[mode_], dimensions_[2]);

    //- Ring parameters
    kr_ = ringCoeffs_.lookup<scalar>("k");
    mr_ = ringCoeffs_.lookup<scalar>("m");
    R_ = 152e-6;
    r_ = 0.5 * ring4Sigma_[mode_];

    // Calculated using second moment method from fitted profiles
    List<scalar> D4sigma_ =
        {162.12e-6, 235.39e-6, 283.97e-6, 283.97e-6,
         340.0e-6, 385.18e-6, 413.13e-6, 438.85e-6};
    
    scalar R2sigma_ = 0.5*D4sigma_[mode_];    

    //- Overwrite the dimensions used by the base class
    dimensions_ = vector(R2sigma_, R2sigma_, ds_.z());
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
        return alpha_ / As().value() / dimensions_.z() * Foam::exp(-xs) 
            + (1.0 - alpha_) / Ar().value() / dimensions_.z() * Foam::exp(-xr);
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
Foam::heatSourceModels::nLight::As()
{
    const scalar a = Foam::pow(2.0, 1.0 / ks_);

    const vector s =
        cmptDivide
        (
            vector(ds_.x(), ds_.y(), dimensions_.z()),
            vector(a, a, 1.0)
        );

    const dimensionedScalar As
    (
        "As",
        dimArea,
        s.x() * s.y() * pi * Foam::tgamma(1.0 + 2.0 / ks_)
      * Foam::tgamma(1.0 + 1.0 / ms_) * Foam::tgamma(1.0 + 2.0 / ms_)
      / Foam::tgamma(1.0 + 3.0 / ms_)
    );

    return As;
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::nLight::Ar()
{
    const scalar a = Foam::pow(2.0, 1.0 / kr_);
    
    const vector s =
        cmptDivide
        (
            vector(r_, r_, dimensions_.z()),
            vector(a, a, 1.0)
        );
    
    const dimensionedScalar Ar
    (
        "Ar",
        dimArea,
        4.0 * pi * R_ * s.x() * Foam::tgamma(1.0 + 1.0 / kr_)
      * Foam::tgamma(1.0 + 1.0 / mr_) * Foam::tgamma(1.0 + 1.0 / mr_)
      / Foam::tgamma(1.0 + 2.0 / mr_)
    );
    
    return Ar;
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
