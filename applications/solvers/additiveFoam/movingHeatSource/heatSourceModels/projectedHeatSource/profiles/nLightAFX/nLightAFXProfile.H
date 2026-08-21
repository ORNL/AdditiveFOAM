/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
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

Class
    Foam::heatSourceProfiles::nLightAFXProfile

Description
    Two-component nLight AFX planar heat source profile.

SourceFiles
    nLightAFXProfile.C

\*---------------------------------------------------------------------------*/

#ifndef nLightAFXProfile_H
#define nLightAFXProfile_H

#include "heatSourceProfile.H"

namespace Foam
{
namespace heatSourceProfiles
{

class nLightAFXProfile
:
    public heatSourceProfile
{
    // Private Data

        //- Fraction of power applied to Gaussian component 1
        scalar alpha_;

        //- Radius of Gaussian component 0
        scalar r0_;

        //- Standard deviation of Gaussian component 0
        scalar sigma0_;

        //- Radius of Gaussian component 1
        scalar r1_;

        //- Standard deviation of Gaussian component 1
        scalar sigma1_;

        //- J_0(0)
        scalar J0_;

        //- J_1(0)
        scalar J1_;

        //- Planar integral
        scalar integral_;

        //- Area-equivalent second-moment diameter
        scalar D4Sigma_;

        //- Planar profile bounds
        boundBox profileBb_;


    // Private Member Functions

        //- Return J(R) = integral_R^infinity I(r)*r dr
        inline scalar J
        (
            const scalar R,
            const scalar r,
            const scalar sigma
        ) const;

        //- Return K = integral_0^infinity I(r)*r^3 dr
        inline scalar K(const scalar r, const scalar sigma) const;

        //- Guaranteed upper bound where J(R)/J(0) <= epsilon
        inline scalar cutoff
        (
            const scalar r,
            const scalar sigma,
            const scalar J0,
            const scalar epsilon
        ) const;

public:

    //- Runtime type information
    TypeName("nLightAFX");


    // Constructors

        nLightAFXProfile
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const scalar tolerance
        );

    //- Destructor
    virtual ~nLightAFXProfile()
    {}

    // Member Functions

        inline virtual scalar weight
        (
            const scalar x,
            const scalar y
        ) const;

        inline virtual scalar integral() const
        {
            return integral_;
        }

        inline virtual scalar D4Sigma() const
        {
            return D4Sigma_;
        }

        inline virtual const boundBox& profileBb() const
        {
            return profileBb_;
        }

};

} // End namespace heatSourceProfiles
} // End namespace Foam

#endif

// ************************************************************************* //
