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
    Foam::heatSourceProfiles::superGaussianProfile

Description
    Planar super-Gaussian heat source profile.

SourceFiles
    superGaussianProfile.C

\*---------------------------------------------------------------------------*/

#ifndef superGaussianProfile_H
#define superGaussianProfile_H

#include "heatSourceProfile.H"
#include "vector2D.H"

namespace Foam
{
namespace heatSourceProfiles
{

class superGaussianProfile
:
    public heatSourceProfile
{
    // Private Data

        //- Heat source dimensions
        vector2D dimensions_;

        //- Shape coefficient
        scalar k_;

        //- Scaled heat source dimensions
        vector2D s_;

        //- Planar integral
        scalar integral_;

        //- Area-equivalent second-moment diameter
        scalar D4Sigma_;

        //- Planar profile bounds
        boundBox profileBb_;

public:

    //- Runtime type information
    TypeName("superGaussian");


    // Constructors

        superGaussianProfile
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const scalar tolerance
        );

    //- Destructor
    virtual ~superGaussianProfile()
    {}

    // Member Functions

        //- Return the area-equivalent second-moment diameter
        static scalar D4Sigma
        (
            const scalar x,
            const scalar y,
            const scalar k
        );

        virtual scalar weight
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
