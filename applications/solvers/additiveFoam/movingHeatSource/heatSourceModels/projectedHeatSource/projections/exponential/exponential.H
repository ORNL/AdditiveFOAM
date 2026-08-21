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
    Foam::heatSourceProjections::exponential

Description
    One-sided generalized exponential projection applied in the negative z
    direction from the beam plane.

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "heatSourceProjection.H"

namespace Foam
{
namespace heatSourceProjections
{

class exponential
:
    public heatSourceProjection
{
    // Private Data

        //- nSlope in n = nSlope*log2(a) + nIntercept
        scalar nSlope_;

        //- nIntercept in n = nSlope*log2(a) + nIntercept
        scalar nIntercept_;

        //- Projection depth
        scalar d_;

        //- Axial decay exponent
        scalar k_;

        //- Projection integral
        scalar integral_;

        //- Lower projection bound
        scalar zMin_;

public:

    //- Runtime type information
    TypeName("exponential");


    // Constructors

        explicit exponential(const dictionary& dict);

    //- Destructor
    virtual ~exponential()
    {}

    // Member Functions

        virtual void update
        (
            const scalar depth,
            const scalar aspectRatio,
            const scalar tolerance
        );

        inline virtual scalar weight(const scalar z) const;

        inline virtual scalar integral() const
        {
            return integral_;
        }

        inline virtual scalar zMin() const
        {
            return zMin_;
        }
};

} // End namespace heatSourceProjections
} // End namespace Foam

#endif

// ************************************************************************* //
