/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
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
    Foam::heatSourceModels::projectedGaussian

Description
    One-sided projected Gaussian heat source model applied in the negative z
    direction from the beam plane. The beam plane is expected to coincide with
    the material surface.

SourceFiles
    projectedGaussian.C

\*---------------------------------------------------------------------------*/

#ifndef projectedGaussian_H
#define projectedGaussian_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                       Class projectedGaussian
\*---------------------------------------------------------------------------*/

class projectedGaussian
:
    public heatSourceModel
{
    // Private Data

        //- Pointer to mesh information
        const fvMesh& mesh_;

        //- nSlope in n = nSlope*log2(x) + nIntercept
        scalar nSlope_;

        //- nIntercept in n = nSlope*log2(x) + nIntercept
        scalar nIntercept_;

        //- Axial decay exponent
        scalar k_;

public:

    //- Runtime type information
    TypeName("projectedGaussian");


    // Constructors

        //- Construct from components
        projectedGaussian
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~projectedGaussian()
    {}


    // Member Functions

        inline virtual scalar weight(const vector& d);

        inline virtual dimensionedScalar V0();

        //- Read the heatSourceDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
